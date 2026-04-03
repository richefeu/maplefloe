#include "BiaxSystem.hpp"
#include "MapleFloe.hpp"

BiaxSystem::BiaxSystem() {}

void BiaxSystem::init(vec2r &min, vec2r &max, int nbCuts) {
  box.min = min;
  box.max = max;

  nbSlices       = nbCuts;
  t_ypos         = max.y;
  b_ypos         = min.y;
  segmentLengths = (t_ypos - b_ypos) / (double)(nbSlices);

  t_yvel = b_yvel = 0.0;
  t_yacc = b_yacc = 0.0;
  l_xpos.resize(nbSlices, min.x);
  r_xpos.resize(nbSlices, max.x);
  l_xvel.resize(nbSlices, 0.0);
  r_xvel.resize(nbSlices, 0.0);
  l_xacc.resize(nbSlices, 0.0);
  r_xacc.resize(nbSlices, 0.0);
  l_res_force.resize(nbSlices);
  r_res_force.resize(nbSlices);
}

void BiaxSystem::read(std::istream &is) {
  vec2r min, max;
  int nbCuts;
  is >> min >> max >> nbCuts;
  init(min, max, nbCuts);
  is >> mass;
  is >> kn_v >> kt_v >> dampRate_v;
  is >> kn_h >> kt_h >> dampRate_h;

  // std::cout << "mass = " << mass << std::endl;
  // std::cout << "kn_v = " << kn_v << ", kn_v = " << kn_v << ", dampRate_v = " << dampRate_v <<  std::endl;
  // std::cout << "kn_h = " << kn_v << ", kn_h = " << kn_h << ", dampRate_h = " << dampRate_h <<  std::endl;
}

void BiaxSystem::readMode(std::istream &is) {
  std::string token;
  is >> token;
  if (token == "ISO") {
    mode = BIAX_MODE_ISO;
    is >> pressure;
    velocity = 0.0;
  } else if (token == "VP") {
    mode = BIAX_MODE_VP;
    is >> velocity >> pressure;
  } else {
    std::cout << MFLOE_WARN << "Unknown biax mode: " << token << std::endl;
  }
}

void BiaxSystem::readData(std::istream &is) {
  is >> t_ypos >> b_ypos;
  for (size_t s = 0; s < l_xpos.size(); s++) { is >> l_xpos[s]; }
  for (size_t s = 0; s < r_xpos.size(); s++) { is >> r_xpos[s]; }

  is >> t_yvel >> b_yvel;
  for (size_t s = 0; s < l_xvel.size(); s++) { is >> l_xvel[s]; }
  for (size_t s = 0; s < r_xvel.size(); s++) { is >> r_xvel[s]; }

  is >> t_res_force >> b_res_force;
  for (size_t s = 0; s < l_res_force.size(); s++) { is >> l_res_force[s]; }
  for (size_t s = 0; s < r_res_force.size(); s++) { is >> r_res_force[s]; }
}

void BiaxSystem::write(std::ostream &os) {

  for (size_t s = 0; s < l_xpos.size(); s++) {
    if (box.min.x > l_xpos[s]) box.min.x = l_xpos[s];
  }
  for (size_t s = 0; s < r_xpos.size(); s++) {
    if (box.max.x < r_xpos[s]) box.max.x = r_xpos[s];
  }
  box.min.y = b_ypos;
  box.max.y = t_ypos;

  os << "BiaxSystem " << box.min << ' ' << box.max << ' ' << nbSlices << std::endl;
  os << ' ' << mass << std::endl;
  os << ' ' << kn_v << ' ' << kt_v << ' ' << dampRate_v << std::endl;
  os << ' ' << kn_h << ' ' << kt_h << ' ' << dampRate_h << std::endl;

  if (mode == BIAX_MODE_ISO) {
    os << "BiaxMode ISO" << ' ' << pressure << std::endl;
  } else if (mode == BIAX_MODE_VP) {
    os << "BiaxMode VP" << ' ' << velocity << ' ' << pressure << std::endl;
  }

  os << "BiaxData" << std::endl;

  os << t_ypos << ' ' << b_ypos;
  for (size_t s = 0; s < l_xpos.size(); s++) { os << ' ' << l_xpos[s]; }
  for (size_t s = 0; s < r_xpos.size(); s++) { os << ' ' << r_xpos[s]; }
  os << std::endl;

  os << t_yvel << ' ' << b_yvel;
  for (size_t s = 0; s < l_xvel.size(); s++) { os << ' ' << l_xvel[s]; }
  for (size_t s = 0; s < r_xvel.size(); s++) { os << ' ' << r_xvel[s]; }
  os << std::endl;

  os << t_res_force << ' ' << b_res_force;
  for (size_t s = 0; s < l_res_force.size(); s++) { os << ' ' << l_res_force[s]; }
  for (size_t s = 0; s < r_res_force.size(); s++) { os << ' ' << r_res_force[s]; }
  os << std::endl;
}

void BiaxSystem::resetForces() {
  t_res_force.reset();
  b_res_force.reset();
  t_yacc = 0.0;
  b_yacc = 0.0;
  for (size_t s = 0; s < l_res_force.size(); s++) {
    l_res_force[s].reset();
    l_xacc[s] = 0.0;
  }
  for (size_t s = 0; s < r_res_force.size(); s++) {
    r_res_force[s].reset();
    r_xacc[s] = 0.0;
  }
}

// première version sans viscosité et sans frottement
void BiaxSystem::computeContactForces(MFloe &floe) {

  double dn{0.0};
  double H = t_ypos - b_ypos;

  size_t nDriven = floe.Drivings.size();
  for (size_t i = nDriven; i < floe.FloeElements.size(); i++) {
    double R  = floe.FloeElements[i].radius;
    vec2r pos = floe.FloeElements[i].pos;
    // vec2r vel = floe.FloeElements[i].vel;

    dn = (t_ypos - pos.y) - R; // top
    if (dn < 0.0) {
      double fy = -kn_v * dn;
      t_res_force.y += fy;
      floe.FloeElements[i].force.y -= fy;
    } else {
      dn = (pos.y - b_ypos) - R; // bottom
      if (dn < 0.0) {
        double fy = kn_v * dn;
        b_res_force.y += fy;
        floe.FloeElements[i].force.y -= fy;
      }
    }

    int s = static_cast<int>(floor(nbSlices * (pos.y - b_ypos) / H));
    if (s < 0 || s >= nbSlices) continue;

    dn = (pos.x - l_xpos[s]) - R; // left
    if (dn < 0.0) {
      double fx = kn_h * dn;
      l_res_force[s].x += fx;
      floe.FloeElements[i].force.x -= fx;
    } else {
      dn = (r_xpos[s] - pos.x) - R; // right
      if (dn < 0.0) {
        double fx = -kn_h * dn;
        r_res_force[s].x += fx;
        floe.FloeElements[i].force.x -= fx;
      }
    }
  }
}

void BiaxSystem::computeAccelerations() {
  if (mode == BIAX_MODE_ISO) {
    double ltop = r_xpos[nbSlices - 1] - l_xpos[nbSlices - 1];
    double lbot = r_xpos[0] - l_xpos[0];

    t_yacc = (t_res_force.y - pressure * ltop) / mass;
    b_yacc = (b_res_force.y + pressure * lbot) / mass;
  } else {
    t_yacc = 0.0;
    b_yacc = 0.0;
  }

  for (size_t s = 0; s < l_res_force.size(); s++) {
    l_xacc[s] = (l_res_force[s].x + pressure * segmentLengths) / mass;
  }
  for (size_t s = 0; s < r_res_force.size(); s++) {
    r_xacc[s] = (r_res_force[s].x - pressure * segmentLengths) / mass;
  }
}

void BiaxSystem::drivePosVel1(double dt, double dt_2, double dt2_2) {

  if (mode == BIAX_MODE_VP) {
    t_yvel = -0.5 * velocity;
    b_yvel = 0.5 * velocity;
    t_yacc = 0.0;
    b_yacc = 0.0;
  }

  // top
  t_ypos += dt * t_yvel + dt2_2 * t_yacc;
  t_yvel += dt_2 * t_yacc;

  // bottom
  b_ypos += dt * b_yvel + dt2_2 * b_yacc;
  b_yvel += dt_2 * b_yacc;

  segmentLengths = (t_ypos - b_ypos) / (double)nbSlices;

  // left
  for (size_t i = 0; i < l_xpos.size(); i++) {
    l_xpos[i] += dt * l_xvel[i] + dt2_2 * l_xacc[i];
    l_xvel[i] += dt_2 * l_xacc[i];
  }

  // right
  for (size_t i = 0; i < r_xpos.size(); i++) {
    r_xpos[i] += dt * r_xvel[i] + dt2_2 * r_xacc[i];
    r_xvel[i] += dt_2 * r_xacc[i];
  }
}

void BiaxSystem::driveVel2(double dt_2) {

  if (mode == BIAX_MODE_ISO) {
    t_yvel += dt_2 * t_yacc; // top
    b_yvel += dt_2 * b_yacc; // bottom
  }

  // left
  for (size_t i = 0; i < l_xpos.size(); i++) { l_xvel[i] += dt_2 * l_xacc[i]; }

  // right
  for (size_t i = 0; i < r_xpos.size(); i++) { r_xvel[i] += dt_2 * r_xacc[i]; }
}
