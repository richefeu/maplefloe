#include "Loading.hpp"
#include "SPICE.hpp"

// Remark: if the creation fails, nullptr is returned
Loading *Loading::create(const std::string &token) {
  if (token == "ShearPV") {
    return new ShearPV();
  } else if (token == "ShearVV") {
    return new ShearVV();
  }

  return nullptr;
}

Loading::~Loading() {}

// ========== ShearVV ==========

ShearVV::ShearVV() {}

void ShearVV::read(std::istream &is) {
  is >> vx >> vy;
}

void ShearVV::write(std::ostream &os) {
  os << "ShearVV " << vx << ' ' << vy << std::endl;
}

void ShearVV::velocityVerlet_halfStep1() {
  for (size_t m = 0; m < box->top.pos.size(); ++m) {
    box->top.pos[m].x += box->dt * vx;
    box->top.pos[m].y += box->dt * vy;

    if (box->top.pos[m].x > box->xmax) { box->top.pos[m].x -= box->xmax - box->xmin; }
    if (box->top.pos[m].x < box->xmin) { box->top.pos[m].x += box->xmax - box->xmin; }
  }
}

// ========== ShearPV ==========

ShearPV::ShearPV() {}

void ShearPV::read(std::istream &is) {
  is >> pressure >> velocity;
}

void ShearPV::write(std::ostream &os) {
  os << "ShearPV " << pressure << ' ' << velocity << std::endl;
}

void ShearPV::init() {
  top_mass = 0.0;
  for (size_t m = 0; m < box->top.Idx.size(); ++m) {
    size_t idx = box->bottom.Idx[m];
    top_mass += box->Particles[idx].mass;
  }
  top_accy = 0.0;
}

void ShearPV::velocityVerlet_halfStep1() {
  double dt2_2 = 0.5 * box->dt * box->dt;
  double dt_2  = 0.5 * box->dt;
  for (size_t m = 0; m < box->top.pos.size(); ++m) {
    box->top.pos[m].x += box->dt * velocity;
    box->top.pos[m].y += box->dt * top_vy + dt2_2 * top_accy;

    if (box->top.pos[m].x > box->xmax) { box->top.pos[m].x -= box->xmax - box->xmin; }
    if (box->top.pos[m].x < box->xmin) { box->top.pos[m].x += box->xmax - box->xmin; }

    top_vy += dt_2 * top_accy;
  }
}

void ShearPV::velocityVerlet_halfStep2() {
  double dt_2 = 0.5 * box->dt;
  for (size_t m = 0; m < box->top.pos.size(); ++m) { top_vy += dt_2 * top_accy; }
}

void ShearPV::forceDrivenAcceleration() {
  // il faut faire des sommes sur les croix pour calculer top_fy
  // mais pour le moment ce n'est pas stocker...
  top_accy = -pressure * (box->xmax - box->xmin) / top_mass;
}
