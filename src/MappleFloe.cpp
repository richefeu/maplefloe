// ....

#include "MappleFloe.hpp"

MFloe::MFloe() {}

// ---------------------------------------------------------
// Print header
// ---------------------------------------------------------
void MFloe::head() {
  std::cout << '\n';
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⣶⣄⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣷⣴⣿⡄⠀⠀⠀⠀⠀⢀⡀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠰⣶⣾⣿⣿⣿⣿⣿⡇⠀⢠⣷⣤⣶⣿⡇⠀⠀⠀MappleFloe - " << MFLOE_VERSION << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠙⣿⣿⣿⣿⣿⣿⣿⣀⣿⣿⣿⣿⣿⣧⣀⠀⠀Université Grenoble Alpes" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⣷⣦⣀⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀" << std::endl;
  std::cout << "⠀⠀⢲⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠁⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠙⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠁⠀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠚⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠿⠿⠂⠀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠉⠙⢻⣿⣿⡿⠛⠉⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠋⠁⠀⠀⠀⠸⡄⠀⠀⠀⠀⠀⠀⠀⠀" << std::endl;
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⡀⠀⠀⠀     " << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------
// Saves the current configuration data to a file with
// the name 'conf<i>'.
// ---------------------------------------------------------
void MFloe::saveConf(int i) {
  char fname[256];
  sprintf(fname, "conf%d", i);
  saveConf(fname);
}

// ---------------------------------------------------------
// Saves the current configuration data to a file with
// the specified name.
// ---------------------------------------------------------
void MFloe::saveConf(const char *name) {
  std::ofstream conf(name);

  conf << "MFloe " << MFLOE_VERSION << std::endl; // format: progName version-date

  conf << std::setprecision(6);
  conf << "t " << t << std::endl;
  conf << "tmax " << tmax << std::endl;
  conf << "dt " << dt << std::endl;
  conf << "interClose " << interClose << std::endl;
  conf << "interOut " << interOut << std::endl;
  conf << "interHist " << interHist << std::endl;
  conf << "dVerlet " << dVerlet << std::endl;
  conf << "zgravNorm " << zgravNorm << std::endl;
  conf << "iconf " << iconf << std::endl;

  conf << "kn " << kn << std::endl;
  conf << "kt " << kt << std::endl;
  conf << "mu " << mu << std::endl;
  conf << "Gc " << Gc << std::endl;

  if (!Drivings.empty()) {
    conf << "nDriven " << Drivings.size() << std::endl;
    for (size_t d = 0; d < Drivings.size(); ++d) { Drivings[d]->write(conf); }
  }

  conf << "FloeElements " << FloeElements.size() << std::endl;
  conf << std::setprecision(15);
  for (size_t i = 0; i < FloeElements.size(); i++) {
    conf << FloeElements[i].pos << ' ' << FloeElements[i].zpos << ' ' << FloeElements[i].vel << ' '
         << FloeElements[i].zvel << ' ' << FloeElements[i].acc << ' ' << FloeElements[i].zacc << ' '
         << FloeElements[i].rot << ' ' << FloeElements[i].vrot << ' ' << FloeElements[i].arot << ' '
         << FloeElements[i].radius << ' ' << FloeElements[i].height << ' ' << FloeElements[i].inertia << ' '
         << FloeElements[i].mass << ' ';
    conf << std::endl;
  }

  size_t nbInteractions = 0;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1e-20 && Interactions[i].isBonded == false) { continue; }
    ++nbInteractions;
  }
  conf << "Interactions " << nbInteractions << std::endl;
  conf << std::setprecision(15);
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1e-20 && Interactions[i].isBonded == false) { continue; }
    conf << Interactions[i].i << ' ' << Interactions[i].j << ' ' << Interactions[i].isBonded << ' '
         << Interactions[i].fn << ' ' << Interactions[i].fnb << ' ' << Interactions[i].ft << ' ' << Interactions[i].ftb
         << ' ' << Interactions[i].fs << ' ' << Interactions[i].fsb << ' ' << Interactions[i].A << ' '
         << Interactions[i].coverage << ' ' << Interactions[i].dn0 << std::endl;
  }
}

// ---------------------------------------------------------
// A way to load 'conf<i>' files just by specifying i
// ---------------------------------------------------------
void MFloe::loadConf(int i) {
  char fname[256];
  sprintf(fname, "conf%d", i);
  loadConf(fname);
}

// ---------------------------------------------------------
// Loading a configuration file by giving its name path
// ---------------------------------------------------------
void MFloe::loadConf(const char *name) {
  std::ifstream conf(name);
  if (!conf.is_open()) { std::cout << MFLOE_WARN << "Cannot read " << name << std::endl; }

  // A warning for the placement of procesing commands
  auto warn_if_wrong_place = [&](const std::string &token) {
    if (FloeElements.empty()) {
      std::cout << MFLOE_WARN << "the command " << token << " should be placed after the definition of FloeElements"
                << std::endl;
    }
  };

  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "MFloe") { std::cout << MFLOE_WARN << "This seems not to be an input file for MappleFloe!" << std::endl; }
  std::string date;
  conf >> date;
  if (date != MFLOE_VERSION) { std::cout << MFLOE_WARN << "The version-date should be " << MFLOE_VERSION << std::endl; }

  std::string token;
  conf >> token;
  while (conf.good()) {
    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(conf, token); // ignore the rest of the current line
      conf >> token;        // next token
      continue;
    } else if (token == "t") {
      conf >> t;
    } else if (token == "tmax") {
      conf >> tmax;
    } else if (token == "dt") {
      conf >> dt;
    } else if (token == "interClose") {
      conf >> interClose;
    } else if (token == "interOut") {
      conf >> interOut;
    } else if (token == "interHist") {
      conf >> interHist;
    } else if (token == "dVerlet") {
      conf >> dVerlet;
    } else if (token == "zgravNorm") {
      conf >> zgravNorm;
    } else if (token == "iconf") {
      conf >> iconf;
    } else if (token == "kn") {
      conf >> kn;
    } else if (token == "kt") {
      conf >> kt;
    } else if (token == "mu") {
      conf >> mu;
    } else if (token == "Gc") {
      conf >> Gc;  
    } else if (token == "activationTime") {
      conf >> activationTime;  
    } else if (token == "healingTime") {
      conf >> healingTime;  
    } else if (token == "coverage0") {
      conf >> coverage0;  
    } 
    else if (token == "nDriven") {
      size_t nDriven;
      conf >> nDriven;
      Drivings.clear();
      std::string keyword;
      for (size_t d = 0; d < nDriven; ++d) {
        conf >> keyword;
        Driving *drv = Driving::create(keyword);
        if (drv != nullptr) {
          drv->read(conf);
          Drivings.push_back(drv);
        }
      }
    } else if (token == "FloeElements") {
      size_t nb;
      conf >> nb;
      FloeElements.clear();
      FloeElement P;
      for (size_t i = 0; i < nb; i++) {
        conf >> P.pos >> P.zpos >> P.vel >> P.zvel >> P.acc >> P.zacc >> P.rot >> P.vrot >> P.arot >> P.radius >>
            P.height >> P.inertia >> P.mass;

        FloeElements.push_back(P);
      }
    } else if (token == "Interactions") {
      size_t nb;
      conf >> nb;
      Interactions.clear();
      Interaction I;
      for (size_t k = 0; k < nb; k++) {
        conf >> I.i >> I.j >> I.isBonded >> I.fn >> I.fnb >> I.ft >> I.ftb >> I.fs >> I.fsb >> I.A >> I.coverage >>
            I.dn0;
        Interactions.push_back(I);
      }
    }

    // Processing commands that should be placed after
    // the definition of FloeElements (at the very end of the conf-file preferably).
    // This commands generally added in input-files but not saved in the conf-files
    else if (token == "computeMasseProperties") {
      warn_if_wrong_place(token);
      double density;
      conf >> density;
      computeMasseProperties(density);
    } else if (token == "activateBonds") {
      warn_if_wrong_place(token);
      double distanceMaxForGluing;
      conf >> distanceMaxForGluing;
      activateBonds(distanceMaxForGluing);
    }

    // Unknown token
    else {
      std::cout << MFLOE_WARN << "Unknown token: " << token << std::endl;
    }

    conf >> token;
  } // End of parsing loop

  // At this point, it is normally not possible to have no elements
  if (FloeElements.empty()) { std::cout << MFLOE_WARN << "No FloeElements after parsing" << std::endl; }

  // precompute things ========================================
  updateBoundLimits();
}

// ---------------------------------------------------------
// Compute the axis aligned bounding box in the sea-plane
// ---------------------------------------------------------
void MFloe::updateBoundLimits() {
  double xmin = FloeElements[0].pos.x - FloeElements[0].radius;
  double xmax = FloeElements[0].pos.x + FloeElements[0].radius;
  double ymin = FloeElements[0].pos.y - FloeElements[0].radius;
  double ymax = FloeElements[0].pos.y + FloeElements[0].radius;
  for (size_t i = 1; i < FloeElements.size(); ++i) {
    double x = FloeElements[i].pos.x + FloeElements[i].radius;
    if (x > xmax) { xmax = x; }
    x = FloeElements[i].pos.x - FloeElements[i].radius;
    if (x < xmin) { xmin = x; }

    double y = FloeElements[i].pos.y + FloeElements[i].radius;
    if (y > ymax) { ymax = y; }
    y = FloeElements[i].pos.y - FloeElements[i].radius;
    if (y < ymin) { ymin = y; }
  }

  aabb.min.set(xmin, ymin);
  aabb.max.set(xmax, ymax);
}

void MFloe::computeMasseProperties(double density) {
  for (size_t i = 0; i < FloeElements.size(); ++i) {
    double R = FloeElements[i].radius;
    double H = FloeElements[i].height;
    double V = M_PI * R * R * H;

    FloeElements[i].mass    = V * density;
    FloeElements[i].inertia = 0.5 * FloeElements[i].mass * R * R;
  }
}

// ---------------------------------------------------------
// Preprocessing function to activate the bonds
// ---------------------------------------------------------
void MFloe::activateBonds(double dmax) {
  // In case the conf-file has no interactions, the neighbor list is updated
  updateNeighbors(dmax);
  
  for (size_t e = 0 ; e < FloeElements.size(); e++) {
    FloeElements[e].Z = 0;
  }
  

  for (size_t k = 0; k < Interactions.size(); k++) {
    size_t i = Interactions[k].i;
    size_t j = Interactions[k].j;

    vec2r branch = FloeElements[j].pos - FloeElements[i].pos;

    double branchLen2 = norm2(branch);
    double sum        = dmax + FloeElements[i].radius + FloeElements[j].radius;
    if (branchLen2 <= sum * sum) {
      Interactions[k].isBonded = true;
      Interactions[k].dn0      = sqrt(branchLen2) - (FloeElements[i].radius + FloeElements[j].radius);      
      Interactions[k].coverage = 1.0;
      FloeElements[i].Z += 1;
      FloeElements[j].Z += 1;
    } // endif
  } // end loop over interactions
  
  
  // calcul des surfaces
  for (size_t k = 0; k < Interactions.size(); k++) {
    if (!Interactions[k].isBonded) {
      continue;
    }
    
    size_t i = Interactions[k].i;
    size_t j = Interactions[k].j;
    
    
    double h = std::min(FloeElements[i].height, FloeElements[j].height);
    
    double li = 2.0 * FloeElements[i].radius;
    if (FloeElements[i].Z > 2) {
      li *= sin(M_PI / FloeElements[i].Z);
    }
    double lj = 2.0 * FloeElements[j].radius;
    if (FloeElements[j].Z > 2) {
      lj *= sin(M_PI / FloeElements[j].Z);
    }
    double l = std::min(li, lj);
      
    Interactions[k].A = h * l;
    
  }
  
}

// ---------------------------------------------------------
// Periodically log useful informations
// on how the computation goes
// ---------------------------------------------------------
void MFloe::screenLog() {
  return;
  std::cout << std::endl;
  std::cout << "————————————————————————————————————————————————————————————————" << std::endl;
  std::cout << " iconf = " << iconf << "/" << iconfMaxEstimated << ", time = " << std::setprecision(10) << t
            << std::endl;
  /*
  std::cout << " Stress-particles:  " << std::endl;
  Sig.fancyPrint(mat4r::ColoredBrackets | mat4r::Scientific, 12);
  std::cout << " Stress-connection: " << std::endl;
  SigConnect.fancyPrint(mat4r::ColoredBrackets | mat4r::Scientific, 12);
  */
  // ...
  std::cout << "————————————————————————————————————————————————————————————————" << std::endl;
}

// ---------------------------------------------------------
// O(N^2) algorithm for updating the known neighbors by each
// FloeElement
// ---------------------------------------------------------
void MFloe::updateNeighbors(double dmax) {
  // store ft because the list will be erased before being rebuilt
  std::vector<Interaction> storedInteractions(Interactions.size());
  std::copy(Interactions.begin(), Interactions.end(), storedInteractions.begin());

  // now clear the list and rebuild it
  Interactions.clear();
  for (size_t i = 0; i < FloeElements.size(); ++i) {
    for (size_t j = i + 1; j < FloeElements.size(); ++j) {

      vec2r branch = FloeElements[j].pos - FloeElements[i].pos;

      double sum = dmax + FloeElements[i].radius + FloeElements[j].radius;
      if (norm2(branch) <= sum * sum) { Interactions.push_back(Interaction(i, j)); }
    }
  }

  // retrieve the embbeded values from Interactions that were present
  size_t k{0}, kold{0};
  for (; k < Interactions.size(); ++k) {
    while (kold < storedInteractions.size() && storedInteractions[kold].i < Interactions[k].i) { ++kold; }
    if (kold == storedInteractions.size()) { break; }

    while (kold < storedInteractions.size() && storedInteractions[kold].i == Interactions[k].i &&
           storedInteractions[kold].j < Interactions[k].j) {
      ++kold;
    }
    if (kold == storedInteractions.size()) { break; }

    if (storedInteractions[kold].i == Interactions[k].i && storedInteractions[kold].j == Interactions[k].j) {
      Interactions[k].copy(storedInteractions[kold]);
      ++kold;
    } else {
      /*combineParameters(k);*/
    }
  }

  // finish parameter-combinations in case of loop-break
  /*
  if (k < Interactions.size()) {
    if (k > 0) { --k; }
    for (; k < Interactions.size(); ++k) { combineParameters(k); }
  }
  */
}

// ---------------------------------------------------------
// The integration loop (velocity-Verlet scheme)
// ---------------------------------------------------------

/*
toutes les particules sont libre de mouvement (on ne peut pas bloquer imposer une vitesse).
Les statégies de chargement porterons idéalement sur l'imposition de force (éviter d'imposer une vitesse).
Par exemple, le vent pourrait être un champ de vitesse qu'on transforme en champs de forces en tenant compte de la
surface des particules.
*/

void MFloe::integrate() {
  double dt_2  = 0.5 * dt;
  double dt2_2 = 0.5 * dt * dt;

  iconfMaxEstimated = iconf + floor((tmax - t) / interHist);

  std::ofstream time_data_file("time_data.txt");

  size_t nDriven = Drivings.size();

  saveConf(iconf); // save before any computation
  screenLog();
  ++iconf;
  interHistC = 0.0;

  std::cout << MFLOE_INFO << "Beginning iterations." << std::endl;
  while (t < tmax) {

    for (size_t i = 0; i < nDriven; ++i) {
      Drivings[i]->set(t);
      FloeElements[i].vel.set(Drivings[i]->vx, Drivings[i]->vy);
      FloeElements[i].zvel = Drivings[i]->vz;
      FloeElements[i].vrot = Drivings[i]->vrot;

      FloeElements[i].pos += dt * FloeElements[i].vel;
      FloeElements[i].zpos += dt * FloeElements[i].zvel;
      FloeElements[i].rot += dt * FloeElements[i].vrot;
    }

    for (size_t i = nDriven; i < FloeElements.size(); i++) {
      FloeElements[i].pos += dt * FloeElements[i].vel + dt2_2 * FloeElements[i].acc;
      FloeElements[i].zpos += dt * FloeElements[i].zvel + dt2_2 * FloeElements[i].zacc;
      FloeElements[i].vel += dt_2 * FloeElements[i].acc;
      FloeElements[i].zvel += dt_2 * FloeElements[i].zacc;

      FloeElements[i].rot += dt * FloeElements[i].vrot + dt2_2 * FloeElements[i].arot;
      FloeElements[i].vrot += dt_2 * FloeElements[i].arot;
    }

    accelerations();

    for (size_t i = nDriven; i < FloeElements.size(); i++) {
      FloeElements[i].vel += dt_2 * FloeElements[i].acc;
      FloeElements[i].zvel += dt_2 * FloeElements[i].zacc;
      FloeElements[i].vrot += dt_2 * FloeElements[i].arot;
    }

    t += dt;

    interCloseC += dt;
    if (interCloseC > interClose - dt_2) {
      updateNeighbors(dVerlet);
      interCloseC = 0.0;
    }

    interOutC += dt;
    if (interOutC > interOut - dt_2) {

#if 1
      // ****** Version pour debug des modèles d'interaction ******
      if (!Interactions.empty()) {
        time_data_file << t << ' ' << Interactions[0].isBonded << ' ' << Interactions[0].fnb << ' '
                       << Interactions[0].ftb << ' ' << Interactions[0].fsb << ' ' << Interactions[0].fn << ' '
                       << Interactions[0].ft << ' ' << Interactions[0].fs << std::endl;
      }
#else
      time_data_file << t << ' ' << std::endl;
#endif
      interOutC = 0.0;
    }

    interHistC += dt;
    if (interHistC > interHist - dt_2) {
      saveConf(iconf);
      screenLog();
      ++iconf;
      interHistC = 0.0;
    }
  }

  return;
}

// ---------------------------------------------------------
// Compute the accelerations
// ---------------------------------------------------------
void MFloe::accelerations() {
  // Set forces and moments to zero
  for (size_t i = 0; i < FloeElements.size(); ++i) {
    FloeElements[i].force.reset();
    FloeElements[i].zforce = 0.0;
    FloeElements[i].moment = 0.0;
    FloeElements[i].acc.reset();
    FloeElements[i].zacc = 0.0; // gravity will be added at the end (or never, we need to discuss this point)
    FloeElements[i].arot = 0.0;
  }

  computeForcesAndMoments();

  // Finally, compute the accelerations (translation and rotation)
  for (size_t i = 0; i < FloeElements.size(); i++) {
    FloeElements[i].acc  = FloeElements[i].force / FloeElements[i].mass;
    FloeElements[i].zacc = FloeElements[i].zforce / FloeElements[i].mass /* - zgravNorm */;
    FloeElements[i].arot = FloeElements[i].moment / FloeElements[i].inertia;
  }
}

// ---------------------------------------------------------------
// The function to compute force interactions between FloeElements
// ---------------------------------------------------------------
void MFloe::computeForcesAndMoments() {

  for (size_t k = 0; k < Interactions.size(); ++k) {

    size_t i = Interactions[k].i;
    size_t j = Interactions[k].j;

    vec2r branch = FloeElements[j].pos - FloeElements[i].pos;

    vec2r unit_n = branch;
    double len   = unit_n.normalize();
    vec2r relVel = FloeElements[j].vel - FloeElements[i].vel; // Does not account for rotation (yet)

    vec2r unit_t(-unit_n.y, unit_n.x);
    double dn = len - FloeElements[i].radius - FloeElements[j].radius;
    // double vn = relVel * unit_n;

    double Li   = FloeElements[i].radius + 0.5 * dn;
    double Lj   = FloeElements[j].radius + 0.5 * dn;
    double vijt = relVel * unit_t - FloeElements[i].vrot * Li - FloeElements[j].vrot * Lj;
    double vijs = FloeElements[j].zvel - FloeElements[i].zvel;

    // ===================================
    // ICE-BONDED
    // ===================================
    if (Interactions[k].isBonded == true) { // i and j are bonded

      // calculate the elastic bonded forces
      Interactions[k].fnb = -kn * (dn - Interactions[k].dn0);
      Interactions[k].ftb = Interactions[k].ftb - kt * dt * vijt;
      Interactions[k].fsb = Interactions[k].fsb - kt * dt * vijs;

      //  ===================================
      //  BREAKAGE OF ICE-BONDS
      //  ===================================

      double fcn  = sqrt(2.0 * Interactions[k].coverage * Interactions[k].A * kn * Gc);
      double fct2 = 2.0 * Interactions[k].coverage * Interactions[k].A * kt * Gc;

      double crit = (-Interactions[k].fnb / fcn +
                     (Interactions[k].ftb * Interactions[k].ftb + Interactions[k].fsb * Interactions[k].fsb) / fct2) -
                    1.0;

      if (crit >= 0.0) {
        Interactions[k].isBonded = false;
        Interactions[k].t0 = -1.0;
        // cancel the bonding forces
        Interactions[k].fnb = 0.0;
        Interactions[k].ftb = 0.0;
        Interactions[k].fsb = 0.0;
        
        // cancel the contact force, that will be updated soon
        Interactions[k].fn = 0.0;
        Interactions[k].ft = 0.0; // integration resarts
        Interactions[k].fs = 0.0; // integration resarts
        
        // update the Z of each element
        FloeElements[i].Z -= 1;
        if (FloeElements[i].Z < 0) {
          std::cout << "Tabaarnak ! Z est devenu négatif" << std::endl;
          FloeElements[i].Z = 0;
        }
        FloeElements[j].Z -= 1;
        if (FloeElements[j].Z < 0) {
          std::cout << "Tabaarnak ! Z est devenu négatif" << std::endl;
          FloeElements[j].Z = 0;
        }
      }

    }
    // ===================================
    // NON-BONDED CONTACT
    // ===================================
    else if (dn < 0.0) { // it means that i and j are in contact (but not bonded)
      
      if (Interactions[k].t0 == -1.0) { // this means that the contact is 'new'
      std::cout << "Yop collage!!!! t = " << t << std::endl;
        Interactions[k].t0 = t;
        std::cout << " t0 = " << Interactions[k].t0 << std::endl;
      } 
      else if (Interactions[k].t0 > 0.0 && (t - Interactions[k].t0) > activationTime) {
        std::cout << "*************\n";
        std::cout << "t = " << t << std::endl;
        std::cout << "Interactions[k].t0 = " << Interactions[k].t0 << std::endl;
        std::cout << "activationTime = " << activationTime << std::endl;
        Interactions[k].isBonded = true;
        Interactions[k].coverage = 1.0; // pour tester, il faudra le mettre à coverage0
        Interactions[k].fnb = 0.0;
        Interactions[k].ftb = 0.0;
        Interactions[k].fsb = 0.0;
        Interactions[k].fn = 0.0;
        Interactions[k].ft = 0.0;
        Interactions[k].fs = 0.0;
        //Interactions[k].t0 = t; // pour le healing
        Interactions[k].dn0 = dn;
        continue; // next interaction k
      }

      // Elastic normal repulsion force
      Interactions[k].fn = -kn * dn;
      if (Interactions[k].fn < 0.0) { Interactions[k].fn = 0.0; }

      // Tangential force (friction)
      double ft    = Interactions[k].ft - kt * (dt * vijt);
      double limit = mu * Interactions[k].fn; // remember that fn > 0 because dn < 0
      if (fabs(ft) > limit) { ft = (ft > 0.0) ? limit : -limit; }
      Interactions[k].ft = ft;

      // out-of-sea-plane tangential force (friction)
      double fs = Interactions[k].fs - kt * (dt * vijs);
      if (fabs(fs) > limit) { fs = (fs > 0.0) ? limit : -limit; }
      Interactions[k].fs = fs;
    } else {
      // dn > 0 and without bond
      Interactions[k].t0 = -1.0;
    }

    // Resultant force and moment
    // (in sea-plane)
    vec2r f = (Interactions[k].fn + Interactions[k].fnb) * unit_n + (Interactions[k].ft + Interactions[k].ftb) * unit_t;
    FloeElements[i].force -= f;
    FloeElements[j].force += f;
    FloeElements[i].moment -= (Interactions[k].ft + Interactions[k].ftb) * Li;
    FloeElements[j].moment -= (Interactions[k].ft + Interactions[k].ftb) * Lj;
    // (out-of-sea-plane)
    FloeElements[i].zforce -= (Interactions[k].fs + Interactions[k].fsb);
    FloeElements[j].zforce += (Interactions[k].fs + Interactions[k].fsb);

  } // End loop over interactions
}
