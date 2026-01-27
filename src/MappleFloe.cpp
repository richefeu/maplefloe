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
  std::cout << "⠀⠀⠀⠀⠀⠀⠀⠙⣿⣿⣿⣿⣿⣿⣿⣀⣿⣿⣿⣿⣿⣧⣀⠀⠀" << std::endl;
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
  // conf << "gravity " << gravity << std::endl;
  conf << "iconf " << iconf << std::endl;
  // conf << "Xperiod " << xmin << ' ' << xmax << std::endl;

  // conf << "bottom " << bottom.Idx.size() << ' ' << bottom.K << std::endl;
  // conf << std::setprecision(15);
  // for (size_t i = 0; i < bottom.Idx.size(); i++) { conf << bottom.Idx[i] << ' ' << bottom.pos[i] << std::endl; }

  // conf << std::setprecision(6);
  // conf << "top " << top.Idx.size() << ' ' << top.K << std::endl;
  // conf << std::setprecision(15);
  // for (size_t i = 0; i < top.Idx.size(); i++) { conf << top.Idx[i] << ' ' << top.pos[i] << std::endl; }

  // conf << "Loading ";
  // if (Load != nullptr) {
  //   Load->write(conf);
  // } else {
  //   std::cout << SPICE_WARN << "Loading is not defined" << std::endl;
  //   conf << "!! TO BE ADDED !!" << std::endl;
  // }

  // LA il faut faire le ménage

  conf << "FloeElements " << FloeElements.size() << std::endl;
  conf << std::setprecision(15);
  for (size_t i = 0; i < FloeElements.size(); i++) {
    conf << FloeElements[i].pos << ' ' << FloeElements[i].vel << ' ' << FloeElements[i].acc << ' ' << FloeElements[i].rot << ' '
         << FloeElements[i].vrot << ' ' << FloeElements[i].arot << ' ' << FloeElements[i].radius << ' ' << FloeElements[i].inertia
         << ' ' << FloeElements[i].mass << ' ';
    //conf << FloeElements[i].normalStiffness << ' ' << FloeElements[i].tangentialStiffness << ' '
    //     << FloeElements[i].normalViscDampingRate << ' ' << FloeElements[i].friction << ' ' << FloeElements[i].rollingFriction
    //     << ' ' << FloeElements[i].adhesion << ' ' << FloeElements[i].GcGlue << std::endl;
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
         /*<< Interactions[i].isSameMaterialBond */ << ' ' << Interactions[i].fn << ' ' << Interactions[i].fnb << ' '
         << Interactions[i].ft << ' ' << Interactions[i].ftb << ' ' << Interactions[i].damp << ' ' << Interactions[i].Gs
         << ' ' << Interactions[i].dn0 << std::endl;
  }
}

// ---------------------------------------------------------
// A fast way to load 'conf<i>' files (just by specifying i)
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
  if (prog != "MappleFloe") { std::cout << MFLOE_WARN << "This is not file for MappleFloe executable!" << std::endl; }
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
    } else if (token == "iconf") {
      conf >> iconf;
    } /*else if (token == "Loading") {
      std::string loadingName;
      conf >> loadingName;
      if (Load != nullptr) {  // Load already exist
        getline(conf, token); // ignore the rest of the current line
        conf >> token;        // next token
        continue;
      }
      if (loadingName[0] == '!') {
        std::cout << SPICE_WARN << "Remember to define the Loading" << std::endl;
        getline(conf, token); // ignore the rest of the current line
        conf >> token;        // next token
        continue;
      }
      Load = Loading::create(loadingName.c_str());
      if (Load == nullptr) {
        std::cout << SPICE_WARN << "Could not create the Loading" << std::endl;
      } else {
        Load->box = this;
        Load->read(conf);
      }

    } */
    else if (token == "FloeElements") {
      size_t nb;
      conf >> nb;
      FloeElements.clear();
      FloeElement P;
      for (size_t i = 0; i < nb; i++) {
        conf >> P.pos >> P.vel >> P.acc >> P.rot >> P.vrot >> P.arot >> P.radius >> P.inertia >> P.mass;
        //conf >> P.normalStiffness >> P.tangentialStiffness >> P.normalViscDampingRate >> P.friction >>
        //    P.rollingFriction >> P.adhesion >> P.GcGlue;

        FloeElements.push_back(P);
      }
    } else if (token == "Interactions") {
      size_t nb;
      conf >> nb;
      Interactions.clear();
      Interaction I;
      for (size_t k = 0; k < nb; k++) {
        conf >> I.i >> I.j >> I.isBonded >> /*I.isSameMaterialBond >>*/ I.fn >> I.fnb >> I.ft >> I.ftb >> I.damp >> I.Gs >>
            I.dn0;
        Interactions.push_back(I);
      }
    }

    // Processing commands that should be placed after
    // the definition of FloeElements (at the very end of the conf-file preferably).
    // This commands generally added in input-files but not saved in the conf-files
    else if (token == "activateBonds") {
      warn_if_wrong_place(token);
      //bool isSameMaterial;
      double distanceMaxForGluing;
      conf /*>> isSameMaterial*/ >> distanceMaxForGluing;
      activateBonds(/*isSameMaterial,*/ distanceMaxForGluing);
    }

    // Unknown token
    else {
      std::cout << MFLOE_WARN << "Unknown token: " << token << std::endl;
    }

    conf >> token;
  }

  // precompute things ========================================
  updateBoundLimits();
  // applyEmbeddedDataProfiles(); // set the property profiles that have been eventually (re-)set
  // combineParameters();

  // some checks
  if (FloeElements.empty()) { std::cout << MFLOE_WARN << "No FloeElements" << std::endl; }
  // if (Load == nullptr) { std::cout << MFLOE_WARN << "No Loading defined!" << std::endl; }
}

// ---------------------------------------------------------
// Compute ymin and ymax
// ---------------------------------------------------------
void MFloe::updateBoundLimits() {
  // TODO: Il faut calculer les xmin et xmax aussi
  // et mettre cela dans un AABB
  ymin = FloeElements[0].pos.y - FloeElements[0].radius;
  ymax = FloeElements[0].pos.y + FloeElements[0].radius;
  for (size_t i = 1; i < FloeElements.size(); ++i) {
    double y = FloeElements[i].pos.y + FloeElements[i].radius;
    if (y > ymax) { ymax = y; }
    y = FloeElements[i].pos.y - FloeElements[i].radius;
    if (y < ymin) { ymin = y; }
  }
}

// ---------------------------------------------------------
// Preprocessing function to activate the bonds
// ---------------------------------------------------------
void MFloe::activateBonds(/*bool sameMaterial,*/ double dmax) {
  // TODO: FONCTION A ADAPTER

  // In case the conf-file has no interactions, the neighbor list is updated
  updateNeighbors(dmax);

  //double Lperiod = xmax - xmin;

  std::cout << " routine to activate bonds" << std::endl;

  for (size_t k = 0; k < Interactions.size(); k++) {
    size_t i = Interactions[k].i;
    size_t j = Interactions[k].j;

    vec2r branch = FloeElements[j].pos - FloeElements[i].pos;
    //branch.x += getBranchShift(branch.x, Lperiod);

    double branchLen2 = norm2(branch);
    double sum        = dmax + FloeElements[i].radius + FloeElements[j].radius;
    if (branchLen2 <= sum * sum) {

      // switch to a cemented/bonded link
      Interactions[k].isBonded           = true;
      //Interactions[k].isSameMaterialBond = sameMaterial;
      // TODO Use Gc and min diameter to define a threshold Wmax

      std::cout << " activate bond = " << k << std::endl;

      // double dn = sqrt(branchLen2) - (FloeElements[i].radius + FloeElements[j].radius);
      // if (dn >= 0.0) Interactions[k].dn0 = dn;
      // else Interactions[k].dn0 = 0.0;
      Interactions[k].dn0 = sqrt(branchLen2) - (FloeElements[i].radius + FloeElements[j].radius);

    } // endif
  } // end loop over interactions
}

// ---------------------------------------------------------
// Periodically log useful informations
// on how the computation goes
// ---------------------------------------------------------
void MFloe::screenLog() {
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
  // TODO: adapter

  // store ft because the list will be erased before being rebuilt
  std::vector<Interaction> storedInteractions(Interactions.size());
  std::copy(Interactions.begin(), Interactions.end(), storedInteractions.begin());

  // now clear the list and rebuild it
  Interactions.clear();
  // double Lperiod = xmax - xmin;
  for (size_t i = 0; i < FloeElements.size(); ++i) {
    for (size_t j = i + 1; j < FloeElements.size(); ++j) {

      vec2r branch = FloeElements[j].pos - FloeElements[i].pos;
      // branch.x += getBranchShift(branch.x, Lperiod);

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
Par exemple, le vent pourrait être un champ de vitesse qu'on transforme en champs de forces en tenant compte de la surface des particules.
*/


void MFloe::integrate() {
  double dt_2  = 0.5 * dt;
  double dt2_2 = 0.5 * dt * dt;

  iconfMaxEstimated = iconf + floor((tmax - t) / interHist);

  std::ofstream stressOut("stress.out.txt");

  // double Lperiod = xmax - xmin;

  // Load->init();

  saveConf(iconf); // save before any computation
  screenLog();
  ++iconf;
  interHistC = 0.0;

  std::cout << MFLOE_INFO << "Beginning iterations." << std::endl;
  while (t < tmax) {

    // Load->servo();

    // Load->velocityVerlet_halfStep1();
    for (size_t i = 0; i < FloeElements.size(); i++) {
      FloeElements[i].pos += dt * FloeElements[i].vel + dt2_2 * FloeElements[i].acc;
      FloeElements[i].vel += dt_2 * FloeElements[i].acc;

      FloeElements[i].rot += dt * FloeElements[i].vrot + dt2_2 * FloeElements[i].arot;
      FloeElements[i].vrot += dt_2 * FloeElements[i].arot;
    }

    accelerations();

    // Load->velocityVerlet_halfStep2();
    for (size_t i = 0; i < FloeElements.size(); i++) {
      FloeElements[i].vel += dt_2 * FloeElements[i].acc;
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
      stressOut << t << " " << Sig /*<< " " << SigConnect*/ << std::endl;
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
    FloeElements[i].moment = 0.0;
    FloeElements[i].acc    = gravity;
    FloeElements[i].arot   = 0.0;
  }
  Sig.reset();
  //SigConnect.reset();

  computeForcesAndMoments();
  // computeFarConnectionForces();

  // updateBoundLimits();
  // double invV = 1.0 / ((xmax - xmin) * (ymax - ymin));
  // Sig *= invV;
  // SigConnect *= invV;

  // Finally compute the accelerations (translation and rotation)
  // Load->forceDrivenAcceleration();
  for (size_t i = 0; i < FloeElements.size(); i++) {
    FloeElements[i].acc  = gravity + FloeElements[i].force / FloeElements[i].mass;
    FloeElements[i].arot = FloeElements[i].moment / FloeElements[i].inertia;
  }
}

// ---------------------------------------------------------
// The function to compute force interactions between FloeElements
// ---------------------------------------------------------
void MFloe::computeForcesAndMoments() {
  // double Lperiod = xmax - xmin;
  for (size_t k = 0; k < Interactions.size(); ++k) {

    size_t i = Interactions[k].i;
    size_t j = Interactions[k].j;

    vec2r branch = FloeElements[j].pos - FloeElements[i].pos;
    // branch.x += getBranchShift(branch.x, Lperiod);

    vec2r unit_n = branch;
    double len   = unit_n.normalize();
    vec2r relVel = FloeElements[j].vel - FloeElements[i].vel; // Does not account for rotation yet

    vec2r unit_t(-unit_n.y, unit_n.x);
    double dn = len - FloeElements[i].radius - FloeElements[j].radius;
    double vn = relVel * unit_n;

    double Li   = FloeElements[i].radius + 0.5 * dn;
    double Lj   = FloeElements[j].radius + 0.5 * dn;
    double vijt = relVel * unit_t - FloeElements[i].vrot * Li - FloeElements[j].vrot * Lj;

    // ===================================
    // BOND
    // ===================================
    if (Interactions[k].isBonded == true) { // i and j are bonded

      // FIXME: il faut ajouter les variables dnb et dtb

      // calculate the bonded forces
      Interactions[k].fnb = -Interactions[k].kn * (dn - Interactions[k].dn0);
      // Interactions[k].ftb = Interactions[k].ft - Interactions[k].kt * dt * vijt;

    } else { // i and j not bonded

      Interactions[k].fnb = 0;
      Interactions[k].ftb = 0;
    }

    // ===================================
    // CONTACT
    // ===================================
    // TODO : IL FAUT RETIRER LA NOTION DE isSameMaterialBond
    // frictional contact (possibly in addition to bond)

      if (Interactions[k].isBonded == false && dn < 0.0) { // it means that i and j are in contact

        double fne = -Interactions[k].kn * dn;   // elastic normal force
        double fnv = -Interactions[k].damp * vn; // viscous normal force

        Interactions[k].fn = fne + fnv;
        if (Interactions[k].fn < 0.0) { Interactions[k].fn = 0.0; }

        // Tangential force (friction)
        double ft    = Interactions[k].ft - Interactions[k].kt * (dt * vijt);
        double ftest = Interactions[k].mu * Interactions[k].fn; // remember that fn >= 0
        if (fabs(ft) > ftest) { ft = (ft > 0.0) ? ftest : -ftest; }
        Interactions[k].ft = ft;

        // Adhesion
        // TODO:
        // Interactions[k].fn -= fadh;
      }
    

    // ===================================
    // BREAKAGE OF BONDS
    // ===================================
    // TODO if (isBonded) ...

    // Resultant force and moment
    vec2r f = (Interactions[k].fn + Interactions[k].fnb) * unit_n + (Interactions[k].ft + Interactions[k].ftb) * unit_t;
    FloeElements[i].force -= f;
    FloeElements[j].force += f;
    FloeElements[i].moment -= (Interactions[k].ft + Interactions[k].ftb) * Li;
    FloeElements[j].moment -= (Interactions[k].ft + Interactions[k].ftb) * Lj;

    // Internal stress
    Sig.xx += f.x * branch.x;
    Sig.xy += f.x * branch.y;
    Sig.yx += f.y * branch.x;
    Sig.yy += f.y * branch.y;

  } // Loop over interactions
}
