// ....
#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Interaction.hpp"
#include "FloeElement.hpp"
#include "Driving.hpp"

#include "toofus/vec2.hpp"
#include "toofus/mat4.hpp"
#include "toofus/AABB_2D.hpp"

#define MFLOE_VERSION "2026.dev"
#define MFLOE_WARN "\033[0m\033[31m\033[1m\033[4mWARN\033[24m\033[39m\033[0m: "
#define MFLOE_INFO "\033[0m\033[32m\033[1m\033[4mINFO\033[24m\033[39m\033[0m: "

class MFloe {
public:
  std::vector<FloeElement> FloeElements;
  std::vector<Interaction> Interactions;
  
  std::vector<Driving*> Drivings;

  // parameters
  double t{0.0};
  double tmax{5.0};
  double dt{1e-6};

  double interLookC{0.0};
  double interCloseC{0.0}, interClose{0.01}, dVerlet{0.0};
  double interOutC{0.0}, interOut{0.1};
  double interHistC{0.0}, interHist{0.25};

  int iconf{0};
  int iconfMaxEstimated{0};
  double zgravNorm{9.81};
   
  double activationTime{0.0};
  double healingTime{0.0};
  double coverage0{0.0};

  AABB_2D aabb;

  bool verbose{false};
  
  // interaction shared parameters
  double kn{1e6}; // contact/bond normal stiffness
  double kt{1e6}; // contact/bond tangent stiffness
  double mu{0.5}; // contact friction coefficient
  double Gc{2.0}; // ...
  
  MFloe();
  void head();

  // Core functions for the computations
  void integrate();
  void accelerations();
  void computeForcesAndMoments();
  void updateNeighbors(double dmax);

  // Save and Load the configuration files (conf-files)
  void saveConf(int i);
  void saveConf(const char *name);
  void loadConf(int i);
  void loadConf(const char *name);

  // Functions for updating relevant details
  void updateBoundLimits();

  void screenLog();

  // pre-processing functions
  void computeMasseProperties(double density);
  void activateBonds(double dmax);
};
