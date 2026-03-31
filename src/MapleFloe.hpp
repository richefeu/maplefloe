// ...
#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Driving.hpp"
#include "FloeElement.hpp"
#include "Interaction.hpp"
#include "MFloeAssembler.hpp"
// #include "STF.hpp"

#include "toofus/AABB_2D.hpp"
#include "toofus/mat4.hpp"
#include "toofus/profiler.hpp"
#include "toofus/toofus-gate/termcolor/termcolor.hpp"
#include "toofus/vec2.hpp"

#define MFLOE_VERSION "2026.dev"
#define MFLOE_WARN termcolor::red << termcolor::bold << "Tabaarnack !" << termcolor::reset << ": "
#define MFLOE_INFO termcolor::cyan << "INFO" << termcolor::reset << ": "

#define BOLD_ termcolor::bold << termcolor::color<104>
#define NORMAL_ termcolor::reset

#ifdef DEBUG
#define MFLOE_SHOW(V)                                                                                        \
  if (debugEnabled)                                                                                          \
  std::cout << termcolor::color<155> << #V << "" << termcolor::reset << " = " << V << std::flush << std::endl
#else
#define MFLOE_SHOW(V)
#endif

#define YIELD_PARABOLA 0
#define YIELD_HEAVISIDE 1
#define YIELD_ELLIPSE 2
#define YIELD_ELLIPSE_ASYM 3

class MFloe {
public:
  std::vector<FloeElement> FloeElements;
  std::vector<Interaction> Interactions;

  std::vector<Driving *> Drivings;

  // Spatio Temporal Fields
  // STF<double> *waveFieldModel{nullptr};
  // STF<vec2r> *windFieldModel{nullptr};

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

  int healingEnabled{1};
  double activationTime{0.0}; // required contact duration for changing to a healing bonded
  double healingTime{0.0};    // reference duration that tune the healing rate
  double coverage0{0.0};      // healingRatio or healingProgress

  int NbBondsInit{0}; // a reference initial number of bonds (computed when bonds are activated)

  std::ofstream breakage_data_file;

  AABB_2D aabb;

  bool verbose{false};
  bool debugEnabled{true};

  int yieldSurfaceModel{YIELD_PARABOLA};

  // interaction shared parameters
  double kn{1e6};           // contact/bond normal stiffness
  double kt{1e6};           // contact/bond tangent stiffness
  double mu{0.5};           // contact friction coefficient
  double Gc{2.0};           // critical energy release rate (per unit crack area)
  double GcComprRatio{1.0}; // Gc_compression / Gc_traction

  // default parameters used for pre-processing
  double defaultFloeHeight{1.0};
  double defaultFloeMassDensity{1027.0};

  // viscous drag forces
  double DragCoef{0.5};
  double seaMassDensity{1000.0};

  MFloe();
  void head();

  // Core functions for the computations
  void initialChecks();
  void integrate();
  void accelerations();
  void computeForcesAndMoments();
  void addDragForces();
  void updateNeighbors(double dmax);

  // Save and Load the configuration files (conf-files)
  void saveConf(int i);
  void saveConf(const char *name);
  void loadConf(int i);
  void loadConf(const char *name);

  // Functions for updating relevant details
  void updateBoundLimits();

  int NbBondsCurrent{0};
  int NbHealingCurrent{0};
  double HealingArea{0.0};
  double HealingAreaCovered{0.0};
  double Ktrans{0.0};
  double Krot{0.0};
  mat4r VStress; // this is the moment tensor (V * Stress)
  void computeMechanicalIndicators();
  void screenLog();

  // pre-processing functions
  void computeMassProperties(double density);
  void activateBonds(double dmax);
};
