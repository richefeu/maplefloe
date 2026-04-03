#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class MFloe;

#include "toofus/AABB_2D.hpp"
#include "toofus/vec2.hpp"

#define BIAX_MODE_ISO 0
#define BIAX_MODE_VP 1

class BiaxSystem {
public:
  AABB_2D box;
  int nbSlices{1};
  double segmentLengths{0.0};
  double t_ypos;
  double b_ypos;
  double t_yvel;
  double b_yvel;
  double t_yacc;
  double b_yacc;
  std::vector<double> l_xpos;
  std::vector<double> r_xpos;
  std::vector<double> l_xvel;
  std::vector<double> r_xvel;
  std::vector<double> l_xacc;
  std::vector<double> r_xacc;

  int mode{BIAX_MODE_ISO};
  double pressure{0.0};
  double velocity{0.0};

  vec2r t_res_force;
  vec2r b_res_force;
  std::vector<vec2r> l_res_force;
  std::vector<vec2r> r_res_force;

  double mass{1.0};
  double kn_v{1000.0}, kt_v{1000.0}, dampRate_v{0.0};
  double kn_h{1000.0}, kt_h{1000.0}, dampRate_h{0.0};

  BiaxSystem();
  void init(vec2r &min, vec2r &max, int nbCuts);
  void read(std::istream &is);
  void readMode(std::istream &is);
  void readData(std::istream &is);
  void write(std::ostream &os);
  
  void resetForces();
  void computeContactForces(MFloe &floe);
  void computeAccelerations();
  void drivePosVel1(double dt, double dt_2, double dt2_2);
  void driveVel2(double dt_2);
};