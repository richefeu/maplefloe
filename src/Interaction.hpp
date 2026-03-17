#pragma once

#include <cmath>
#include <cstddef>
#include <functional>

#include "FloeElement.hpp"
class MFloe;

struct Interaction {
  size_t i{0};
  size_t j{0};

  // flags
  bool isBonded{false}; // Indicates if the particle has any bonds or connections to other particles

  // forces
  double fn{0.0};  // contact normal force
  double fnb{0.0}; // bond normal force
  double ft{0.0};  // contact tangential force (in plane)
  double ftb{0.0}; // bond tangential force (in plane)
  double fs{0.0};  // contact tangential force (out-of-plane)
  double fsb{0.0}; // bond tangential force (out-of-plane)

  // double meff{0.0};     // effective mass
  double A{0.0};        // bonding surface area
  double coverage{0.0}; // Ice coverage ratio
  double dn0{0.0};      // initial distance for bonds; fnb(dn0) = 0

  double t0{-1.0}; // il faudra un bond commentaire pour expliquer

  Interaction();
  Interaction(size_t I, size_t J);

  void copy(Interaction &I);
  void computeBondedArea(double h_i, double radius_i, int Nbonds_i, double h_j, double radius_j, int Nbonds_j);

  static std::function<double(Interaction &, MFloe &)> ruptureCriterion[4];
};
