#pragma once

#include <cstddef>

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

  // embedded parameters
  double meff{0.0}; // effective mass
  // double kn{0.0};   // normal stiffness
  // double kt{0.0};   // tangential stiffness
  // double mu{0.0};   // Coulomb friction coefficient
  // double muR{0.0};  // rolling friction coefficient
  // double fadh{0.0}; // adhesion force
  // double damp{0.0}; // viscuous damping rate
  double A{0.0};        // surface
  double coverage{0.0}; // Ice coverage ratio
  double Gc{0.0};       // Griffith parameter
  double dn0{0.0};      // initial distance for bonds; fnb(dn0) = 0

  Interaction();
  Interaction(size_t I, size_t J);

  void copy(Interaction &I);
};
