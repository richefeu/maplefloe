#pragma once

#include <cstddef>
//#include "toofus/vec2.hpp"

struct Interaction {
  size_t i{0};
  size_t j{0};

  // flags
  bool isBonded{false};           // Indicates if the particle has any bonds or connections to other particles
  //bool isSameMaterialBond{false}; // true means bonded by the same material, otherwise false

  // forces
  double fn{0.0};  // contact normal force
  double fnb{0.0}; // bond normal force
  double ft{0.0};  // contact tangential force
  double ftb{0.0}; // bond tangential force

  // embedded parameters
  double meff{0.0}; // effective mass
  double kn{0.0};   // normal stiffness
  double kt{0.0};   // tangential stiffness
  double mu{0.0};   // Coulomb friction coefficient
  double muR{0.0};  // rolling friction coefficient
  double fadh{0.0}; // adhesion force
  double damp{0.0}; // viscuous damping rate
  double Gs{0.0};   //
  double dn0{0.0};  // initial distance for bonds; fnb(dn0) = 0

  Interaction();
  Interaction(size_t I, size_t J);

  void copy(Interaction &I);
};
