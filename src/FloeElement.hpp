#pragma once

#include "toofus/vec2.hpp"

// A disk tile of ice
struct FloeElement {

  vec2r pos; // Position
  vec2r vel; // Velocity
  vec2r acc; // Acceleration

  // out-of-plane kinematics
  double zpos{0.0};
  double zvel{0.0};
  double zacc{0.0};

  double rot{0.0};  // Angular position
  double vrot{0.0}; // Angular velocity
  double arot{0.0}; // Angular acceleration

  double radius{0.0};
  double height{0.0};
  
  double inertia{0.0};
  double mass{0.0};

  vec2r force;        // resultant force
  double zforce{0.0}; // out-of-plane resultant force
  double moment{0.0}; // resultant moment

  FloeElement(); // Ctor
};
