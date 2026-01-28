#include "Driving.hpp"

// Remark: if the creation fails, nullptr is returned
Driving *Driving::create(const std::string &token) {
  if (token == "imposedVelocities") {
    return new imposedVelocities();
  } /* else if (token == "____") {
    return new ____();
  }*/

  return nullptr;
}

Driving::~Driving() {}

// ========== imposedVelocities ==========

imposedVelocities::imposedVelocities() {}

void imposedVelocities::read(std::istream &is) {
  is >> vx >> vy >> vz >> vrot;
}

void imposedVelocities::write(std::ostream &os) {
  os << "imposedVelocities " << vx << ' ' << vy << ' ' << vz << ' ' << vrot << std::endl;
}

void imposedVelocities::set(double) {}
