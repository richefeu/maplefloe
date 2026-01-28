#pragma once

#include <iostream>

// Mother class for any lodaing condition
class Driving {
public:
  double vx{0.0};
  double vy{0.0};
  double vz{0.0};
  double vrot{0.0};
  
  static Driving *create(const std::string &token);

  virtual void read(std::istream &is) = 0;
  virtual void write(std::ostream &os) = 0;
  virtual void set(double) {}

  //Drive() = delete; // deactivated Ctor
  virtual ~Driving(); // virtual Dtor
};

class imposedVelocities : public Driving {
public:
  
  imposedVelocities();

  virtual void read(std::istream &is);
  virtual void write(std::ostream &os);
  virtual void set(double);
};


