#pragma once

#include <iostream>

class SPICE;

// Mother class for any lodaing condition
class Loading {
public:
  SPICE *box;

  static Loading *create(const std::string &token);

  virtual void read(std::istream &is) = 0;
  virtual void write(std::ostream &os) = 0;
  virtual void init() {}
  virtual void servo() {}
  virtual void velocityVerlet_halfStep1(){};
  virtual void velocityVerlet_halfStep2(){};
  virtual void forceDrivenAcceleration(){};

  //Loading() = delete; // deactivated Ctor
  virtual ~Loading(); // virtual Dtor
};

class ShearVV : public Loading {
public:
  double vx{0.0};
  double vy{0.0};

  ShearVV();

  virtual void read(std::istream &is);
  virtual void write(std::ostream &os);
  virtual void velocityVerlet_halfStep1();
};

class ShearPV : public Loading {
public:
  double pressure{0.0};
  double velocity{0.0};

  ShearPV();

  virtual void read(std::istream &is);
  virtual void write(std::ostream &os);
  virtual void init();
  virtual void velocityVerlet_halfStep1();
  virtual void velocityVerlet_halfStep2();
  virtual void forceDrivenAcceleration();
private:
  double top_mass{0.0};
  //double top_fy{0.0};
  double top_vy{0.0};
  double top_accy{0.0};
  
};
