#pragma once

/*
Même si cette application s'appelle MFloeAssembler, elle peut fonctionner
indépendament de MapleFloe.
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

#include "toofus/AABB_2D.hpp"
#include "toofus/vec2.hpp"
#include "toofus/geoPack2D.hpp"

#define MFLOE_PAC_WARN "\033[0m\033[31m\033[1m\033[4mTabaarnack !\033[24m\033[39m\033[0m: "
#define MFLOE_PAC_INFO "\033[0m\033[32m\033[1m\033[4mINFO\033[24m\033[39m\033[0m: "

struct Disk {
  vec2r pos;
  double R{0.0};
};

class MFloeAssembler {
public:
  std::vector<vec2r> zone;
  std::vector<Disk> disks;
  AABB_2D aabb_zone;
  std::string method{"none"};
  std::string output_filename{"disks.txt"};

  // radius range
  double Rmin{-1.0}, Rmax{-1.0};

  MFloeAssembler();
  void read(const char *filename);
  void run();
  
  // packing methods
  void pack_regular_triangles();
  void pack_regular_squares();
  void pack_poisson_sampling();
  
  void keepOnlyZone();
  
  int PoissonSampling_k{5000};
};