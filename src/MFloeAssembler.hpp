#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

struct Disk {
  double x{0.0}, y{0.0}, R{0.0};
};

struct ZoneVertex {
  double x{0.0}, y{0.0};
};

class MFloeAssembler {
public:
  std::vector<ZoneVertex> zone;
  std::vector<Disk> disks;

  // area bounding box
  double xmin{0.0}, xmax{0.0};
  double ymin{0.0}, ymax{0.0};

  // radius range
  double Rmin{0.0}, Rmax{0.0};

  MFloeAssembler();
  void read(const char *filename);
  void run();
};