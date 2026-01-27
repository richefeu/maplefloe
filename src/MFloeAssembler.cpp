#include "MFloeAssembler.hpp"

// Ctor
MFloeAssembler::MFloeAssembler() {}

// Read the command text-file
void MFloeAssembler::read(const char *filename) {
  // TODO
  // std::cout << "filename = " << filename << std::endl;

  std::ifstream file(filename);
  if (!file.is_open()) { std::cout /* << MFLOE_WARN*/ << "Cannot read " << filename << std::endl; }

  std::string token;
  file >> token;
  while (file.good()) {
    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(file, token); // ignore the rest of the current line
      file >> token;        // next token
      continue;
    } else if (token == "radiusRange") {
      file >> Rmin >> Rmax;
      std::cout << "Radius in range [" << Rmin << ", " << Rmax << "]" << std::endl;
    } else if (token == "Rectangle") {
      file >> xmin >> xmax >> ymin >> ymax;
      zone.clear();
      zone.push_back({xmin, ymin});
      zone.push_back({xmax, ymin});
      zone.push_back({xmax, ymax});
      zone.push_back({xmin, ymax});
      // std::cout << "Radius in range [" << Rmin << ", " << Rmax << "]" << std::endl;
    }
    // Unknown token
    else {
      std::cout << /*MFLOE_WARN <<*/ "Unknown token: " << token << std::endl;
    }

    file >> token;
  }
}

// Does the job of assembling
void MFloeAssembler::run() {

  // ...
  disks.clear();

  // Pack the disks on a rectangular zone (packing strategy will be improved after)
  double dx = 2.0 * Rmax;
  double dy = sqrt(3.0) * Rmax;
  double cx = xmin;
  double cy = ymin;

  Disk D;
  D.R     = Rmax;
  int odd = 0;
  while (cy <= ymax) {
    
    while (cx <= xmax) {
      D.x = cx;
      D.y = cy;
      disks.push_back(D);

      cx += dx;
    }
    cy += dy;
    odd = 1 - odd;
    cx  = xmin + odd * Rmax;
  }

  std::cout << "disks.size() = " << disks.size() << std::endl;
  std::ofstream txt("disks.txt");
  for (auto d : disks) {
    std::cout << d.x << " " << d.y << " " << d.R << std::endl;
    txt << d.x << " " << d.y << " " << d.R << std::endl;
  }

  // System.save(...);
}
