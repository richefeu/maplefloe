#include "MFloeAssembler.hpp"

// Ctor
MFloeAssembler::MFloeAssembler() {
  srand(time(NULL));
}

// Read the command text-file
void MFloeAssembler::read(const char *filename) {
  std::ifstream file(filename);
  if (!file.is_open()) { std::cout << MFLOE_PAC_WARN << "Cannot read " << filename << std::endl; }
  read(file);
}

// Read the command text-file
void MFloeAssembler::read(std::istream &is) {

  std::string token;
  is >> token;
  while (is.good()) {
    if (token[0] == '/' || token[0] == '#' || token[0] == '!') {
      getline(is, token); // ignore the rest of the current line
      is >> token;        // next token
      continue;
    } else if (token == "radiusRange") {
      is >> Rmin >> Rmax;
      std::cout << MFLOE_PAC_INFO << "Radius in range [" << Rmin << ", " << Rmax << "]" << std::endl;
    } else if (token == "radius") {
      is >> Rmin;
      Rmax = Rmin;
      std::cout << MFLOE_PAC_INFO << "Single radius is " << Rmin << std::endl;
    } /*else if (token == "height") {
      is >> height;
      std::cout << MFLOE_PAC_INFO << "common height is " << height << std::endl;
    } */ else if (token == "rectangleZone") {
      is >> aabb_zone.min.x >> aabb_zone.min.y >> aabb_zone.max.x >> aabb_zone.max.y;
      zone.clear();
      zone.push_back({aabb_zone.min.x, aabb_zone.min.y});
      zone.push_back({aabb_zone.max.x, aabb_zone.min.y});
      zone.push_back({aabb_zone.max.x, aabb_zone.max.y});
      zone.push_back({aabb_zone.min.x, aabb_zone.max.y});
      std::cout << MFLOE_PAC_INFO << "Rectangle zone is from (" << aabb_zone.min.x << ", " << aabb_zone.min.y
                << ") to (" << aabb_zone.max.x << ", " << aabb_zone.max.y << ")" << std::endl;
    } else if (token == "anyZone") {
      size_t nbVertices = 0;
      is >> nbVertices;
      zone.clear();
      vec2r vertex;
      for (size_t v = 0; v < nbVertices; v++) {
        is >> vertex;
        zone.push_back(vertex);
      }
      aabb_zone.set_single(zone[0]);
      for (size_t v = 1; v < zone.size(); v++) { aabb_zone.add(zone[v]); }
      std::cout << MFLOE_PAC_INFO << "zone shape is from (" << aabb_zone.min.x << ", " << aabb_zone.min.y << ") to ("
                << aabb_zone.max.x << ", " << aabb_zone.max.y << ")" << std::endl;
    } else if (token == "method") {
      is >> method;
      std::cout << MFLOE_PAC_INFO << "Packing method is " << method << std::endl;
    } else if (token == "PoissonSampling_k") {
      is >> PoissonSampling_k;
      std::cout << MFLOE_PAC_INFO << "PoissonSampling_k is " << PoissonSampling_k << std::endl;
    } else if (token == "output_filename") {
      is >> output_filename;
      std::cout << MFLOE_PAC_INFO << "Output filename is " << output_filename << std::endl;
    }
    // Unknown token
    else {
      std::cout << MFLOE_PAC_WARN << "Unknown token: " << token << std::endl;
    }

    is >> token;
  }
}

void MFloeAssembler::pack_regular_triangles() {
  disks.clear();
  Disk D;

  double l = 0.0;
  for (double y = aabb_zone.min.y + Rmax; y <= aabb_zone.max.y - Rmax; y += 1.73205081 * Rmax) {
    for (double x = aabb_zone.min.x + Rmax + l * Rmax; x <= aabb_zone.max.x - Rmax; x += 2 * Rmax) {
      D.R     = Rmin + (Rmax - Rmin) * rand() / (double)RAND_MAX;
      D.pos.x = x;
      D.pos.y = y;
      disks.push_back(D);
    }
    l = 1.0 - l;
  }

  std::cout << "Number of disks after Poisson sampling: " << disks.size() << std::endl;
  keepOnlyZone();
  std::cout << "Number of disks after keeping only zone: " << disks.size() << std::endl;
}

void MFloeAssembler::pack_regular_squares() {
  disks.clear();
  Disk D;

  for (double y = aabb_zone.min.y + Rmax; y <= aabb_zone.max.y - Rmax; y += 2 * Rmax) {
    for (double x = aabb_zone.min.x + Rmax; x <= aabb_zone.max.x - Rmax; x += 2 * Rmax) {
      D.R     = Rmin + (Rmax - Rmin) * rand() / (double)RAND_MAX;
      D.pos.x = x;
      D.pos.y = y;
      disks.push_back(D);
    }
  }

  std::cout << "Number of disks after Poisson sampling: " << disks.size() << std::endl;
  keepOnlyZone();
  std::cout << "Number of disks after keeping only zone: " << disks.size() << std::endl;
}

void MFloeAssembler::pack_poisson_sampling() {

  GeoPack2D GP(Rmin, Rmax, PoissonSampling_k, aabb_zone.min.x, aabb_zone.max.x, aabb_zone.min.y, aabb_zone.max.y);
  GP.seedTime();
  GP.exec();

  Disk D;
  for (size_t i = 0; i < GP.sample.size(); i++) {
    D.pos.x = GP.sample[i].x;
    D.pos.y = GP.sample[i].y;
    D.R     = GP.sample[i].r;
    disks.push_back(D);
  }

  std::cout << "Number of disks after Poisson sampling: " << disks.size() << std::endl;
  keepOnlyZone();
  std::cout << "Number of disks after keeping only zone: " << disks.size() << std::endl;
}

void MFloeAssembler::keepOnlyZone() {
  std::vector<Disk> okDisks;

  vec2r ray_dir(0.937, 0.349);

  for (size_t d = 0; d < disks.size(); d++) {
    bool inside = false;
    for (size_t i = 0, j = zone.size() - 1; i < zone.size(); j = i++) {
      const vec2r &A = zone[i];
      const vec2r &B = zone[j];
      const vec2r &P = disks[d].pos;
      if (((A.y > P.y) != (B.y > P.y)) && (P.x < (B.x - A.x) * (P.y - A.y) / (B.y - A.y) + A.x)) { inside = !inside; }
    }
    if (inside) { okDisks.push_back(disks[d]); }
  }

  disks.clear();
  for (size_t i = 0; i < okDisks.size(); i++) { disks.push_back(okDisks[i]); }
}

// Does the job of assembling
void MFloeAssembler::run() {

  if (zone.empty()) {
    std::cout << MFLOE_PAC_WARN << "You need to define a zone for packing" << std::endl;
    return;
  }

  if (Rmin < 0.0 || Rmax < 0.0) {
    std::cout << MFLOE_PAC_WARN << "You need to define radius range" << std::endl;
    return;
  }

  if (method == "none") {
    std::cout << MFLOE_PAC_WARN << "You need to define a method for packing" << std::endl;
    return;
  } else if (method == "regular-triangles") {
    pack_regular_triangles();
  } else if (method == "regular-squares") {
    pack_regular_squares();
  } else if (method == "Poisson-sampling") {
    pack_poisson_sampling();
  } else {
    std::cout << MFLOE_PAC_WARN << "This packing method '" << method << "' is unknown" << std::endl;
    return;
  }

  std::cout << MFLOE_PAC_INFO << "Number of packed disks = " << disks.size() << std::endl;
  std::ofstream txt(output_filename);
  for (auto d : disks) { txt << d.pos.x << " " << d.pos.y << " " << d.R << std::endl; }
}
