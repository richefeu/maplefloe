#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <random>
#include <stdexcept>
#include <vector>

struct Disk {
  double x, y, r;
};

struct Triangle {
  size_t id0, id1, id2;
  bool actif;
};

class PackingAlgorithm {
public:
  std::vector<Disk> disks;
  std::list<Triangle> triangles;
  double width, height;
  double Dmin, Dmax;
  std::mt19937 gen;

  // Paramètres pour la descente de gradient (il faudra le mettre en paramètre de la classe)
  double learning_rate{0.5};
  int max_iter{1000};
  double tol{1e-6};

  PackingAlgorithm(double w, double h, double minD, double maxD)
      : width(w), height(h), Dmin(minD), Dmax(maxD), gen(std::time(nullptr)) {}

  double genererRayon() {
    std::uniform_real_distribution<double> dist(0.5 * Dmin, 0.5 * Dmax);
    return dist(gen);
  }

  // Fonction pour calculer la distance entre deux points
  double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
  }

  // Fonction pour vérifier si trois points sont colinéaires
  bool sontColineaires(double x1, double y1, double x2, double y2, double x3, double y3, double epsilon = 1e-6) {
    double aire = 0.5 * fabs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)));
    return aire < epsilon;
  }

  // Fonction d'erreur à minimiser
  double erreur(double x, double y, double r, const Disk &d1, const Disk &d2, const Disk &d3) {
    double e1 = fabs(distance(x, y, d1.x, d1.y) - (r + d1.r));
    double e2 = fabs(distance(x, y, d2.x, d2.y) - (r + d2.r));
    double e3 = fabs(distance(x, y, d3.x, d3.y) - (r + d3.r));
    return e1 + e2 + e3;
  }

  // Fonction pour calculer l'aire signée d'un triangle formé par trois points
  double aireSignee(double x1, double y1, double x2, double y2, double x3, double y3) {
    return 0.5 * ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));
  }

  // Fonction pour vérifier si un point est à l'intérieur d'un triangle
  bool insideTriangle(const Disk &d, const Disk &d1, const Disk &d2, const Disk &d3) {
    double x = d.x, y = d.y;
    double x1 = d1.x, y1 = d1.y;
    double x2 = d2.x, y2 = d2.y;
    double x3 = d3.x, y3 = d3.y;

    // Calcul des aires signées
    double aire1 = aireSignee(x, y, x1, y1, x2, y2);
    double aire2 = aireSignee(x, y, x2, y2, x3, y3);
    double aire3 = aireSignee(x, y, x3, y3, x1, y1);

// Le point est à l'intérieur si les aires ont toutes le même signe (ou sont nulles)
    bool a1 = aire1 >= 0;
    bool a2 = aire2 >= 0;
    bool a3 = aire3 >= 0;

    bool b1 = aire1 <= 0;
    bool b2 = aire2 <= 0;
    bool b3 = aire3 <= 0;

    return (a1 && a2 && a3) || (b1 && b2 && b3);
  }

  int overlapOther(Disk &d) {
    double cx = d.x;
    double cy = d.y;

    for (size_t i = 0; i < disks.size(); ++i) {
      double dx   = disks[i].x - cx;
      double dy   = disks[i].y - cy;
      double d2   = dx * dx + dy * dy;
      double sumR = d.r + disks[i].r;
      if (d2 < sumR * sumR) { return i; }
    }
    return -1;
  }

  // Fonction pour placer un disque tangent à trois autres
  Disk place3(const Disk &d1, const Disk &d2, const Disk &d3) {
    // Vérification de la colinéarité
    if (sontColineaires(d1.x, d1.y, d2.x, d2.y, d3.x, d3.y)) {
      std::cerr << "Les trois disques sont colinéaires, impossible de placer un disque tangent.\n";
    }

    // Initialisation : centre du triangle
    double x = (d1.x + d2.x + d3.x) / 3.0;
    double y = (d1.y + d2.y + d3.y) / 3.0;
    double r = std::min({d1.r, d2.r, d3.r}) * 0.5;

    int iter = 0;
    for (; iter < max_iter; ++iter) {
      // Calcul des gradients numériques
      double dx = (erreur(x + tol, y, r, d1, d2, d3) - erreur(x - tol, y, r, d1, d2, d3)) / (2 * tol);
      double dy = (erreur(x, y + tol, r, d1, d2, d3) - erreur(x, y - tol, r, d1, d2, d3)) / (2 * tol);
      double dr = (erreur(x, y, r + tol, d1, d2, d3) - erreur(x, y, r - tol, d1, d2, d3)) / (2 * tol);

      // Mise à jour des paramètres
      x -= learning_rate * dx;
      y -= learning_rate * dy;
      r -= learning_rate * dr;

      // Vérification de la convergence
      if (erreur(x, y, r, d1, d2, d3) < tol) { break; }
    }
    std::cout << "iter = " << iter << std::endl;

    // Vérification finale de la tangence
    double e = erreur(x, y, r, d1, d2, d3);
    if (e > 0.1) {
      std::cerr << "Warning : Le disque calculé n'est pas tangent aux trois disques. Erreur : " << e << std::endl;
    }

    // Retourne le disque tangent
    return {x, y, r};
  }

  bool placerDisqueDansTriangle(std::list<Triangle>::iterator it) {
    // Triangle &t = *it;
    size_t id0 = it->id0;
    size_t id1 = it->id1;
    size_t id2 = it->id2;

    // Calculer le diamètre en fonction des 3 disques aux sommets
    const Disk &d1 = disks[id0];
    const Disk &d2 = disks[id1];
    const Disk &d3 = disks[id2];

    Disk d4 = place3(d1, d2, d3);

    double diam = 2 * d4.r;
    if (diam > Dmax) {

      double cx = 0.3333 * (d1.x + d2.x + d3.x);
      double cy = 0.3333 * (d1.y + d2.y + d3.y);

      double dst1 = distance(cx, cy, d1.x, d1.y) - d1.r;
      double dst2 = distance(cx, cy, d2.x, d2.y) - d2.r;
      double dst3 = distance(cx, cy, d3.x, d3.y) - d3.r;
      double Rmin = 0.5 * Dmin;
      if (dst1 < Rmin || dst2 < Rmin || dst3 < Rmin) {
        triangles.erase(it);
        return false;
      }

      d4.r = genererRayon();
      if (dst1 < d4.r) { d4.r = dst1; }
      if (dst2 < d4.r) { d4.r = dst2; }
      if (dst3 < d4.r) { d4.r = dst3; }

    } else if (diam < Dmin) {
      triangles.erase(it);
      return false;
    }

    int idOverlaped = overlapOther(d4);
    if (idOverlaped >= 0) {
      triangles.erase(it);
      return false;
    }

    // Ok (added)
     if (insideTriangle(d4, d1, d2, d3)){
       size_t id_disk_added = static_cast<int>(disks.size());
       subdiviserTriangle(it, id_disk_added);
       disks.push_back(d4);
       return true;   
     }
     
     return false;
    
  }

  void subdiviserTriangle(std::list<Triangle>::iterator it, size_t id_disk_added) {
    size_t id0 = it->id0;
    size_t id1 = it->id1;
    size_t id2 = it->id2;
    it         = triangles.erase(it);

    triangles.insert(it, {id0, id1, id_disk_added, true});
    triangles.insert(it, {id1, id2, id_disk_added, true});
    triangles.insert(it, {id2, id0, id_disk_added, true});
  }

  // pour devel
  void initialiserTriangulation() {
    double padding = Dmax; // Marge pour éviter que les disques ne dépassent
    disks          = {{padding, padding, genererRayon()},
                      {width - padding, padding, genererRayon()},
                      {width - padding, height - padding, genererRayon()}};

    // Créer 2 triangles pour former le rectangle
    triangles = {{0, 1, 2, true}};
  }

  // pour devel
  void genererSVG(const std::string &filename) {
    std::ofstream out(filename);
    if (!out) {
      std::cerr << "Erreur: impossible d'ouvrir " << filename << std::endl;
      return;
    }

    // Calculer la hauteur maximale
    double maxHeight = 0;
    for (const auto &d : disks) {
      if (d.y + d.r > maxHeight) { maxHeight = d.y + d.r; }
    }

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    out << "<svg width=\"" << width * 100 << "\" height=\"" << height * 100 << "\" viewBox=\"0 0 " << width << " "
        << height << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    // Dessiner le rectangle
    out << "<rect x=\"0\" y=\"0\" width=\"" << width << "\" height=\"" << height
        << "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.05\" />\n";

    // Dessiner les triangles (structure de l'algorithme) en gris transparent
    int id = 0;
    for (const auto &t : triangles) {
      const Disk &d1 = disks[t.id0];
      const Disk &d2 = disks[t.id1];
      const Disk &d3 = disks[t.id2];

      std::string fill = t.actif ? "rgba(200,200,200,0.2)" : "rgba(200,200,200,0.1)";
      out << "<polygon points=\"" << d1.x << "," << d1.y << " " << d2.x << "," << d2.y << " " << d3.x << "," << d3.y
          << "\" fill=\"" << fill << "\" stroke=\"gray\" stroke-width=\"0.01\" />\n";

      // Ajouter le numéro du triangle
      double cx = (d1.x + d2.x + d3.x) / 3;
      double cy = (d1.y + d2.y + d3.y) / 3;
      out << "<text x=\"" << cx << "\" y=\"" << cy << "\" font-size=\"0.3\" fill=\"gray\" text-anchor=\"middle\">T"
          << "T" << id++ << "</text>\n";
    }

    // Dessiner les disques (résultat final)
    id = 0;
    for (const auto &d : disks) {
      out << "<circle cx=\"" << d.x << "\" cy=\"" << d.y << "\" r=\"" << d.r
          << "\" fill=\"blue\" stroke=\"black\" stroke-width=\"0.01\" />\n";
      out << "<text x=\"" << d.x << "\" y=\"" << d.y << "\" font-size=\"" << d.r / 2
          << "\" text-anchor=\"middle\" fill=\"white\">" << id++ << "</text>\n";
    }

    // Légende
    out << "<text x=\"10\" y=\"" << height - 0.5 << "\" font-size=\"0.5\" fill=\"black\">";
    out << "Disks: " << disks.size() << " | Active triangles: ";
    int actifs = std::count_if(triangles.begin(), triangles.end(), [](const Triangle &t) { return t.actif; });
    out << actifs << " | Structure (gray triangles)</text>\n";

    out << "</svg>\n";
    out.close();
  }

  void exec() {

    // Initialisation du générateur de nombres aléatoires
    std::srand(std::time(0));

    int count = 0;
    while (!triangles.empty() && count < 1000) {
      // Tirage aléatoire d'un itérateur dans la liste
      int random_index = std::rand() % triangles.size();
      auto it          = triangles.begin();
      std::advance(it, random_index);

      // Action sur le triangle tiré au sort
      placerDisqueDansTriangle(it);

      count++;
    }
  }
};

int main() {

  double width  = 10.0;
  double height = 10.0;
  double Dmin   = 0.2;
  double Dmax   = 1;

  PackingAlgorithm pac(width, height, Dmin, Dmax);

  pac.initialiserTriangulation();
  pac.exec();

  pac.genererSVG("packing.svg");
  std::cout << "SVG image generated: packing.svg" << std::endl;

  return 0;
}
