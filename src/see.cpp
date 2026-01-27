#include "see.hpp"

void printHelp() {
  std::cout << std::endl;
  std::cout << "+         load next configuration file" << std::endl;
  std::cout << "-         load previous configuration file" << std::endl;
  std::cout << "=         fit the view" << std::endl;
  std::cout << "a/A       particle transparency" << std::endl;
  // std::cout << "b/B       ghost particle transparency" << std::endl;
  std::cout << "c         show/hide periodic cell" << std::endl;
  std::cout << "C         show/hide contacts" << std::endl;
  std::cout << "f         show/hide normal forces" << std::endl;
  std::cout << "g         show/hide ghost particles" << std::endl;
  std::cout << "h         print this help" << std::endl;
  std::cout << "i         print information" << std::endl;
  std::cout << "m         show/hide polar plot" << std::endl;
  std::cout << "n         go to file (see terminal to enter the file number)" << std::endl;
  std::cout << "o         show/hide particle orientations" << std::endl;
  std::cout << "p         show/hide particles" << std::endl;
  std::cout << "P         switch pipe display" << std::endl;
  std::cout << "q         quit" << std::endl;
  std::cout << "s/S       tune vector sizes" << std::endl;
  std::cout << "w/W       tune displayed ghost width" << std::endl;
  std::cout << "x         show/hide surrouding material 'spider-map'" << std::endl;
  // std::cout << "x         xxxx" << std::endl;
  std::cout << std::endl;
  std::cout << "0         particles colored with light gray" << std::endl;
  std::cout << "1         particles colored with pressure" << std::endl;
  std::cout << std::endl;
}

void printInfo() {
  int V  = glutGet(GLUT_VERSION);
  int V1 = V / 10000;
  int V2 = (V - V1 * 10000) / 100;
  int V3 = V - V1 * 10000 - V2 * 100;
  std::cout << "glut version " << V1 << "." << V2 << "." << V3 << "\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case '0': {
    setColorOption(COLOR_NONE);
  } break;

  case '1': { // particle pressures
    setColorOption(COLOR_RADIUS);
  } break;

  case '2': { // particle pressures
    setColorOption(COLOR_VELOCITY_MAGNITUDE);
  } break;

  case 'a': {
    alpha_particles = std::max(0.0f, alpha_particles - 0.05f);
  } break;
  case 'A': {
    alpha_particles = std::min(1.0f, alpha_particles + 0.05f);
  } break;

  case 'b': {
    show_background = 1 - show_background;
  } break;

  case 'c': {
    show_contacts = 1 - show_contacts;
  } break;

  case 'f': {
    show_forces = 1 - show_forces;
  } break;

  case 'h': {
    printHelp();
  } break;

  case 'i': {
    printInfo();
  } break;

  case 'n': {
    std::cout << "Go to file number ";
    int conNumTry;
    std::cin >> conNumTry;
    try_to_readConf(conNumTry, Conf, confNum);
  } break;

  case 'o': {
    showOrientations = 1 - showOrientations;
  } break;

  case 'p': {
    show_particles = 1 - show_particles;
  } break;

  case 'q': {
    exit(0);
  } break;

  case 'S': {
    vScale *= 1.05;
  } break;

  case 's': {
    vScale *= 0.95;
    if (vScale < 0.0) vScale = 1.0;
  } break;

  case '-': {
    if (confNum > 0) { try_to_readConf(confNum - 1, Conf, confNum); }
    updateTextLine();
  } break;

  case '+': {
    try_to_readConf(confNum + 1, Conf, confNum);
    updateTextLine();
  } break;

  case '=': {
    fit_view();
    reshape(width, height);
  } break;
  };

  glutPostRedisplay();
}

void updateTextLine() {
  textZone.addLine("conf %d,  t %0.4g s", confNum, Conf.t);
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    // display();
    glutPostRedisplay();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON: {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
        mouse_mode = PAN;
      } else {
        mouse_mode = ROTATION;
      }
    } break;
    case GLUT_MIDDLE_BUTTON: {
      mouse_mode = ZOOM;
    } break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING) { return; }

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case ZOOM: {
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    double ddx = (worldBox.max.x - worldBox.min.x) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x += ddx;
    worldBox.min.y -= ddy;
    worldBox.max.y += ddy;
  } break;

  case PAN: {
    double ddx = (worldBox.max.x - worldBox.min.x) * dx;
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x -= ddx;
    worldBox.min.y += ddy;
    worldBox.max.y += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  glutPostRedisplay();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  //glTools::clearBackground((bool)show_background);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (show_particles == 1) { drawParticles(); }
  //if (show_ghosts == 1) { drawGhosts(); }
  //if (showConnectors == 1) { drawConnectorSprings(); }
  if (show_contacts == 1) { drawContacts(); }
  if (show_forces == 1) { drawForces(); }
  if (show_velocity_field == 1) { drawVelocityField(); }
  //if (show_period == 1) { drawPeriod(); }

  if (particle_color_option > COLOR_NONE) { colorBar.show(width, height, colorTable); }

  // profile devel
  //if (show_profile == 1) { selectProfileToDraw(); }

  textZone.draw();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  // FIXME: woldBox fit (use the Boundaries also ?)
  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {

    vec2r pos = Conf.FloeElements[i].pos;
    double R  = Conf.FloeElements[i].radius;

    double xmin = pos.x - R;
    double xmax = pos.x + R;
    double ymin = pos.y - R;
    double ymax = pos.y + R;

    worldBox.min.x = std::min(worldBox.min.x, xmin);
    worldBox.min.y = std::min(worldBox.min.y, ymin);
    worldBox.max.x = std::max(worldBox.max.x, xmax);
    worldBox.max.y = std::max(worldBox.max.y, ymax);
  }

  reshape(width, height);
}

void reshape(int w, int h) {
  width  = w;
  height = h;

  double left   = worldBox.min.x;
  double right  = worldBox.max.x;
  double bottom = worldBox.min.y;
  double top    = worldBox.max.y;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW     = 0.1 * worldW;
  double dH     = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;

  if (worldW > worldH) {
    worldH = worldW * ((GLfloat)height / (GLfloat)width);
    top    = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * ((GLfloat)width / (GLfloat)height);
    right  = 0.5 * (left + right + worldW);
    left   = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
}


// set and compute the colors of particles
void setColorOption(int option) {
  particle_color_option = option;
  color_values.resize(Conf.FloeElements.size(), 0.0);

  switch (particle_color_option) {
  case COLOR_NONE: {
    // nothing to do
  } break;

  case COLOR_RADIUS: {
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) { color_values[i] = Conf.FloeElements[i].radius; }
    colorTable.setMinMax(Rmin, Rmax);
    colorBar.setTitle("Radius");
  } break;

  case COLOR_VELOCITY_MAGNITUDE: {
    double Vmin = 0.0;
    double Vmax = 0.0;
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
      double V        = sqrt(Conf.FloeElements[i].vel * Conf.FloeElements[i].vel);
      color_values[i] = V;
      if (V < Vmin) { Vmin = V; }
      if (V > Vmax) { Vmax = V; }
    }
    colorTable.setMinMax(Vmin, Vmax);
    colorBar.setTitle("Vel. Mag.");
  } break;

  /*
  case COLOR_NORMAL_STIFFNESS: {
    double normalStiffnessMin = Conf.Particles[0].normalStiffness;
    double normalStiffnessMax = normalStiffnessMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double normalStiffness = Conf.Particles[i].normalStiffness;
      color_values[i]        = normalStiffness;
      if (normalStiffness < normalStiffnessMin) { normalStiffnessMin = normalStiffness; }
      if (normalStiffness > normalStiffnessMax) { normalStiffnessMax = normalStiffness; }
    }
    colorTable.setMinMax(normalStiffnessMin, normalStiffnessMax);
    colorBar.setTitle("Normal stiffness");
  } break;

  case COLOR_TANGENTIAL_STIFFNESS: {
    double tangentialStiffnessMin = Conf.Particles[0].tangentialStiffness;
    double tangentialStiffnessMax = tangentialStiffnessMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double tangentialStiffness = Conf.Particles[i].tangentialStiffness;
      color_values[i]            = tangentialStiffness;
      if (tangentialStiffness < tangentialStiffnessMin) { tangentialStiffnessMin = tangentialStiffness; }
      if (tangentialStiffness > tangentialStiffnessMax) { tangentialStiffnessMax = tangentialStiffness; }
    }
    colorTable.setMinMax(tangentialStiffnessMin, tangentialStiffnessMax);
    colorBar.setTitle("Tangential stiffness");
  } break;

  case COLOR_NORMAL_VISC_DAMPING_RATE: {
    double normalViscDampingRateMin = Conf.Particles[0].normalViscDampingRate;
    double normalViscDampingRateMax = normalViscDampingRateMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double normalViscDampingRate = Conf.Particles[i].normalViscDampingRate;
      color_values[i]              = normalViscDampingRate;
      if (normalViscDampingRate < normalViscDampingRateMin) { normalViscDampingRateMin = normalViscDampingRate; }
      if (normalViscDampingRate > normalViscDampingRateMax) { normalViscDampingRateMax = normalViscDampingRate; }
    }
    colorTable.setMinMax(normalViscDampingRateMin, normalViscDampingRateMax);
    colorBar.setTitle("Normal viscous damping rate");
  } break;

  case COLOR_FRICTION_COEFFICIENT: {
    double frictionMin = Conf.Particles[0].friction;
    double frictionMax = frictionMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double friction = Conf.Particles[i].friction;
      color_values[i] = friction;
      if (friction < frictionMin) { frictionMin = friction; }
      if (friction > frictionMax) { frictionMax = friction; }
    }
    colorTable.setMinMax(frictionMin, frictionMax);
    colorBar.setTitle("Friction");
  } break;

  case COLOR_ROLLING_FRICTION: {
    double rollingFrictionMin = Conf.Particles[0].rollingFriction;
    double rollingFrictionMax = rollingFrictionMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double rollingFriction = Conf.Particles[i].rollingFriction;
      color_values[i]        = rollingFriction;
      if (rollingFriction < rollingFrictionMin) { rollingFrictionMin = rollingFriction; }
      if (rollingFriction > rollingFrictionMax) { rollingFrictionMax = rollingFriction; }
    }
    colorTable.setMinMax(rollingFrictionMin, rollingFrictionMax);
    colorBar.setTitle("Rolling friction");
  } break;

  case COLOR_ADHESION: {
    double adhesionMin = Conf.Particles[0].adhesion;
    double adhesionMax = adhesionMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double adhesion = Conf.Particles[i].adhesion;
      color_values[i] = adhesion;
      if (adhesion < adhesionMin) { adhesionMin = adhesion; }
      if (adhesion > adhesionMax) { adhesionMax = adhesion; }
    }
    colorTable.setMinMax(adhesionMin, adhesionMax);
    colorBar.setTitle("Adhesion");
  } break;

  case COLOR_GC_GLUE: {
    double GcGlueMin = Conf.Particles[0].GcGlue;
    double GcGlueMax = GcGlueMin;
    for (size_t i = 0; i < Conf.Particles.size(); ++i) {
      double GcGlue   = Conf.Particles[i].GcGlue;
      color_values[i] = GcGlue;
      if (GcGlue < GcGlueMin) { GcGlueMin = GcGlue; }
      if (GcGlue > GcGlueMax) { GcGlueMax = GcGlue; }
    }
    colorTable.setMinMax(GcGlueMin, GcGlueMax);
    colorBar.setTitle("Gc glue");
  } break;
*/
  
  default:
    break;
  }
}

void setColor(int i, GLfloat alpha) {
  colorRGBA col;
  if (particle_color_option == COLOR_NONE) {
    col.r = col.g = col.b = 200;
  } else {
    colorTable.getRGB(color_values[i], &col);
  }
  glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, alpha);
}



void drawParticles() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
    vec2r pos = Conf.FloeElements[i].pos;
    double R  = Conf.FloeElements[i].radius;

    setColor(i, alpha_particles);
    glBegin(GL_POLYGON);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
    glBegin(GL_LINE_LOOP);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    if (showOrientations) {
      double rot = Conf.FloeElements[i].rot;
      glBegin(GL_LINES);
      glVertex2f(pos.x, pos.y);
      glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
      glEnd();
    }
  }
}

void drawVelocityField() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(1.0f);
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

  double v2max = Conf.FloeElements[0].vel * Conf.FloeElements[0].vel;
  for (size_t i = 1; i < Conf.FloeElements.size(); ++i) {
    double v2 = Conf.FloeElements[i].vel * Conf.FloeElements[i].vel;
    if (v2 > v2max) { v2max = v2; }
  }
  double vmax     = sqrt(v2max);
  double factor   = 2.0 * Rmean / vmax;
  double headSize = 0.5 * Rmean;

  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
    vec2r pos = Conf.FloeElements[i].pos;
    vec2r vel = Conf.FloeElements[i].vel;
    vec2r posHead(pos.x + factor * vel.x, pos.y + factor * vel.y);
    vec2r n = posHead - pos;
    n.normalize();
    vec2r t(-n.y * 0.6, n.x * 0.6);

    glBegin(GL_LINES);
    glVertex2f(pos.x, pos.y);
    glVertex2f(posHead.x, posHead.y);

    glVertex2f(posHead.x, posHead.y);
    glVertex2f(posHead.x + headSize * (-n.x + t.x), posHead.y + headSize * (-n.y + t.y));

    glVertex2f(posHead.x, posHead.y);
    glVertex2f(posHead.x + headSize * (-n.x - t.x), posHead.y + headSize * (-n.y - t.y));
    glEnd();
  }
}

void drawContacts() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(1.5f);
  // double Xmid = 0.5 * (Conf.xmin + Conf.xmax);
  double L     = Conf.xmax - Conf.xmin;
  double halfL = 0.5 * L;

  // grain-grain
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i   = Conf.Interactions[k].i;
    size_t j   = Conf.Interactions[k].j;
    vec2r posi = Conf.FloeElements[i].pos;
    vec2r sij  = Conf.FloeElements[j].pos - Conf.FloeElements[i].pos;
    if (sij.x > halfL) {
      sij.x -= L;
    } else if (sij.x < -halfL) {
      sij.x += L;
    }
    vec2r posj = posi + sij;
    glVertex2f(posi.x, posi.y);
    glVertex2f(posj.x, posj.y);
  }
  glEnd();
}

// remove periodic things
void drawForces() {
  if (mouse_mode != NOTHING) { return; }

  double L     = Conf.xmax - Conf.xmin;
  double halfL = 0.5 * L;

  // grain-grain
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i   = Conf.Interactions[k].i;
    size_t j   = Conf.Interactions[k].j;
    vec2r posi = Conf.FloeElements[i].pos;
    vec2r sij  = Conf.FloeElements[j].pos - Conf.FloeElements[i].pos;
    if (sij.x > halfL) {
      sij.x -= L;
    } else if (sij.x < -halfL) {
      sij.x += L;
    }
    vec2r posj = posi + sij;

    // Calculate the width of the rectangle
    GLfloat width = Rmean * (Conf.Interactions[k].fn / fnMax);
    // GLfloat width = 0.0001 * Conf.Interactions[k].fn;

    // Calculate the direction vector and the perpendicular vector
    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    // Calculate the four corners of the rectangle
    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    // Draw the filled rectangle
    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }
}

bool try_to_readConf(int num, MFloe &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << file_name << std::endl;
    OKNum = num;
    CF.loadConf(file_name);
    CF.accelerations();
    preComputations();
    setColorOption(particle_color_option);
    return true;
  } else {
    std::cout << file_name << " does not exist" << std::endl;
  }
  return false;
}

void preComputations() {
  if (Conf.FloeElements.empty()) { return; }

  Rmin  = Conf.FloeElements[0].radius;
  Rmax  = Rmin;
  Rmean = 0.0;
  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
    double R = Conf.FloeElements[i].radius;
    Rmean += R;
    if (R < Rmin) { Rmin = R; }
    if (R > Rmax) { Rmax = R; }
  }
  Rmean /= static_cast<double>(Conf.FloeElements.size());

  fnMax = 0.0;
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    if (Conf.Interactions[k].fn > fnMax) { fnMax = Conf.Interactions[k].fn; }
  }
}

void menu(int num) {
  switch (num) {

  case 0: {
    exit(0);
  } break;
  case 1: {
    Conf.updateNeighbors(Conf.dVerlet);
  } break;
  case 2: {
    std::vector<Interaction> storedInteractions(Conf.Interactions.size());
    std::copy(Conf.Interactions.begin(), Conf.Interactions.end(), storedInteractions.begin());
    Conf.Interactions.clear();
    for (size_t k = 0; k < storedInteractions.size(); k++) {
      if (fabs(storedInteractions[k].fn) < 1e-20 && storedInteractions[k].isBonded == false) { continue; }
      Conf.Interactions.push_back(storedInteractions[k]);
    }
  } break;

  // submenu10
  case 10: {
    show_particles = 1 - show_particles;
  } break;
  
  case 13: {
    show_forces = 1 - show_forces;
  } break;
  case 14: {
    show_contacts = 1 - show_contacts;
  } break;
  case 15: {
    showOrientations = 1 - showOrientations;
  } break;
  
  case 17: {
    show_velocity_field = 1 - show_velocity_field;
  } break;

  // submenu40
  case 40: {
    setColorOption(COLOR_NONE);
  } break;
  case 41: {
    setColorOption(COLOR_RADIUS);
  } break;
  case 42: {
    setColorOption(COLOR_VELOCITY_MAGNITUDE);
  } break;
  
  // submenu70
  case 70: {
    double Z     = 0.0;
    double Zprox = 0.0;
    for (size_t k = 0; k < Conf.Interactions.size(); k++) {
      if (Conf.Interactions[k].fn > 0.0) { Z += 2.0; } // FIXME: not ok for glued interactions
      Zprox += 2.0;
    }
    double N = static_cast<double>(Conf.FloeElements.size());
    if (N > 1.0) {
      Z /= N;
      Zprox /= N;
    }
    textZone.addLine("Z = %.4f, Zprox = %.4f", Z, Zprox);
  } break;

  // submenu100
  // ...

  }; // end switch

  glutPostRedisplay();
}

void buildMenu() {
  int submenu10 = glutCreateMenu(menu);
  glutAddMenuEntry("Show particles", 10);
  //glutAddMenuEntry("Show ghost-particles", 11);
  //glutAddMenuEntry("Show period", 12);
  glutAddMenuEntry("Show forces", 13);
  glutAddMenuEntry("Show contacts", 14);
  glutAddMenuEntry("Show orientations", 15);
  //glutAddMenuEntry("Show connectors", 16);
  glutAddMenuEntry("Show velocity field", 17);

  int submenu40 = glutCreateMenu(menu);
  glutAddMenuEntry("None", 40);
  glutAddMenuEntry("Radius", 41);
  glutAddMenuEntry("Velocity Magnitude", 42);
  //glutAddMenuEntry("Normal stiffness", 43);
  //glutAddMenuEntry("Tangential stiffness", 44);
  //glutAddMenuEntry("Normal viscous damping rate", 45);
  //glutAddMenuEntry("Friction coefficient", 46);
  //glutAddMenuEntry("Rolling friction", 47);
  //glutAddMenuEntry("Adhesion", 48);
  //glutAddMenuEntry("Gc glue", 49);

  int submenu70 = glutCreateMenu(menu);
  glutAddMenuEntry("Connectivity", 70);

  //int submenu100 = glutCreateMenu(menu);
  //glutAddMenuEntry("None", 100);
  //glutAddMenuEntry("Radius", 101);
  //glutAddMenuEntry("Velocity X", 102);

  glutCreateMenu(menu); // Main menu
  glutAddSubMenu("Display options", submenu10);
  glutAddSubMenu("Color particles", submenu40);
  glutAddSubMenu("Extract data", submenu70);
  //glutAddSubMenu("Profiles", submenu100);
  glutAddMenuEntry("Quit", 0);
  glutAddMenuEntry("Update the list of neighbors", 1);
  glutAddMenuEntry("Clean the list of neighbors", 2);
}

// =====================================================================
// Main function
// =====================================================================
int main(int argc, char *argv[]) {

  if (argc == 1) {
    confNum = 0;
    std::cout << "Current Configuration: ";
    try_to_readConf(confNum, Conf, confNum);
  } else if (argc == 2) {
    if (fileTool::containsOnlyDigits(argv[1])) {
      confNum = std::atoi(argv[1]);
      std::cout << "Current Configuration: ";
      try_to_readConf(confNum, Conf, confNum);
    } else {
      Conf.loadConf(argv[1]);
      setColorOption(particle_color_option);
      std::cout << "Current Configuration: " << argv[1] << std::endl;
      confNum = Conf.iconf;
    }
  }

  mouse_mode = NOTHING;

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  int X0 = (glutGet(GLUT_SCREEN_WIDTH) - width) / 2;
  int Y0 = (glutGet(GLUT_SCREEN_HEIGHT) - height) / 2;
  glutInitWindowPosition(X0, Y0);
  glutInitWindowSize(width, height);

  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // ==== Other initialisations
  glText::init();
  updateTextLine();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
