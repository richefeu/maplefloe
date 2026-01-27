#pragma once

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/freeglut.h>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>

#include "MappleFloe.hpp"

// toofus
#include "toofus/AABB.hpp"
#include "toofus/ColorTable.hpp"
#include "toofus/fileTool.hpp"
#include "toofus/glTools.hpp"

// the conf files
MFloe Conf;
int confNum = 1;

AABB worldBox;

int main_window;

// flags
int show_background     = 1; // not used
int show_particles      = 1;
int show_ghosts         = 0;
int show_period         = 1;
int show_forces         = 0;
int show_velocity_field = 0;
int show_contacts       = 0;
int showOrientations    = 0;
int showConnectors      = 1;

// particle coloring
#define COLOR_NONE 0
#define COLOR_RADIUS 1
#define COLOR_VELOCITY_MAGNITUDE 2
#define COLOR_PRESSURE 3
#define COLOR_NORMAL_STIFFNESS 4
#define COLOR_TANGENTIAL_STIFFNESS 5
#define COLOR_NORMAL_VISC_DAMPING_RATE 6
#define COLOR_FRICTION_COEFFICIENT 7
#define COLOR_ROLLING_FRICTION 8
#define COLOR_ADHESION 9
#define COLOR_GC_GLUE 10
int particle_color_option{COLOR_NONE};
std::vector<double> color_values;
ColorTable colorTable;
glColorBar colorBar;
void setColorOption(int option);     // this will compute the colors
void setColor(int i, GLfloat alpha); // this will set a color depending on the selected option
GLfloat alpha_particles = 1.0f;

// profiles
#define PROFILE_RADIUS 0
#define PROFILE_VELOCITY_X 1
int show_profile = 1;
int profile_option{PROFILE_VELOCITY_X};
std::string profile_title{"Velocity X"};
void plotProfile(std::vector<double> &prof, std::vector<double> &var);
void selectProfileToDraw();

// vectors
double arrowSize  = 0.0005;
double arrowAngle = 0.7;
double vScale     = 0.01;

double forceTubeFactor = 1.0;

// window
int width      = 800;
int height     = 800;
float wh_ratio = (float)width / (float)height;
glTextZone textZone(3, &width, &height);

// precomputed data
double Rmin{0.0};
double Rmax{0.0};
double Rmean{0.0};
double fnMax{0.0};

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];

void printInfo();
void preComputations();

// Drawing functions
void drawForces();
void drawContacts();
void drawBox();
void drawParticles();
void drawGhosts();
void drawPeriod();
void drawConnectorSprings();
void drawVelocityField();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
void printHelp();
void fit_view();
bool try_to_readConf(int num, MFloe &CF, int &OKNum);
void updateTextLine();
