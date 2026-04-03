#pragma once

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <unordered_map>
#include <vector>

#include "MapleFloe.hpp"

#define SEE_SHOW(V)                                                                                          \
  std::cout << termcolor::color<155> << #V << "" << termcolor::reset << " = " << V << std::flush << std::endl

// toofus
#include "toofus/AABB.hpp"
#include "toofus/ColorTable.hpp"
#include "toofus/fileTool.hpp"
#include "toofus/glTools.hpp"
#include "toofus/mat4.hpp"

#include "toofus/toofus-gate/toml++/toml.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "toofus/toofus-gate/stb/stb_image_write.h"

// the conf files
MFloe Conf;
int confNum = 1;

// Avoid systematic redraws
bool needsRedraw{true};

// the zone displayed
AABB worldBox;
int fit_at_loading{1};

// display flags
int show_background{0};
int show_particles{1};
int show_forces{0};
int show_velocity_field{0};
int show_orientations{0};
int show_bonding{0};
int show_breakage_mode{0};
int show_biax_system{1};

// particle coloring
#define COLOR_NONE 0
#define COLOR_RADIUS 1
#define COLOR_VELOCITY_MAGNITUDE 2
#define COLOR_PRESSURE 3
#define COLOR_DAMAGE 4

std::vector<double> color_values;
ColorTable colorTable;
glColorBar colorBar;
int update_color_bar{1};
void setColorOption(int option);     // this will compute the particle colors
void setColor(int i, GLfloat alpha); // this will set a color depending on the selected option

ColorTable BreakageColors;
struct BreakageDataInfo {
  double t;
  size_t i;
  size_t j;
  double fnb;
  double ftb;
  double fsb;
  double coverage;
  double A;
};
std::vector<BreakageDataInfo> breakageDataInfo;

// particle stress matrices
std::vector<mat4r> Sigma;
void computeStressMatrices();

// particle initial number of bonds
std::vector<int> NbBondsRef;
void getReferenceNumberOfBonds();

// particle options
int particle_color_option{COLOR_NONE};
GLfloat alpha_particles{0.5f};
int contour_particles{1};

// arrow options
GLfloat arrow_line_width{1.0f};
double arrow_head_ratio{0.5};
double arrow_scale{1.0};

// window
int width{800};
int height{800};
float wh_ratio = (float)width / (float)height;
glTextZone textZone(3, &width, &height);
vec3i bottom_color{3, 65, 252};
vec3i top_color{136, 149, 189};

// precomputed data
double Rmin{0.0};
double Rmax{0.0};
double Rmean{0.0};
double fnMax{0.0};

// mouse
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];

void printInfo();
void preComputations();
void captureScreenshot(const char *filename);

// option configuration file
void readTomlOptions();
void saveTomlOptions();

// Drawing functions
void drawWorldBox();
void drawBiaxSystem();
void drawForces();
void drawParticles();
void drawConnectorSprings();
void drawBonding();
void drawBreakageModes();
void drawVelocityField();

// Callback functions
void keyboard(GLFWwindow *window, int key, int scancode, int action, int mods);
void mouse_button(GLFWwindow *window, int button, int action, int mods);
void cursor_pos(GLFWwindow *window, double xpos, double ypos);
void display(GLFWwindow *window);
void reshape(GLFWwindow *window, int width, int height);

// Helper functions
void printHelp();
void fit_view(GLFWwindow *window);
bool try_to_readConf(int num, MFloe &CF, int &OKNum);
void updateTextLine();