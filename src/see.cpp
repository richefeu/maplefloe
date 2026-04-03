#include "see.hpp"

void readTomlOptions() {
  if (!fileTool::fileExists("see-options.txt")) {
    saveTomlOptions();
    return;
  }

  toml::table tbl = toml::parse_file("see-options.txt");

  if (tbl.contains("display")) {
    show_particles      = tbl["display"]["show_particles"].value_or(show_particles);
    show_velocity_field = tbl["display"]["show_velocity_field"].value_or(show_velocity_field);
    show_orientations   = tbl["display"]["show_orientations"].value_or(show_orientations);
    show_bonding        = tbl["display"]["show_bonding"].value_or(show_bonding);
    show_background     = tbl["display"]["show_background"].value_or(show_background);
    show_forces         = tbl["display"]["show_forces"].value_or(show_forces);

    if (tbl["display"].as_table()->contains("bottom_color")) {
      auto color_array = tbl["display"]["bottom_color"].as_array();
      if (color_array && color_array->size() >= 3) {
        bottom_color.x = (*color_array)[0].value_or(bottom_color.x);
        bottom_color.y = (*color_array)[1].value_or(bottom_color.y);
        bottom_color.z = (*color_array)[2].value_or(bottom_color.z);
      }
    }

    if (tbl["display"].as_table()->contains("top_color")) {
      auto color_array = tbl["display"]["top_color"].as_array();
      if (color_array && color_array->size() >= 3) {
        top_color.x = (*color_array)[0].value_or(top_color.x);
        top_color.y = (*color_array)[1].value_or(top_color.y);
        top_color.z = (*color_array)[2].value_or(top_color.z);
      }
    }
  }

  if (tbl.contains("particles")) {
    alpha_particles       = tbl["particles"]["alpha_particles"].value_or(alpha_particles);
    contour_particles     = tbl["particles"]["contour_particles"].value_or(contour_particles);
    particle_color_option = tbl["particles"]["particle_color_option"].value_or(particle_color_option);
  }

  if (tbl.contains("arrows")) {
    arrow_line_width = tbl["arrows"]["arrow_line_width"].value_or(arrow_line_width);
    arrow_head_ratio = tbl["arrows"]["arrow_head_ratio"].value_or(arrow_head_ratio);
    arrow_scale      = tbl["arrows"]["arrow_scale"].value_or(arrow_scale);
  }

  if (tbl.contains("window")) {
    width  = tbl["window"]["width"].value_or(width);
    height = tbl["window"]["height"].value_or(height);
  }

  if (tbl.contains("view")) {
    fit_at_loading = tbl["view"]["fit_at_loading"].value_or(1);
    worldBox.min.x = tbl["view"]["xmin"].value_or(worldBox.min.x);
    worldBox.max.x = tbl["view"]["xmax"].value_or(worldBox.max.x);
    worldBox.min.y = tbl["view"]["ymin"].value_or(worldBox.min.y);
    worldBox.max.y = tbl["view"]["ymax"].value_or(worldBox.max.y);
  }
}

void saveTomlOptions() {
  // it overwrites the file is it exists;

  fit_at_loading = 0; // so that the next time the view is saved

  // clang-format off
  auto tbl = toml::table{
      {"display",
       toml::table{
         {"show_particles", show_particles},
         {"show_velocity_field", show_velocity_field},
         {"show_orientations", show_orientations},
         {"show_bonding", show_bonding},
         {"show_background", show_background},
         {"show_forces", show_forces},
         {"bottom_color", toml::array{bottom_color.x, bottom_color.y, bottom_color.z}},
         {"top_color", toml::array{top_color.x, top_color.y, top_color.z}},
       }
      },
      {"window",
       toml::table{
         {"width", width},
         {"height", height},
       }
      },
      {"view",
       toml::table{
         {"fit_at_loading", fit_at_loading},
         {"xmin", worldBox.min.x},
         {"xmax", worldBox.max.x},
         {"ymin", worldBox.min.y},
         {"ymax", worldBox.max.y},
       }
      },
      {"particles",
       toml::table{
         {"alpha_particles", alpha_particles},
         {"contour_particles", contour_particles},
         {"particle_color_option", particle_color_option},  
       }
      },
      {"arrows",
       toml::table{
         {"arrow_line_width", arrow_line_width},
         {"arrow_head_ratio", arrow_head_ratio},
         {"arrow_scale", arrow_scale},
       }
      },
  };
  // clang-format on

  std::ofstream file("see-options.txt");
  file << tbl << std::endl;
}

void printHelp() {
  std::cout << std::endl;
  std::cout << "->        load next configuration file" << std::endl;
  std::cout << "<-        load previous configuration file" << std::endl;
  std::cout << "=         fit the view" << std::endl;
  std::cout << "a/A       particle transparency" << std::endl;
  std::cout << "h         print this help" << std::endl;
  std::cout << "q         quit" << std::endl;
}

void printInfo() {
  const char *version = glfwGetVersionString();
  std::cout << "GLFW version " << version << "\n";
}

void captureScreenshot(const char *filename) {
  // Allouer un buffer pour stocker les pixels
  unsigned char *pixels = new unsigned char[width * height * 3];

  // Lire les pixels depuis le framebuffer
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

  // Inverser les lignes (OpenGL stocke les pixels de bas en haut)
  unsigned char *flipped_pixels = new unsigned char[width * height * 3];
  for (int y = 0; y < height; y++) {
    memcpy(flipped_pixels + (height - 1 - y) * width * 3, pixels + y * width * 3, width * 3);
  }

  // Sauvegarder l'image
  stbi_write_png(filename, width, height, 3, flipped_pixels, width * 3);

  // Libérer la mémoire
  delete[] pixels;
  delete[] flipped_pixels;
}

void keyboard(GLFWwindow *window, int key, int /*scancode*/, int action, int mods) {
  if (action != GLFW_PRESS) return;

  // SEE_SHOW(key);
  // SEE_SHOW(scancode);
  // SEE_SHOW(action);
  // SEE_SHOW(mods);

  switch (key) {

  case GLFW_KEY_SPACE: {
    if (mods == GLFW_MOD_SHIFT) {
      readTomlOptions();
      glfwSetWindowSize(window, width, height);
      std::cout << MFLOE_INFO << "Loaded 'see-options.txt'" << std::endl;
      textZone.addLine("Loaded 'see-options.txt'");
    } else {
      saveTomlOptions();
      std::cout << MFLOE_INFO << "Saved 'see-options.txt'" << std::endl;
      textZone.addLine("Saved 'see-options.txt'");
    }
  } break;

  case GLFW_KEY_0: {
    setColorOption(COLOR_NONE);
    textZone.addLine("particles not colored");
  } break;

  case GLFW_KEY_1: {
    setColorOption(COLOR_RADIUS);
    textZone.addLine("color particles by radius");
  } break;

  case GLFW_KEY_2: {
    setColorOption(COLOR_VELOCITY_MAGNITUDE);
    textZone.addLine("color particles by velocity magnitude");
  } break;

  case GLFW_KEY_3: {
    setColorOption(COLOR_PRESSURE);
    textZone.addLine("color particles by pressure");
  } break;

  case GLFW_KEY_4: {
    setColorOption(COLOR_DAMAGE);
    textZone.addLine("color particles by damage (1 - Nb/Nb0)");
  } break;

  case GLFW_KEY_Q: { // 'a/A'
    if (mods == GLFW_MOD_SHIFT) {
      alpha_particles = std::min(1.0f, alpha_particles + 0.05f);
    } else {
      alpha_particles = std::max(0.0f, alpha_particles - 0.05f);
    }
    SEE_SHOW(alpha_particles);
    textZone.addLine("alpha_particles = %g", alpha_particles);
  } break;

  case GLFW_KEY_B: {
    if (mods == GLFW_MOD_SHIFT) {
      show_background = 1 - show_background;
      textZone.addLine("show_background = %d", show_background);
    } else {
      show_bonding = 1 - show_bonding;
      SEE_SHOW(show_bonding);
      textZone.addLine("show_bonding = %d", show_bonding);

      BreakageColors.setMinMax(-1.0f, 1.0f);
    }
  } break;

  case GLFW_KEY_C: {
    contour_particles = 1 - contour_particles;
  } break;

  case GLFW_KEY_F: {
    show_forces = 1 - show_forces;
    SEE_SHOW(show_forces);
    textZone.addLine("show_forces = %d", show_forces);
  } break;

  case GLFW_KEY_H: {
    printHelp();
  } break;

  case GLFW_KEY_I: {
    printInfo();
  } break;

  case GLFW_KEY_SEMICOLON: { // 'm'
    show_breakage_mode = 1 - show_breakage_mode;
    SEE_SHOW(show_breakage_mode);
    textZone.addLine("show_breakage_mode = %d", show_breakage_mode);
  } break;

  case GLFW_KEY_O: {
    show_orientations = 1 - show_orientations;
    SEE_SHOW(show_orientations);
    textZone.addLine("show_orientations = %d", show_orientations);
  } break;

  case GLFW_KEY_P: {
    show_particles = 1 - show_particles;
    SEE_SHOW(show_particles);
    textZone.addLine("show_particles = %d", show_particles);
  } break;

  case GLFW_KEY_A: { // 'q/Q'
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  } break;

  case GLFW_KEY_S: {
    if (mods == GLFW_MOD_SHIFT) {
      arrow_scale *= 1.05;
    } else {
      arrow_scale *= 0.95;
      if (arrow_scale < 0.0) arrow_scale = 1.0;
    }
    SEE_SHOW(arrow_scale);
    textZone.addLine("arrow_scale = %g", arrow_scale);
  } break;

  case GLFW_KEY_V: {
    show_velocity_field = 1 - show_velocity_field;
    SEE_SHOW(show_velocity_field);
    textZone.addLine("show_velocity_field = %d", show_velocity_field);
  } break;

  case GLFW_KEY_W: { // 'z'
    float minC = colorTable.getMin();
    float maxC = colorTable.getMax();
    colorTable.setMinMax(minC * 0.9, maxC * 0.9);
  } break;

  case GLFW_KEY_X: {
    if (mods == GLFW_MOD_SHIFT) {
      do {
        char filename[256];
        snprintf(filename, 256, "screenshot%d.png", Conf.iconf);
        std::cout << MFLOE_INFO << filename << " saved" << std::endl;
        updateTextLine();
        display(window);
        captureScreenshot(filename);
      } while (try_to_readConf(confNum + 1, Conf, confNum) == true);

      /*
ffmpeg -start_number 0 -framerate 10 -i 'screenshot%d.png' -c:v libx264 -pix_fmt yuv420p -crf 23 -preset slow
output.mp4
      */
    } else {
      std::cout << MFLOE_INFO << "screenshot.png saved" << std::endl;
      captureScreenshot("screenshot.png");
    }
  } break;

  case GLFW_KEY_Z: { // 'w/W'
    update_color_bar = 1 - update_color_bar;
    colorBar.setLock(!static_cast<bool>(update_color_bar));
    SEE_SHOW(update_color_bar);
    textZone.addLine("update_color_bar = %d", update_color_bar);
  } break;

  case GLFW_KEY_UP: {
    textZone.increase_nbLine();
  } break;
  case GLFW_KEY_DOWN: {
    textZone.decrease_nbLine();
  } break;

  case GLFW_KEY_LEFT: {
    if (mods == GLFW_MOD_SHIFT) {
      if (try_to_readConf(0, Conf, confNum)) updateTextLine();
    } else {
      if (confNum > 0) {
        if (try_to_readConf(confNum - 1, Conf, confNum)) updateTextLine();
      }
    }
  } break;
  case GLFW_KEY_RIGHT: {
    if (try_to_readConf(confNum + 1, Conf, confNum)) updateTextLine();
  } break;

  case GLFW_KEY_SLASH: { // '='
    fit_view(window);
    textZone.addLine("fitted view");
  } break;
  };

  needsRedraw = true;
}

void mouse_button(GLFWwindow *window, int button, int action, int mods) {
  double x, y;
  glfwGetCursorPos(window, &x, &y);

  if (action == GLFW_RELEASE) {
    mouse_mode = NOTHING;
  } else if (action == GLFW_PRESS) {
    mouse_start[0] = static_cast<int>(x);
    mouse_start[1] = static_cast<int>(y);

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
      if (mods == GLFW_MOD_SHIFT) {
        mouse_mode = PAN;
      } else {
        mouse_mode = ROTATION;
      }
    } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
      mouse_mode = ZOOM;
    }
  }

  needsRedraw = true;
}

void cursor_pos(GLFWwindow *window, double xpos, double ypos) {
  if (mouse_mode == NOTHING) { return; }

  double dx = (xpos - mouse_start[0]) / static_cast<double>(width);
  double dy = (ypos - mouse_start[1]) / static_cast<double>(height);

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
    double Lx  = worldBox.max.x - worldBox.min.x;
    double Ly  = worldBox.max.y - worldBox.min.y;
    double L   = (Lx + Ly);
    double ddx = L * dx;
    double ddy = L * dy;

    worldBox.min.x -= ddx;
    worldBox.max.x -= ddx;
    worldBox.min.y += ddy;
    worldBox.max.y += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = static_cast<int>(xpos);
  mouse_start[1] = static_cast<int>(ypos);

  reshape(window, width, height); // ???
  needsRedraw = true;
}

void display(GLFWwindow *window) {
  glTools::clearBackground(show_background, bottom_color.x, bottom_color.y, bottom_color.z, top_color.x,
                           top_color.y, top_color.z);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (show_breakage_mode) { drawBreakageModes(); }
  if (show_velocity_field == 1) { drawVelocityField(); }
  if (show_forces) { drawForces(); }
  if (show_bonding == 1) { drawBonding(); }
  if (show_particles == 1) { drawParticles(); }
  if (show_biax_system == 1) { drawBiaxSystem(); }

  textZone.draw();
  if (particle_color_option > COLOR_NONE) { colorBar.show(width, height, colorTable); }

  glFlush();
  glfwSwapBuffers(window);
}

void fit_view(GLFWwindow *window) {
  worldBox.min.x = 1e12;
  worldBox.max.x = -1e12;
  worldBox.min.y = 1e12;
  worldBox.max.y = -1e12;
  size_t nDriven = Conf.Drivings.size();
  for (size_t i = nDriven; i < Conf.FloeElements.size(); ++i) {
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

  reshape(window, width, height);
}

void reshape(GLFWwindow *window, int w, int h) {
  if (window) { glfwGetFramebufferSize(window, &w, &h); }

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
  needsRedraw = true;
}

void setColorOption(int option) {
  needsRedraw = true;

  particle_color_option = option;
  color_values.resize(Conf.FloeElements.size(), 0.0);

  switch (particle_color_option) {
  case COLOR_NONE: {
    // nothing to do
  } break;

  case COLOR_RADIUS: {
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) { color_values[i] = Conf.FloeElements[i].radius; }
    if (update_color_bar == 1) colorTable.setMinMax(Rmin, Rmax);
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
    if (update_color_bar == 1) colorTable.setMinMax(Vmin, Vmax);
    colorBar.setTitle("Velocity magnitude");
  } break;

  case COLOR_PRESSURE: {
    computeStressMatrices();

    double Pmin = 0.0;
    double Pmax = 0.0;
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
      double P        = 0.5 * Sigma[i].trace();
      color_values[i] = P;
      if (P < Pmin) { Pmin = P; }
      if (P > Pmax) { Pmax = P; }
    }
    double Pabs = std::max(fabs(Pmin), fabs(Pmax));
    if (update_color_bar == 1) colorTable.setMinMax(-Pabs, Pabs);
    colorBar.setTitle("Pressure");

  } break;

  case COLOR_DAMAGE: {
    if (NbBondsRef.size() != Conf.FloeElements.size()) {
      std::cout << MFLOE_WARN << "NbBondsRef.size() != Conf.FloeElements.size()" << std::endl;
      break;
    }

    std::vector<int> NbBonds(Conf.FloeElements.size(), 0);
    for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
      if (Conf.Interactions[k].isBonded == false) continue;
      NbBonds[Conf.Interactions[k].i] += 1;
      NbBonds[Conf.Interactions[k].j] += 1;
    }

    double Dmin = 0.0;
    double Dmax = 0.0;
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
      double D        = 1.0 - static_cast<double>(NbBonds[i]) / static_cast<double>(NbBondsRef[i]);
      color_values[i] = D;
      if (D < Dmin) { Dmin = D; }
      if (D > Dmax) { Dmax = D; }
    }
    if (update_color_bar == 1) colorTable.setMinMax(Dmin, Dmax);
    else colorTable.setMinMax(0.0f, 1.0f);
    colorBar.setTitle("Damage");
  } break;

  default:
    needsRedraw = false;
    break;
  }
}

void setColor(int i, GLfloat alpha) {
  colorRGBA col;
  if (particle_color_option == COLOR_NONE) {
    col.r = 157;
    col.g = 230;
    col.b = 245;
  } else {
    colorTable.getRGB(color_values[i], &col);
  }
  glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, alpha);
}

void drawWorldBox() {
  glLineWidth(1.0f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  glBegin(GL_LINE_LOOP);
  glVertex2f(Conf.aabb.min.x, Conf.aabb.min.y);
  glVertex2f(Conf.aabb.max.x, Conf.aabb.min.y);
  glVertex2f(Conf.aabb.max.x, Conf.aabb.max.y);
  glVertex2f(Conf.aabb.min.x, Conf.aabb.max.y);
  glEnd();
}

void drawBiaxSystem() {
  if (Conf.biaxSystem == nullptr) return;

  glLineWidth(2.0f);
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);

  glBegin(GL_LINES);
  glVertex2f(Conf.aabb.min.x, Conf.biaxSystem->t_ypos);
  glVertex2f(Conf.aabb.max.x, Conf.biaxSystem->t_ypos);

  glVertex2f(Conf.aabb.min.x, Conf.biaxSystem->b_ypos);
  glVertex2f(Conf.aabb.max.x, Conf.biaxSystem->b_ypos);

  for (size_t i = 0; i < Conf.biaxSystem->l_xpos.size(); i++) {
    glVertex2f(Conf.biaxSystem->l_xpos[i],
               Conf.biaxSystem->b_ypos + (i + 0.05) * Conf.biaxSystem->segmentLengths);
    glVertex2f(Conf.biaxSystem->l_xpos[i],
               Conf.biaxSystem->b_ypos + (i + 0.95) * Conf.biaxSystem->segmentLengths);
  }

  for (size_t i = 0; i < Conf.biaxSystem->r_xpos.size(); i++) {
    glVertex2f(Conf.biaxSystem->r_xpos[i],
               Conf.biaxSystem->b_ypos + (i + 0.05) * Conf.biaxSystem->segmentLengths);
    glVertex2f(Conf.biaxSystem->r_xpos[i],
               Conf.biaxSystem->b_ypos + (i + 0.95) * Conf.biaxSystem->segmentLengths);
  }
  glEnd();
}

void drawParticles() {
  if (mouse_mode != NOTHING) {
    drawWorldBox();

    glLineWidth(1.0f);
    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
    for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
      vec2r pos = Conf.FloeElements[i].pos;
      double R  = Conf.FloeElements[i].radius;
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.15 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();
    }

    return;
  }

  glLineWidth(1.0f);
  size_t nDriven = Conf.Drivings.size();

  // OpenGL trace ce qui est avant en dernier
  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
    vec2r pos = Conf.FloeElements[i].pos;
    double R  = Conf.FloeElements[i].radius;

    glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
    if (show_orientations) {
      double rot = Conf.FloeElements[i].rot;
      glBegin(GL_LINES);
      glVertex2f(pos.x, pos.y);
      glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
      glEnd();
    }

    // contour
    if (contour_particles) {
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();
    }

    // interieur
    if (i >= nDriven) setColor(i, alpha_particles);
    else glColor4f(0.8f, 0.8f, 0.8f, 1.0f);
    glBegin(GL_POLYGON);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();
  }
}

void drawVelocityField() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(arrow_line_width);
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);

  size_t i0    = Conf.Drivings.size();
  double v2max = Conf.FloeElements[i0].vel * Conf.FloeElements[i0].vel;
  for (size_t i = i0; i < Conf.FloeElements.size(); ++i) {
    double v2 = Conf.FloeElements[i].vel * Conf.FloeElements[i].vel;
    if (v2 > v2max) { v2max = v2; }
  }
  double vmax     = sqrt(v2max);
  double factor   = arrow_scale * 2.0 * Rmean / vmax;
  double headSize = 0.5 * Rmean;

  for (size_t i = i0; i < Conf.FloeElements.size(); ++i) {
    vec2r pos = Conf.FloeElements[i].pos;
    vec2r vel = Conf.FloeElements[i].vel;
    vec2r posHead(pos.x + factor * vel.x, pos.y + factor * vel.y);
    vec2r n = posHead - pos;
    n.normalize();
    vec2r t(-n.y * arrow_head_ratio, n.x * arrow_head_ratio);

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

void drawBreakageModes() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(1.0f);
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);

  for (size_t k = 0; k < breakageDataInfo.size(); k++) {
    if (breakageDataInfo[k].t >= Conf.t) break;

    double fnb    = breakageDataInfo[k].fnb;
    double fn0Neg = sqrt(2.0 * breakageDataInfo[k].coverage * breakageDataInfo[k].A * Conf.kn * Conf.Gc);
    double fn0Pos = fn0Neg * Conf.GcComprRatio;

    switch (Conf.yieldSurfaceModel) {
    case YIELD_PARABOLA: {
      fnb /= fn0Neg;
    } break;
    case YIELD_HEAVISIDE: {
      fnb /= fn0Neg;
    } break;
    case YIELD_ELLIPSE: {
      fnb /= fn0Neg;
    } break;
    case YIELD_ELLIPSE_ASYM: {
      if (fnb < 0.0) {
        fnb /= fn0Neg;
      } else {
        fnb /= fn0Pos;
      }
    } break;
    }

    // colorRGBA col;
    // BreakageColors.getRGB(-fnb, &col);
    // glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);

    if (fnb < 0.0) {
      if (fnb < -0.5) glColor3f(1.0f, 0.0f, 0.0f); // rouge
      else glColor3f(0.82f, 0.44f, 0.8f);          // rose
    } else {
      if (fnb > 0.5) glColor3f(0.0f, 0.0f, 1.0f); // bleu
      else glColor3f(0.29f, 0.95f, 0.88f);        // cyan
    }

    vec2r posi = Conf.FloeElements[breakageDataInfo[k].i].pos;
    vec2r posj = Conf.FloeElements[breakageDataInfo[k].j].pos;
    double Ri  = Conf.FloeElements[breakageDataInfo[k].i].radius;
    double Rj  = Conf.FloeElements[breakageDataInfo[k].j].radius;

    vec2r n    = posj - posi;
    double len = n.normalize();
    double dn  = len - Ri - Rj;

    vec2r pos = posi + (Ri + 0.5 * dn) * n;

    glBegin(GL_POLYGON);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.1 * M_PI) {
      glVertex2f(pos.x + Rmean * cos(angle), pos.y + Rmean * sin(angle));
    }
    glEnd();
  }
}

void drawBonding() {
  if (mouse_mode != NOTHING) { return; }

  glLineWidth(1.0f);
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);

  glBegin(GL_QUADS);
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    if (Conf.Interactions[k].isBonded == false) continue;

    size_t i        = Conf.Interactions[k].i;
    size_t j        = Conf.Interactions[k].j;
    vec2r posi      = Conf.FloeElements[i].pos;
    vec2r posj      = Conf.FloeElements[j].pos;
    double Di       = Conf.FloeElements[i].radius /*+ 0.5 * Conf.Interactions[k].dn0*/;
    double Dj       = Conf.FloeElements[j].radius /*+ 0.5 * Conf.Interactions[k].dn0*/;
    double hij      = 0.5 * (Conf.FloeElements[i].height + Conf.FloeElements[j].height);
    double half_lij = 0.5 * Conf.Interactions[k].A / hij;

    vec2r nij = posj - posi;
    nij.normalize();
    vec2r tij(-nij.y, nij.x);

    vec2r ci1 = posi + Di * nij - half_lij * tij;
    vec2r ci2 = posi + Di * nij + half_lij * tij;
    // setColor(i, alpha_particles);
    glVertex2f(ci1.x, ci1.y);
    glVertex2f(ci2.x, ci2.y);

    vec2r cj1 = posj - Dj * nij - half_lij * tij;
    vec2r cj2 = posj - Dj * nij + half_lij * tij;
    // setColor(j, alpha_particles);
    glVertex2f(cj2.x, cj2.y);
    glVertex2f(cj1.x, cj1.y);
  }
  glEnd();
}

void drawForces() {
  if (mouse_mode != NOTHING) { return; }

  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i   = Conf.Interactions[k].i;
    size_t j   = Conf.Interactions[k].j;
    vec2r posi = Conf.FloeElements[i].pos;
    vec2r posj = Conf.FloeElements[j].pos;

    GLfloat width = Rmean * (fabs(Conf.Interactions[k].fn) / fnMax);

    if (Conf.Interactions[k].fn > 0.0) {
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    } else {
      glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
    }

    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }
}

void computeStressMatrices() {
  Sigma.clear();
  Sigma.resize(Conf.FloeElements.size());
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i     = Conf.Interactions[k].i;
    size_t j     = Conf.Interactions[k].j;
    vec2r unit_n = Conf.FloeElements[j].pos - Conf.FloeElements[i].pos;
    double l     = unit_n.normalize();
    double dn    = l - Conf.FloeElements[i].radius - Conf.FloeElements[j].radius;
    vec2r Bi     = (Conf.FloeElements[i].radius + 0.5 * dn) * unit_n;
    vec2r Bj     = -(Conf.FloeElements[j].radius + 0.5 * dn) * unit_n;
    vec2r unit_t(-unit_n.y, unit_n.x);
    vec2r fji = (Conf.Interactions[k].fn + Conf.Interactions[k].fnb) * unit_n +
                (Conf.Interactions[k].ft + Conf.Interactions[k].ftb) * unit_t;
    vec2r fij = -fji;

    Sigma[i].xx += fji.x * Bi.x;
    Sigma[i].xy += fji.x * Bi.y;
    Sigma[i].yx += fji.y * Bi.x;
    Sigma[i].yy += fji.y * Bi.y;

    Sigma[j].xx += fij.x * Bj.x;
    Sigma[j].xy += fij.x * Bj.y;
    Sigma[j].yx += fij.y * Bj.x;
    Sigma[j].yy += fij.y * Bj.y;
  }

  for (size_t i = 0; i < Conf.FloeElements.size(); ++i) {
    double meanxy = 0.5 * (Sigma[i].yx + Sigma[i].xy);
    Sigma[i].yx = Sigma[i].xy = meanxy;
    double V                  = M_PI * Conf.FloeElements[i].radius * Conf.FloeElements[i].radius;
    Sigma[i] *= (1.0 / V);
  }
}

void getReferenceNumberOfBonds() {
  NbBondsRef.clear();
  NbBondsRef.resize(Conf.FloeElements.size(), 0);
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    if (Conf.Interactions[k].isBonded == false) continue;
    NbBondsRef[Conf.Interactions[k].i] += 1;
    NbBondsRef[Conf.Interactions[k].j] += 1;
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

void updateTextLine() {
  textZone.addLine("#conf %d,  time %0.4g s", confNum, Conf.t);
}

int main(int argc, char *argv[]) {
  INIT_TIMERS();

  printInfo();
  Conf.debugEnabled = false; // disable debug logs even when compiled with 'DEBUG' defined

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

  // the number of bonds per element in the first opened file is the reference number
  getReferenceNumberOfBonds();

  if (fileTool::fileExists("breakage_data_file.txt")) {
    std::cout << "Lecture de breakage_data_file.txt" << std::endl;
    std::ifstream file("breakage_data_file.txt");
    std::string trash;
    getline(file, trash);
    breakageDataInfo.clear();
    BreakageDataInfo DI;
    while (file.good()) {
      file >> DI.t >> DI.i >> DI.j >> DI.fnb >> DI.ftb >> DI.fsb >> DI.coverage >> DI.A;
      breakageDataInfo.push_back(DI);
    }
    std::cout << "breakageDataInfo.size() = " << breakageDataInfo.size() << std::endl;
  }

  mouse_mode = NOTHING;

  readTomlOptions();

  if (!glfwInit()) {
    fprintf(stderr, "Failed to initialize GLFW\n");
    return -1;
  }

  // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

  // Désactive le support Retina (force une fenêtre en résolution 1x)
  glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);

  GLFWwindow *window = glfwCreateWindow(width, height, "CONF VISUALIZER (GLFW)", NULL, NULL);
  if (!window) {
    fprintf(stderr, "Failed to create GLFW window\n");
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, keyboard);
  glfwSetMouseButtonCallback(window, mouse_button);
  glfwSetCursorPosCallback(window, cursor_pos);
  glfwSetFramebufferSizeCallback(window, reshape);

  int X0 = (glfwGetVideoMode(glfwGetPrimaryMonitor())->width - width) / 2;
  int Y0 = (glfwGetVideoMode(glfwGetPrimaryMonitor())->height - height) / 2;
  glfwSetWindowPos(window, X0, Y0);

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

  glEnable(GL_DEPTH_TEST);

  setColorOption(particle_color_option);
  if (fit_at_loading != 0) fit_view(window);

  while (!glfwWindowShouldClose(window)) {
    glfwWaitEvents();

    if (needsRedraw) {
      reshape(window, width, height);
      display(window);
      needsRedraw = false;
    }
  }

  glfwTerminate();
  return 0;
}
