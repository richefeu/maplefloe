#include "Interaction.hpp"
#include "MapleFloe.hpp"

Interaction::Interaction() : i(0), j(0) {}

Interaction::Interaction(size_t I, size_t J) : i(I), j(J) {}

void Interaction::copy(Interaction &I) {
  isBonded = I.isBonded;

  fn  = I.fn;
  fnb = I.fnb;
  ft  = I.ft;
  ftb = I.ftb;
  fs  = I.fs;
  fsb = I.fsb;

  // meff     = I.meff;
  A        = I.A;
  coverage = I.coverage;
  dn0      = I.dn0;
  t0       = I.t0;
}

// Nbonds needs to be already updated when using this method
void Interaction::computeBondedArea(double h_i, double radius_i, int Nbonds_i, double h_j, double radius_j,
                                    int Nbonds_j) {
  double h = std::min(h_i, h_j);

  double li = 2.0 * radius_i;
  if (Nbonds_i > 2) { li *= sin(M_PI / Nbonds_i); }
  double lj = 2.0 * radius_j;
  if (Nbonds_j > 2) { lj *= sin(M_PI / Nbonds_j); }

  double l = std::min(li, lj);
  A        = h * l;
}

std::function<double(Interaction &, MFloe &)> Interaction::ruptureCriterion[4]{

    // -----------------------------------------------------------
    // ------ PARABOLA
    [](Interaction &I, MFloe &System) -> double {
      double fcn  = sqrt(2.0 * I.coverage * I.A * System.kn * System.Gc);
      double fct2 = 2.0 * I.coverage * I.A * System.kt * System.Gc;

      return (-I.fnb / fcn + (I.ftb * I.ftb + I.fsb * I.fsb) / fct2 - 1.0);
    },

    // -----------------------------------------------------------
    // ------ HEAVISIDE
    [](Interaction &I, MFloe &System) -> double {      
      double fct2 = 2.0 * I.coverage * I.A * System.kt * System.Gc;
      if (I.fnb < 0.0) {
        return ((I.ftb * I.ftb + I.fsb * I.fsb) / fct2 - 1.0);
      }
      double fcn2 = 2.0 * I.coverage * I.A * System.kn * System.Gc;
      return ((I.fnb * I.fnb) / fcn2 + (I.ftb * I.ftb + I.fsb * I.fsb) / fct2 - 1.0);
    },

    // -----------------------------------------------------------
    // ------ ELLIPSE
    [](Interaction &I, MFloe &System) -> double {
      double fcn2 = 2.0 * I.coverage * I.A * System.kn * System.Gc;
      double fct2 = 2.0 * I.coverage * I.A * System.kt * System.Gc;

      return ((I.fnb * I.fnb) / fcn2 + (I.ftb * I.ftb + I.fsb * I.fsb) / fct2 - 1.0);
    },

    // -----------------------------------------------------------
    // ------ ELLIPSE_ASYM
    [](Interaction &I, MFloe &System) -> double {
      double factMinus = (1.0 - System.GcComprRatio);
      double factPlus  = (1.0 + System.GcComprRatio);

      double fcn  = sqrt(2.0 * I.coverage * I.A * System.kn * System.Gc);
      double fct2 = 2.0 * I.coverage * I.A * System.kt * System.Gc;

      double top_n = 2.0 * I.fnb + fcn * factMinus;
      double bot_n = fcn * factPlus;

      double top_ts = 4.0 * System.GcComprRatio * (I.ftb * I.ftb + I.fsb * I.fsb);
      double bot_ts = fct2 * factPlus * factPlus;

      return ((top_n * top_n) / (bot_n * bot_n) + top_ts / bot_ts - 1.0);
    }};
