#include "Interaction.hpp"

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

  meff = I.meff;
  //kn   = I.kn;
  //kt   = I.kt;
  //mu   = I.mu;
  //muR  = I.muR;
  //fadh = I.fadh;

  //damp = I.damp;
  A = I.A;
  coverage = I.coverage;
  //Gc   = I.Gc;
  dn0  = I.dn0;
  t0  = I.t0;
}