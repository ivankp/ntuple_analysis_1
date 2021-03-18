#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#define CATEGORIES (isp)

#elif defined(ANALYSIS_INIT) // =====================================

h_(Hj1_mass)
h_(Hj1_mass,H_pT)
h_(Hj1_mass,j1_pT)
h_(Hj1_mass,j2_pT)
h_(Hj1_mass,H_pT,j1_pT)
h_(Hj1_mass,H_pT,j2_pT)
h_(Hj1_mass,j1_pT,j2_pT)

h_(Hj2_mass)
h_(Hj2_mass,H_pT)
h_(Hj2_mass,j1_pT)
h_(Hj2_mass,j2_pT)
h_(Hj2_mass,H_pT,j1_pT)
h_(Hj2_mass,H_pT,j2_pT)
h_(Hj2_mass,j1_pT,j2_pT)

h_(jj_mass)
h_(jj_mass,H_pT)
h_(jj_mass,j1_pT)
h_(jj_mass,j2_pT)
h_(jj_mass,H_pT,j1_pT)
h_(jj_mass,H_pT,j2_pT)
h_(jj_mass,j1_pT,j2_pT)

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < 1) continue;

const double H_pT = higgs.Pt();
const double j1_pT = jets[0].pt();

const auto   Hj1 = higgs + jets[0];
const double Hj1_mass = Hj1.M();

FILL(Hj1_mass)
FILL(Hj1_mass,H_pT)
FILL(Hj1_mass,j1_pT)
FILL(Hj1_mass,H_pT,j1_pT)

if (njets < 2) continue;

const double j2_pT = jets[1].pt();

const auto   Hj2 = higgs + jets[1];
const double Hj2_mass = Hj2.M();

const auto   jj = jets[0] + jets[1];
const double jj_mass = jj.m();

FILL(Hj1_mass,j2_pT)
FILL(Hj1_mass,H_pT,j2_pT)
FILL(Hj1_mass,j1_pT,j2_pT)

FILL(Hj2_mass)
FILL(Hj2_mass,H_pT)
FILL(Hj2_mass,j1_pT)
FILL(Hj2_mass,j2_pT)
FILL(Hj2_mass,H_pT,j1_pT)
FILL(Hj2_mass,H_pT,j2_pT)
FILL(Hj2_mass,j1_pT,j2_pT)

FILL(jj_mass)
FILL(jj_mass,H_pT)
FILL(jj_mass,j1_pT)
FILL(jj_mass,j2_pT)
FILL(jj_mass,H_pT,j1_pT)
FILL(jj_mass,H_pT,j2_pT)
FILL(jj_mass,j1_pT,j2_pT)

#endif
