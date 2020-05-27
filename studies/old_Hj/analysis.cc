#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT)
h_(j1_pT)
h_(j2_pT)

h_(Hj_mass)
h_(Hj_pT)

h_(Hjj_mass)
h_(Hjj_pT)

#elif defined(ANALYSIS_LOOP) // =====================================

h_H_pT(higgs.Pt());

if (njets < 1) continue; // -----------------------------------------

h_j1_pT(jets[0].pt());

const auto Hj = higgs + jets[0];

h_Hj_mass(Hj.M());
h_Hj_pT(Hj.Pt());

if (njets < 2) continue; // -----------------------------------------

h_j2_pT(jets[1].pt());

const auto Hjj = Hj + jets[1];

h_Hjj_mass(Hjj.M());
h_Hjj_pT(Hjj.Pt());

#endif
