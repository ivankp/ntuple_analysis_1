#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT) h_(H_y)

h_(j1_pT) h_(j1_y)
h_(j2_pT) h_(j2_y)
h_(j3_pT) h_(j3_y)

h_(Hj_mass) h_(Hj_pT) h_(Hj_y)

h_(Hjj_mass) h_(Hjj_pT) h_(Hjj_y)

h_(jj_mass) h_(jj_pT) h_(jj_y)

#elif defined(ANALYSIS_LOOP) // =====================================

h_H_pT(higgs.Pt());
h_H_y(higgs.Rapidity());

if (njets < 1) continue; // -----------------------------------------

h_j1_pT(jets[0].pt());
h_j1_pT(jets[0].rap());

const auto Hj = higgs + jets[0];

h_Hj_mass(Hj.M());
h_Hj_pT(Hj.Pt());
h_Hj_y(Hj.Rapidity());

if (njets < 2) continue; // -----------------------------------------

h_j2_pT(jets[1].pt());
h_j2_y(jets[1].rap());

const auto Hjj = Hj + jets[1];

h_Hjj_mass(Hjj.M());
h_Hjj_pT(Hjj.Pt());
h_Hjj_y(Hjj.Rapidity());

const auto jj = jets[0] + jets[1];

h_jj_mass(jj.m());
h_jj_pT(jj.pt());
h_jj_y(jj.rap());

if (njets < 3) continue; // -----------------------------------------

h_j3_pT(jets[2].pt());
h_j3_y(jets[2].rap());

#endif
