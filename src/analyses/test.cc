#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "test.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#elif defined(ANALYSIS_INIT) // =====================================

// h_(H_pT)
// h_(HT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < njets_born) continue;

// const double H_pT = higgs.Pt();
// h_H_pT(H_pT);

// double HT = H_pT;
// for (const auto& jet : jets) HT += jet.pt();
// h_HT(HT);
//
// h_H_y(higgs.Rapidity());
// h_H_eta(higgs.Eta());
// h_H_phi(higgs.Phi());
// h_H_mass(higgs.M());

#elif defined(ANALYSIS_END) // ======================================

#endif

