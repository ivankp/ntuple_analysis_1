#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#define CATEGORIES (photon_cuts)

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT) h_(j1_pT)

#elif defined(ANALYSIS_LOOP) // =====================================

h_H_pT(higgs.Pt());
h_j1_pT(jets[0].pt());

#endif
