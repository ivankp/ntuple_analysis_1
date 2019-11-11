#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

// #define CATEGORIES (photon_cuts)

#elif defined(ANALYSIS_INIT) // =====================================

h_(yy_pT)
h_(yy_pT_myy_117_133)
h_(yy_pT_myy_121_129)
h_(yy_pT_myy_124_126)
h_(j1_pT)

#elif defined(ANALYSIS_LOOP) // =====================================

const double
  myy = higgs.M(),
  yy_pT = higgs.Pt();

FILL(yy_pT);

if (117.<myy && myy<133.) {
  h_yy_pT_myy_117_133(yy_pT);
  if (121.<myy && myy<129.) {
    h_yy_pT_myy_121_129(yy_pT);
    if (124.<myy && myy<126.) {
      h_yy_pT_myy_124_126(yy_pT);
    }
  }
}

if (njets < 1) continue; // -----------------------------------------

h_j1_pT(jets[0].pt());

#endif
