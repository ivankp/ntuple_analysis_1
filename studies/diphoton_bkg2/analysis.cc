#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#define CATEGORIES (photon_cuts)

#elif defined(ANALYSIS_INIT) // =====================================

h_(m_yy_105_160,pT_yy)
h_(m_yy_121_129,pT_yy)

#elif defined(ANALYSIS_LOOP) // =====================================

const double m_yy = higgs.M();
const double pT_yy = higgs.Pt();

h_m_yy_105_160_pT_yy(m_yy,pT_yy);
h_m_yy_121_129_pT_yy(m_yy,pT_yy);

#endif
