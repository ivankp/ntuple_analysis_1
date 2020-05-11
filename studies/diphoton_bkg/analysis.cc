#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#elif defined(ANALYSIS_INIT) // =====================================

h_(pT_yy)

h_(m_yy)
h_(m_yy,pT_yy)

h_(m_yy_zoom)
h_(m_yy_zoom,pT_yy)

#elif defined(ANALYSIS_LOOP) // =====================================

const double m_yy = higgs.M();
const double pT_yy = higgs.Pt();

FILL(pT_yy)

FILL(m_yy)
FILL(m_yy,pT_yy)

if (100 <= m_yy && m_yy < 200) {
  h_m_yy_zoom(m_yy);
  h_m_yy_zoom_pT_yy(m_yy,pT_yy);
}

#endif
