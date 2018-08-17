#ifndef HIST_HJ
#define HIST_HJ test.cc
#include "hist_Hjets.hh"
#endif

#ifdef HIST_HJ_GLOBAL // ===============================================

#endif
#ifdef HIST_HJ_INIT // =================================================

h_(AA_dR)

#endif
#ifdef HIST_HJ_LOOP // =================================================

const double HM = higgs->M();
if (HM<121) continue;
if (HM>129) continue;

h_AA_dR(A1->DeltaR(*A2));

#endif

