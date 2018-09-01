#ifndef HIST_HJ
#define HIST_HJ test.cc
#include "hist_Hjets.hh"
#endif

#ifdef HIST_HJ_GLOBAL // ===============================================

#endif
#ifdef HIST_HJ_INIT // =================================================

h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

#endif
#ifdef HIST_HJ_LOOP // =================================================

if (njets < njets_required) continue;

const double H_pT = higgs.Pt();
h_H_pT(H_pT);

double HT = H_pT;
for (const auto& jet : jets) HT += jet.pt();
h_HT(HT);

h_H_y(higgs.Rapidity());
h_H_eta(higgs.Eta());
h_H_phi(higgs.Phi());
h_H_mass(higgs.M());

#endif

