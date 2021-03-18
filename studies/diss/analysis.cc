#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

// #define CATEGORIES (photon_cuts)

constexpr double pi = M_PI;
constexpr double twopi = M_PI*2;

double dphi_mpi_pi(double dphi) noexcept {
  while (dphi >= pi) dphi -= twopi;
  while (dphi < -pi) dphi += twopi;
  return dphi;
}
double abs_dphi(double dphi) noexcept {
  return std::abs(dphi_mpi_pi(dphi));
}

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT) h_(j1_pT) h_(j2_pT) h_(j3_pT)

h_(H_eta) h_(j1_eta) h_(j2_eta) h_(j3_eta)
h_(H_rap)

h_(Hj_m)
h_(jj_m)

h_(jj_drap)
h_(jj_dphi)

#elif defined(ANALYSIS_LOOP) // =====================================

h_H_pT(higgs.Pt());
h_H_eta(higgs.Eta());
h_H_rap(higgs.Rapidity());

if (njets < 1) continue; // -----------------------------------------

h_j1_pT(jets[0].pt());
h_j1_eta(jets[0].eta());

h_Hj_m((higgs+jets[1]).M());

if (njets < 2) continue; // -----------------------------------------

h_j2_pT(jets[1].pt());
h_j2_eta(jets[1].eta());

h_jj_m((jets[0]+jets[1]).m());

h_jj_drap(std::abs(jets[0].rap()-jets[1].rap()));
h_jj_dphi(abs_dphi(jets[0].phi()-jets[1].phi())/pi);

if (njets < 3) continue; // -----------------------------------------

h_j3_pT(jets[2].pt());
h_j3_eta(jets[2].eta());

#endif
