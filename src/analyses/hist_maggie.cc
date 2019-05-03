#ifndef ANALYSIS

#define HIST_HJ "hist_maggie.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(H_mass_region,(all)(left)(sig)(right))

#define CATEGORIES (photon_cuts)(isp)(H_mass_region)

TLorentzVector root(const fj::PseudoJet& j) {
  return { j[0], j[1], j[2], j[3] };
}

#elif defined(ANALYSIS_INIT) // =====================================

#define HIST_MAX_D 1

/*
auto h_reserve = [&](unsigned n, auto name_f) {
  auto hs = reserve<hist<1>>(n);
  for (unsigned i=0; i<n; ++i) {
    const auto name = name_f(i);
    hs.emplace_back(name,ra[name]);
  }
  return hs;
};
*/

h_(AA_pTrat)
h_(AA_dpT)
h_(AA_dR)
h_(AA_dy)
h_(AA_dphi)

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < njets_min) continue;
// const unsigned max_njets = std::min(njets,njets_min+1);

const double H_mass = higgs.M();

bin_t::id<H_mass_region>((unsigned)(
  H_mass<121 ? H_mass_region::left  : (
  H_mass>129 ? H_mass_region::right : (
  H_mass_region::sig
))));

const double A1_phi = photons[0].Phi();
const double A2_phi = photons[1].Phi();
const double A1_y   = photons[0].Rapidity();
const double A2_y   = photons[1].Rapidity();

h_AA_pTrat(A1_pT/A2_pT);
h_AA_dpT  (A1_pT-A1_pT);
h_AA_dR   (deltaR(A1_eta,A2_eta,A1_phi,A2_phi));
h_AA_dy   (std::abs(A1_y-A2_y));
h_AA_dphi (dphi(A1_phi,A2_phi));

/*
const double H_pT  = higgs.Pt();
const double H_y   = higgs.Rapidity();
const double H_eta = higgs.Eta();
const double H_phi = higgs.Phi();

const auto jet_vars = jets | [](const auto& jet){
  struct vars { double pT, y, eta, phi; };
  return vars { jet.pt(), jet.rap(), jet.eta(), jet.phi_std() };
};
*/

#elif defined(ANALYSIS_END) // ======================================

#endif

