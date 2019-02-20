#ifndef ANALYSIS

#define OUTPUT_SCRIBE

#define HIST_HJ "hist_Hj_rel.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(VBF_cuts,(all)(with_VBF_cuts))
MAKE_ENUM(Nj,(all)(njets_eq)(njets_gt))

#define CATEGORIES (photon_cuts)(VBF_cuts)(Nj)(isp)

TLorentzVector root(const fj::PseudoJet& j) {
  return { j[0], j[1], j[2], j[3] };
}

#elif defined(ANALYSIS_INIT) // =====================================

#define HIST_MAX_D 1

auto h_reserve = [&](unsigned n, auto name_f) {
  auto hs = reserve<hist<1>>(n);
  for (unsigned i=0; i<n; ++i) {
    const auto name = name_f(i);
    hs.emplace_back(name,ra[name]);
  }
  return hs;
};

h_(AA_pTrat)
h_(AA_dpT)
h_(AA_dR)
h_(AA_dy)
h_(AA_dphi)

#define h_Hj_(X) auto h_Hj_##X = h_reserve(njets_born+1, \
  [](unsigned i){ return cat("H_j",i+1,"_"#X); } );

h_Hj_(pTrat)
h_Hj_(dpT)
h_Hj_(dR)
h_Hj_(dy)
h_Hj_(dphi)
h_Hj_(pT)
h_Hj_(mass)

#define h_jj_(X) h_(j1_j2_##X) h_(jjfb_##X)

h_jj_(pTrat)
h_jj_(dpT)
h_jj_(dR)
h_jj_(dy)
h_jj_(dphi)
h_jj_(pT)
h_jj_(mass)

#define h_Hjj_(X) h_(H_j1j2_##X) h_(H_jjfb_##X)

h_Hjj_(pTrat)
h_Hjj_(dpT)
h_Hjj_(dR)
h_Hjj_(dy)
h_Hjj_(dphi)
h_Hjj_(pT)
h_Hjj_(mass)

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < njets_born) continue;
const unsigned max_njets = std::min(njets,njets_born+1);

bin_t::id<VBF_cuts>(0);
bin_t::id<Nj>(njets > njets_born);

const double A1_phi = photons[0].Phi();
const double A2_phi = photons[1].Phi();
const double A1_y   = photons[0].Rapidity();
const double A2_y   = photons[1].Rapidity();

h_AA_pTrat(A1_pT/A2_pT);
h_AA_dpT  (A1_pT-A1_pT);
h_AA_dR   (deltaR(A1_eta,A2_eta,A1_phi,A2_phi));
h_AA_dy   (std::abs(A1_y-A2_y));
h_AA_dphi (dphi(A1_phi,A2_phi));

const double H_pT  = higgs.Pt();
const double H_y   = higgs.Rapidity();
const double H_eta = higgs.Eta();
const double H_phi = higgs.Phi();

const auto jet_vars = jets | [](const auto& jet){
  struct vars { double pT, y, eta, phi; };
  return vars { jet.pt(), jet.rap(), jet.eta(), jet.phi_std() };
};

for (unsigned j=0; j<max_njets; ++j) {
  const auto& jet = jet_vars[j];
  h_Hj_pTrat[j](H_pT/jet.pT);
  h_Hj_dpT  [j](H_pT-jet.pT);
  h_Hj_dR   [j](deltaR(H_eta,jet.eta,H_phi,jet.phi));
  h_Hj_dy   [j](std::abs(H_y-jet.y));
  h_Hj_dphi [j](dphi(H_phi,jet.phi));
  const auto Hj = higgs + jets[j];
  h_Hj_pT   [j](Hj.Pt());
  h_Hj_mass [j](Hj.M());
}

if (njets < 2) continue;

std::array<decltype(&jet_vars[0]),2> _jj { &jet_vars[0], &jet_vars[1] };
double jj_dy = std::abs(_jj[0]->y - _jj[1]->y);
auto jj = jets[0]+jets[1];
double jj_mass = jj.m();
double jj_pT = jj.pt();

bin_t::id<VBF_cuts>(
  (jj_dy > 2.8) and
  (jj_mass > 400)
);

h_j1_j2_pTrat(_jj[0]->pT/_jj[1]->pT);
h_j1_j2_dpT  (_jj[0]->pT-_jj[1]->pT);
h_j1_j2_dR   (deltaR(_jj[0]->eta,_jj[1]->eta,_jj[0]->phi,_jj[1]->phi));
h_j1_j2_dy   (std::abs(_jj[0]->y - _jj[1]->y));
h_j1_j2_dphi (dphi(_jj[0]->phi,_jj[1]->phi));
h_j1_j2_pT   (jj_pT);
h_j1_j2_mass (jj_mass);

h_H_j1j2_pTrat(H_pT/jj_pT);
h_H_j1j2_dpT  (H_pT-jj_pT);
h_H_j1j2_dR   (higgs.DeltaR(root(jj)));
h_H_j1j2_dy   (std::abs(H_y - jj.rap()));
h_H_j1j2_dphi (dphi(jj.phi_std(),H_phi));
h_H_j1j2_pT   ((higgs+jj).Pt());
h_H_j1j2_mass ((higgs+jj).M());

if (njets > 2) {
  std::array<unsigned,2> fb{0,0};
  for (unsigned i=1; i<max_njets; ++i) {
    const auto& j = jet_vars[i];
    const auto& f = jet_vars[fb[0]];
    const auto& b = jet_vars[fb[1]];
    if (j.eta < b.eta) fb[1] = i;
    if (j.eta > f.eta) fb[0] = i;
  }
  if (fb[0] > fb[1]) std::swap(fb[0],fb[1]);
  if (fb[0]!=0 || fb[1]!=1) {
    jj = jets[fb[0]] + jets[fb[1]];
    _jj = { &jet_vars[fb[0]], &jet_vars[fb[1]] };

    jj_dy = std::abs(_jj[0]->y - _jj[1]->y);
    jj = jets[0]+jets[1];
    jj_mass = jj.m();

    bin_t::id<VBF_cuts>(
      (jj_dy > 2.8) and
      (jj_mass > 400)
    );
  }
}

h_jjfb_pTrat(_jj[0]->pT/_jj[1]->pT);
h_jjfb_dpT  (_jj[0]->pT-_jj[1]->pT);
h_jjfb_dR   (deltaR(_jj[0]->eta,_jj[1]->eta,_jj[0]->phi,_jj[1]->phi));
h_jjfb_dy   (std::abs(_jj[0]->y - _jj[1]->y));
h_jjfb_dphi (dphi(_jj[0]->phi,_jj[1]->phi));
h_jjfb_pT   (jj_pT);
h_jjfb_mass (jj_mass);

h_H_jjfb_dpT  (H_pT-jj_pT);
h_H_jjfb_dR   (higgs.DeltaR(root(jj)));
h_H_jjfb_dy   (std::abs(H_y - jj.rap()));
h_H_jjfb_dphi (dphi(jj.phi_std(),H_phi));
h_H_jjfb_pT   ((higgs+jj).Pt());
h_H_jjfb_mass ((higgs+jj).M());

#elif defined(ANALYSIS_END) // ======================================

#endif

