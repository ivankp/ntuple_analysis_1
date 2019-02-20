#ifndef ANALYSIS

#define OUTPUT_SCRIBE

#define HIST_HJ "hist_Hjets.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#elif defined(ANALYSIS_INIT) // =====================================

#define HIST_MAX_D 1

#define h_A_(X) std::array<hist<1>,2> h_A_##X {{ \
  { "A1_"#X, ra["A1_"#X] }, { "A2_"#X, ra["A2_"#X] } }};

auto h_reserve = [&](unsigned n, auto name_f) {
  auto hs = reserve<hist<1>>(n);
  for (unsigned i=0; i<n; ++i) {
    const auto name = name_f(i);
    hs.emplace_back(name,ra[name]);
  }
  return hs;
};

#define h_j_(X) auto h_jet_##X = h_reserve(njets_born+1, \
  [](unsigned i){ return cat("jet",i+1,"_"#X); } );

h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass) h_(H_mass_hgam)

h_j_(pT) h_j_(y) h_j_(eta) h_j_(phi) h_j_(mass) h_j_(tau)

h_(jets_tau_max) h_(jets_tau_sum)

h_A_(pT) h_A_(y) h_A_(eta) h_A_(phi)

h_(AA_cosTS_Hframe) h_(AA_cosTS_CSframe)
h_(AA_pTt)

auto h_Hjs_mass = h_reserve(njets_born+1,
  [](unsigned i){ return cat('H',i+1,"j_mass"); } );

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < njets_born) continue;
const unsigned max_njets = std::min(njets,njets_born+1);

const double H_pT   = higgs.Pt();
const double H_y    = higgs.Rapidity();
const double H_mass = higgs.M();
h_H_pT(H_pT);
h_H_y(H_y);
h_H_eta(higgs.Eta());
h_H_phi(higgs.Phi());
h_H_mass(H_mass);
h_H_mass_hgam(H_mass);

double HT = H_pT;
double jets_tau_max = 0, jets_tau_sum = 0;
auto Hjs = higgs;

for (unsigned j=0; j<max_njets; ++j) {
  const auto jet_pT = jets[j].pt();
  HT += jet_pT;
  h_jet_pT  [j](jet_pT);
  h_jet_y   [j](jets[j].rapidity());
  h_jet_eta [j](jets[j].eta());
  h_jet_phi [j](jets[j].phi_std());
  h_jet_mass[j](jets[j].m());

  const auto jet_tau = tau(jets[j],H_y);
  if (jet_tau > jets_tau_max) jets_tau_max = jet_tau;
  jets_tau_sum += jet_tau;
  h_jet_tau[j](jet_tau);

  Hjs += jets[j];
  h_Hjs_mass[j](Hjs.M());
}

h_HT(HT);

h_jets_tau_max(jets_tau_max);
h_jets_tau_sum(jets_tau_sum);

h_A_pT [0](A1_pT ); h_A_pT [1](A2_pT );
h_A_eta[0](A1_eta); h_A_eta[1](A2_eta);
for (unsigned i=0; i<2; ++i) {
  h_A_y  [i](photons[i].Rapidity());
  h_A_phi[i](photons[i].Phi());
}

// cos θ* in Collins-Soper frame
h_AA_cosTS_CSframe(
  std::sinh(std::abs(A1_eta-A2_eta)) * A1_pT * A2_pT * 2
  / ( std::sqrt(1.+sq(H_pT/H_mass)) * sq(H_mass) ) );

// cos θ* in Higgs rest frame
const auto Hframe_boost = -higgs.BoostVector();
auto A1_Hframe = photons[0];
A1_Hframe.Boost(Hframe_boost);
h_AA_cosTS_Hframe( std::abs(A1_Hframe.CosTheta()) );

h_AA_pTt(pTt(photons[0],photons[1]));

#endif

