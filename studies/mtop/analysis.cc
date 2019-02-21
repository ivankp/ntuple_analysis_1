#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

struct ang_t {
  double M, cos_theta;
};
ang_t ang(const TLorentzVector& higgs, const TLorentzVector& jet) {
  const auto Q = higgs + jet;
  const double Q2 = Q*Q;

  const TLorentzVector Z(0,0,Q.E(),Q.Pz());
  const auto ell = ((Q*jet)/Q2)*higgs - ((Q*higgs)/Q2)*jet;

  return {
    sqrt(Q2),
    (ell*Z) / sqrt(sq(ell)*sq(Z))
  };
};

MAKE_ENUM(H_rapidity_cut,(all)(central_higgs))
MAKE_ENUM(fat_jet,(all)(nsubjets_1)(nsubjets_2))
MAKE_ENUM(p1_flavor,(all)(p1_g)(p1_q))
MAKE_ENUM(p2_flavor,(all)(p2_g)(p2_q))

#define CATEGORIES (isp)(photon_cuts)(H_rapidity_cut)(fat_jet)(p1_flavor)(p2_flavor)

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT) h_(j1_pT) h_(j2_pT)

h_(j1_nsub)

h_(H_y)
h_(H_cosTheta)
h_(Hj_mass)

h_(j1_mass)
h_(j2_mass)

h_(j1_x)

h_(H_cosTheta,Hj_mass)

h_(Hj_mass,j1_x)
h_(j1_mass,j1_x)

h_(H_pT,Hj_mass)
h_(j1_pT,Hj_mass)

h_(pp_pTrat,Hj_mass)
h_(pp_dphi,Hj_mass)

#elif defined(ANALYSIS_LOOP) // =====================================

// category cuts ====================================================

std::sort( particles.begin(), particles.end(),
  [](const fj::PseudoJet& a, const fj::PseudoJet& b){
    return ( a.pt() > b.pt() );
  });

const auto H_y = higgs.Rapidity();
bin_t::id<H_rapidity_cut>( std::abs(H_y) < 0.1 );

const auto j1_nsub = jets[0].constituents().size();
bin_t::id<fat_jet>(j1_nsub);

if (particles.size()>=1)
  bin_t::id<p1_flavor>(_kf[particles[0].user_index()]==21 ? 1 : 2 );
if (particles.size()>=2)
  bin_t::id<p2_flavor>(_kf[particles[1].user_index()]==21 ? 1 : 2 );

// ==================================================================

h_H_pT(higgs.Pt());
h_j1_pT(jets[0].pt());

h_H_y(H_y);

const auto ang_Hj1 = ang(higgs,
  { jets[0][0], jets[0][1], jets[0][2], jets[0][3] }
);

h_H_cosTheta(ang_Hj1.cos_theta);
h_Hj_mass(ang_Hj1.M);

h_H_cosTheta_Hj_mass(ang_Hj1.cos_theta,ang_Hj1.M);

h_j1_nsub(j1_nsub);

const double j1_mass = jets[0].m();
h_j1_mass(j1_mass);

const double j1_x = j1_mass/(jet_def.R()*jets[0].pt());
h_j1_x(j1_x);

h_Hj_mass_j1_x(ang_Hj1.M,j1_x);
h_j1_mass_j1_x(j1_mass,j1_x);

h_H_pT_Hj_mass(higgs.Pt(),ang_Hj1.M);
h_j1_pT_Hj_mass(jets[0].pt(),ang_Hj1.M);

h_pp_pTrat_Hj_mass(particles[0].pt()/particles[1].pt(),ang_Hj1.M);
h_pp_dphi_Hj_mass(std::abs(particles[0].delta_phi_to(particles[1])),ang_Hj1.M);

if (jets.size() < 2) continue; // --------------------------------

h_j2_pT(jets[1].pt());

#endif
