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

MAKE_ENUM(H_y_cat,(all)(central))

#define CATEGORIES (H_y_cat)(photon_cuts)(isp)

#elif defined(ANALYSIS_INIT) // =====================================

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

#elif defined(ANALYSIS_LOOP) // =====================================

const auto H_y = higgs.Rapidity();
bin_t::id<H_y_cat>( std::abs(H_y) < 0.1 );

h_H_y(H_y);

const auto ang_Hj1 = ang(higgs,
  { jets[0][0], jets[0][1], jets[0][2], jets[0][3] }
);

h_H_cosTheta(ang_Hj1.cos_theta);
h_Hj_mass(ang_Hj1.M);

h_H_cosTheta_Hj_mass(ang_Hj1.cos_theta,ang_Hj1.M);

const auto j1_nsub = jets[0].constituents().size();
h_j1_nsub(j1_nsub);

if (j1_nsub>1) {
  const double j1_mass = jets[0].m();
  h_j1_mass(j1_mass);

  const double j1_x = j1_mass/(jet_def.R()*jets[0].pt());
  h_j1_x(j1_x);

  h_Hj_mass_j1_x(ang_Hj1.M,j1_x);
  h_j1_mass_j1_x(j1_mass,j1_x);
}

if (jets.size() < 2) continue;

h_j2_mass(jets[1].m());

#endif
