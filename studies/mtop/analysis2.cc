#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis2.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(fat_jet,(all)(nsubjets_1)(nsubjets_2))
MAKE_ENUM(p1_flavor,(all)(p1_g)(p1_q))
MAKE_ENUM(p2_flavor,(all)(p2_g)(p2_q))

#define CATEGORIES (isp)(photon_cuts)(fat_jet)(p1_flavor)(p2_flavor)

#elif defined(ANALYSIS_INIT) // =====================================

h_(H_pT) h_(j1_pT) h_(j2_pT)

h_(j1_nsub)

h_(Hj_mass)
h_(Hj_mass,H_pT)
h_(Hj_mass,j1_pT)
h_(Hj_mass,j2_pT)

h_(H_pT,j1_pT)
h_(H_pT,j2_pT)
h_(j1_pT,j2_pT)

#elif defined(ANALYSIS_LOOP) // =====================================

// category cuts ====================================================

std::sort( particles.begin(), particles.end(),
  [](const fj::PseudoJet& a, const fj::PseudoJet& b){
    return ( a.pt() > b.pt() );
  });

const auto j1_nsub = jets[0].constituents().size();
bin_t::id<fat_jet>(j1_nsub);

if (particles.size()>=1)
  bin_t::id<p1_flavor>(_kf[particles[0].user_index()]==21 ? 1 : 2 );
if (particles.size()>=2)
  bin_t::id<p2_flavor>(_kf[particles[1].user_index()]==21 ? 1 : 2 );

// ==================================================================

FILL(j1_nsub)

const auto Hj = higgs + jets[0];

const double
  H_pT = higgs.Pt(),
  j1_pT = jets[0].pt(),
  Hj_mass = Hj.M();

FILL(H_pT)
FILL(j1_pT)
FILL(H_pT,j1_pT)

FILL(Hj_mass)

FILL(Hj_mass,H_pT)
FILL(Hj_mass,j1_pT)

if (jets.size() < 2) continue; // --------------------------------

const double
  j2_pT = jets[1].pt();

FILL(j2_pT)
FILL(H_pT,j2_pT)
FILL(j1_pT,j2_pT)
FILL(Hj_mass,j2_pT)

#endif
