#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(T24_sqrtS,(all)(T24_small)(T24_large))
MAKE_ENUM(swap_34,(all)(no)(yes))

#define CATEGORIES (isp)(photon_cuts)(swap_34)(T24_sqrtS)

TLorentzVector operator>>(TLorentzVector v, const TVector3& b) {
  v.Boost(b);
  return v;
}
TLorentzVector operator>>(const fastjet::PseudoJet& j, const TVector3& b) {
  return TLorentzVector(j[0],j[1],j[2],j[3]) >> b;
}
double angle(const TLorentzVector& a, const TLorentzVector& b) {
  return std::abs(a.Angle(b.Vect()));
}

#elif defined(ANALYSIS_INIT) // =====================================

branch_reader<Double_t> _x1(reader,"x1"), _x2(reader,"x2");

h_(Hj1_mass)
h_(Hj1_mass,j2_pT)
h_(sqrtS,T24)
h_(Q,T24)

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < 2) continue; // -----------------------------------------

// A ----------------------------------------------------------------
const double Hj1_mass = (higgs + jets[0]).M();
const double
  j2_pT = jets[1].pt();

// B ----------------------------------------------------------------
const auto cm_boost = -(higgs + jets[0] + jets[1]).BoostVector();

const double E1 = 6500. * *_x1;
const double E2 = 6500. * *_x2;

const auto cm_i1 = TLorentzVector(0,0,E1,E1) >> cm_boost;
const auto cm_i2 = TLorentzVector(0,0,-E2,E2) >> cm_boost;
const auto cm_higgs = higgs >> cm_boost;
const auto cm_jet1 = jets[0] >> cm_boost;
const auto cm_jet2 = jets[1] >> cm_boost;

const auto& pH = cm_higgs;
for (bool swap34 : {false,true}) {
  bin_t::id<swap_34>( swap34 ? 2 : 1 );

  const auto& p4 = swap34 ? cm_jet1 : cm_jet2;
  const auto& p3 = swap34 ? cm_jet2 : cm_jet1;

  const bool i1_closer_to_4 = angle(p4,cm_i1) < angle(p4,cm_i2);

  const auto& p2 = (i1_closer_to_4 ? cm_i1 : cm_i2);
  const auto& p1 = (i1_closer_to_4 ? cm_i2 : cm_i1);

  const double T24   = std::sqrt(std::abs((p2-p4).M2()));
  const double sqrtS = std::sqrt(std::abs((pH+p3).M2()));
  const double Q     = std::sqrt(std::abs((p1-pH).M2()));

  bin_t::id<T24_sqrtS>( T24 < sqrtS ? 1 : 2 );

  FILL(Hj1_mass)
  FILL(Hj1_mass,j2_pT)

  FILL(sqrtS,T24)
  FILL(Q,T24)
}

#endif
