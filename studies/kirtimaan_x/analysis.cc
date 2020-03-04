#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(swap_45,(all)(no)(yes))

#define CATEGORIES (isp)(photon_cuts)(swap_45)

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

h_(x_0_3)
h_(x_1_3)
h_(x_2_3)
h_(x_3_3)
h_(x_4_3)
h_(x_5_3)
h_(x_6_3)
h_(x_7_3)
h_(x_8_3)
h_(x_9_3)
h_(x_10_3)

const double mh2 = sq(125.), mt2 = sq(172.3);

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < 2) continue; // -----------------------------------------

const auto cm_boost = -(higgs + jets[0] + jets[1]).BoostVector();

const double E1 = 6500. * *_x1;
const double E2 = 6500. * *_x2;

const auto k1 = TLorentzVector(0,0,E1,E1) >> cm_boost;
const auto k2 = TLorentzVector(0,0,-E2,E2) >> cm_boost;
const auto k3 = higgs >> cm_boost;
const auto j1 = jets[0] >> cm_boost;
const auto j2 = jets[1] >> cm_boost;

#define def_s(i,j) const auto s##i##j = (k##i + k##j).M2();
#define def_t(i,j) const auto t##i##j = (k##i - k##j).M2();

def_s(1,2)
def_t(2,3)

for (bool swap45 : {false,true}) {
  bin_t::id<swap_45>( swap45 ? 2 : 1 );

  const auto& k4 = swap45 ? j2 : j1;
  const auto& k5 = swap45 ? j1 : j2;

  def_s(3,4)
  def_s(4,5)
  def_t(1,5)

  h_x_0_3(mh2 - 4*mt2 - s12 + s45 - t23);
  h_x_1_3(s34 - 4*mt2);
  h_x_2_3(mh2 - 4*mt2 + s12 - s34 - s45);
  h_x_3_3(s45 - 4*mt2);
  h_x_4_3(s12 - 4*mt2);
  h_x_5_3(4*mt2 + s12 - s34 + t15);
  h_x_6_3(4*mt2 + s45 + t15 - t23);
  h_x_7_3(mh2 - 4*mt2 - s34 + t15 - t23);
  h_x_8_3(mh2*(s12*s45 - 2*(s12+s45)*mt2) + mh2*mh2*mt2 + sq(s12-s45)*mt2);
  h_x_9_3(mh2*(2*mt2*(s12-s34+s45+2*t15-t23)+(s12-s34+t15)*(s45+t15-t23))
          + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23));
  h_x_10_3(t15*mh2*mh2 + mt2*sq(s34+t23) - t15*mh2*(4*mt2+s34-t15+t23));
}

#endif
