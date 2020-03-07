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

template <typename T>
double norm2(const T& k) noexcept { return sq(k[3]) - sq(k[0],k[1],k[2]); }

double ssrt(double x) noexcept {
  return x < 0 ? -std::sqrt(-x) : std::sqrt(x);
}

#elif defined(ANALYSIS_INIT) // =====================================

branch_reader<Double_t> _x1(reader,"x1"), _x2(reader,"x2");

h_(x_0_3_zoom)
h_(x_1_3_zoom)
h_(x_2_3_zoom)
h_(x_3_3_zoom)
h_(x_4_3_zoom)
h_(x_5_3_zoom)
h_(x_6_3_zoom)
h_(x_7_3_zoom)
h_(x_8_3_zoom)
h_(x_9_3_zoom)
h_(x_10_3_zoom)

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

h_(x1)
h_(x2)

h_(s12)
h_(s34)
h_(s45)
h_(t23)
h_(t15)

const double mh2 = sq(125.), mt2 = sq(172.3);

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < 2) continue; // -----------------------------------------

const auto x1 = *_x1;
const auto x2 = *_x2;

FILL(x1);
FILL(x2);

const double E1 = 6500. * x1;
const double E2 = 6500. * x2;

const auto& k1 = TLorentzVector(0,0,E1,E1);
const auto& k2 = TLorentzVector(0,0,-E2,E2);
const auto& k3 = higgs;
const auto& j1 = jets[0];
const auto& j2 = jets[1];

#define def_s(i,j) const auto s##i##j = norm2(k##i + k##j);
#define def_t(i,j) const auto t##i##j = norm2(k##i - k##j);

def_s(1,2)
def_t(2,3)

h_s12(ssrt(s12));
h_t23(ssrt(t23));

for (bool swap45 : {false,true}) {
  bin_t::id<swap_45>(int( swap45 ? swap_45::yes : swap_45::no ));

  const auto& k4 = swap45 ? j2 : j1;
  const auto& k5 = swap45 ? j1 : j2;

  def_s(3,4)
  def_s(4,5)
  def_t(1,5)

  h_s34(ssrt(s34));
  h_s45(ssrt(s45));
  h_t15(ssrt(t15));

  const auto x_0_3 = ssrt(mh2 - 4*mt2 - s12 + s45 - t23);
  const auto x_1_3 = ssrt(s34 - 4*mt2);
  const auto x_2_3 = ssrt(mh2 - 4*mt2 + s12 - s34 - s45);
  const auto x_3_3 = ssrt(s45 - 4*mt2);
  const auto x_4_3 = ssrt(s12 - 4*mt2);
  const auto x_5_3 = ssrt(4*mt2 + s12 - s34 + t15);
  const auto x_6_3 = ssrt(4*mt2 + s45 + t15 - t23);
  const auto x_7_3 = ssrt(mh2 - 4*mt2 - s34 + t15 - t23);
  const auto x_8_3 = ssrt((
      mh2*(s12*s45 - 2*(s12+s45)*mt2) + mh2*mh2*mt2 + sq(s12-s45)*mt2
    ) / (mh2*mt2) );
  const auto x_9_3 = ssrt((
      mh2*(2*mt2*(s12-s34+s45+2*t15-t23)+(s12-s34+t15)*(s45+t15-t23))
      + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23)
    ) / (mh2*mt2) );
  const auto x_10_3 = ssrt((
      t15*mh2*mh2 + mt2*sq(s34+t23) - t15*mh2*(4*mt2+s34-t15+t23)
    ) / (mh2*mt2) );

  FILL(x_0_3)
  FILL(x_1_3)
  FILL(x_2_3)
  FILL(x_3_3)
  FILL(x_4_3)
  FILL(x_5_3)
  FILL(x_6_3)
  FILL(x_7_3)
  FILL(x_8_3)
  FILL(x_9_3)
  FILL(x_10_3)

  h_x_0_3_zoom(x_0_3);
  h_x_1_3_zoom(x_1_3);
  h_x_2_3_zoom(x_2_3);
  h_x_3_3_zoom(x_3_3);
  h_x_4_3_zoom(x_4_3);
  h_x_5_3_zoom(x_5_3);
  h_x_6_3_zoom(x_6_3);
  h_x_7_3_zoom(x_7_3);
  h_x_8_3_zoom(x_8_3);
  h_x_9_3_zoom(x_9_3);
  h_x_10_3_zoom(x_10_3);
}

bin_t::id<swap_45>(int( swap_45::all ));

#endif
