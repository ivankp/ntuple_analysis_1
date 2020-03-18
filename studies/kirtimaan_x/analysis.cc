#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

// #include <boost/preprocessor/repetition/repeat.hpp>

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

#define REPEAT(r, data, elem) \
  BOOST_PP_SEQ_ELEM(0,data)(BOOST_PP_CAT(BOOST_PP_CAT( \
    BOOST_PP_SEQ_ELEM(1,data), elem), BOOST_PP_SEQ_ELEM(2,data) ))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(x)(), \
  (1)(2))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(x_)(_3), \
  (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(x_)(_3_zoom), \
  (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(x_)(_4), \
  (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)\
  (10)(11)(12)(13)(14)(15)(16)(17)(18)(19)\
  (20)(21)(22)(23)(24)(25)(26)(27)(28)(29)\
  (30))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(s)(), \
  (12)(34)(45)(35))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(t)(), \
  (23)(15)(34)(35)(45))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(sqrt_s)(), \
  (12)(34)(45)(35))

BOOST_PP_SEQ_FOR_EACH( REPEAT, (h_)(sqrt_t)(), \
  (23)(15)(34)(35)(45))

#define FILL_zoom(X1) h_##X1##_zoom(X1);

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

FILL(s12)
FILL(t23)

h_sqrt_s12(sqrt(s12));
h_sqrt_t23(sqrt(t23));

for (bool swap45 : {false,true}) {
  bin_t::id<swap_45>(int( swap45 ? swap_45::yes : swap_45::no ));

  const auto& k4 = swap45 ? j2 : j1;
  const auto& k5 = swap45 ? j1 : j2;

  def_s(3,4)
  def_s(3,5)
  def_s(4,5)
  def_t(1,5)
  def_t(3,4)
  def_t(3,5)
  def_t(4,5)

  BOOST_PP_SEQ_FOR_EACH( REPEAT, (FILL)(s)(), \
    (34)(35)(45))
  BOOST_PP_SEQ_FOR_EACH( REPEAT, (FILL)(t)(), \
    (15)(34)(35)(45))

  h_sqrt_s34(ssrt(s34));
  h_sqrt_s35(ssrt(s35));
  h_sqrt_s45(ssrt(s45));
  h_sqrt_t15(ssrt(t15));
  h_sqrt_t34(ssrt(t34));
  h_sqrt_t35(ssrt(t35));
  h_sqrt_t45(ssrt(t45));

  // Triangles
  const auto x_0_3 = mh2 - 4*mt2 - s12 + s45 - t23;
  const auto x_1_3 = s34 - 4*mt2;
  const auto x_2_3 = mh2 - 4*mt2 + s12 - s34 - s45;
  const auto x_3_3 = s45 - 4*mt2;
  const auto x_4_3 = s12 - 4*mt2;
  const auto x_5_3 = 4*mt2 + s12 - s34 + t15;
  const auto x_6_3 = 4*mt2 + s45 + t15 - t23;
  const auto x_7_3 = mh2 - 4*mt2 - s34 + t15 - t23;
  const auto x_8_3 =
    mh2*(s12*s45 - 2*(s12+s45)*mt2) + mh2*mh2*mt2 + sq(s12-s45)*mt2;
  const auto x_9_3 =
    mh2*(2*mt2*(s12-s34+s45+2*t15-t23)+(s12-s34+t15)*(s45+t15-t23))
    + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23);
  const auto x_10_3 =
    t15*mh2*mh2 + mt2*sq(s34+t23) - t15*mh2*(4*mt2+s34-t15+t23);

  BOOST_PP_SEQ_FOR_EACH( REPEAT, (FILL)(x_)(_3), \
    (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10))

  BOOST_PP_SEQ_FOR_EACH( REPEAT, (FILL_zoom)(x_)(_3), \
    (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10))

  // Boxes Double cut
  const auto x_0_4 = s12 - 4*mt2;
  const auto x_1_4 = s34 - 4*mt2;
  const auto x_2_4 = mh2 - 4*mt2 + s12 - s34 - s45;
  const auto x_3_4 = s45 - 4*mt2;
  const auto x_4_4 = 4*mt2 + s45 + t15 - t23;
  const auto x_5_4 = 4*mt2 - t23;
  const auto x_6_4 = 4*mt2 + s12 - s34 + t15;
  const auto x_7_4 = mh2 - 4*mt2 - s12 + s45 - t23;
  const auto x_8_4 = mh2 - 4*mt2 - s34 + t15 - t23;
  const auto x_9_4 = -mh2 + 4*mt2 + s12 - s45 + t23;

  // Boxes Triple cut
  const auto x_10_4 = mh2*(s12*s45 - 2*(s12+s45)*mt2)
                    + mh2*mh2*mt2 + sq(s12-s45)*mt2;
  const auto x_11_4 = mh2*( 2*mt2*(s12-s34+s45+2*t15-t23)
                          + (s12-s34+t15)*(s45+t15-t23) )
                    + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23);
  const auto x_12_4 = t15*mh2*mh2 + mt2*sq(s34+t23)
                    - t15*mh2*(4*mt2+s34-t15+t23);

  // Boxes Quadruple cut
  const auto x_13_4 = (s12-s34)*(mh2-s34) + s45*(s34-4*mh2);
  const auto x_14_4 = 4*mt2*(s34-s12)*(mh2-s34) + s34*s45*(s34-4*mt2);
  const auto x_15_4 = -s45*(mh2*mh2 + sq(-s12+s34+s45))
                    + 2*mh2*(2*(s12-s34)*mt2 + s45*(-s12+s34+s45))
                    + 4*s34*(-s12+s34+s45)*mt2;
  const auto x_16_4 = 4*mh2*mt2*(s45+t15)
                    - t23*(4*(-s12+s34+s45)*mt2 + t23*(s12-s34+t15));
  const auto x_17_4 = mh2*(s45+t15) - 4*mt2*(s12-s34+t15) + t23*(s12-s34-s45);
  const auto x_18_4 = mh2*mh2*(s12-s34+t15) // TODO
  const auto x_19_4 =
  const auto x_20_4 =
  const auto x_21_4 =
  const auto x_22_4 =
  const auto x_23_4 =
  const auto x_24_4 =
  const auto x_25_4 =
  const auto x_26_4 =
  const auto x_27_4 =
  const auto x_28_4 =
  const auto x_29_4 =
  const auto x_30_4 =

  BOOST_PP_SEQ_FOR_EACH( REPEAT, (FILL)(x_)(_4), \
    (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)\
    (10)(11)(12)(13)(14)(15)(16)(17)(18)(19)\
    (20)(21)(22)(23)(24)(25)(26)(27)(28)(29)\
    (30))
}

bin_t::id<swap_45>(int( swap_45::all ));

#endif
