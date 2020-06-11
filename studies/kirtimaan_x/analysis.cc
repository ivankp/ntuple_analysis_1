#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#include <boost/preprocessor/seq/for_each_product.hpp>

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

#define REPEAT_MACRO(r, x) \
  BOOST_PP_SEQ_ELEM(0,x) ( BOOST_PP_CAT(BOOST_PP_CAT( \
    BOOST_PP_SEQ_ELEM(1,x), BOOST_PP_SEQ_ELEM(3,x)), BOOST_PP_SEQ_ELEM(2,x) ))

#define REPEAT(SEQ) BOOST_PP_SEQ_FOR_EACH_PRODUCT( REPEAT_MACRO, SEQ )

REPEAT( ((h_))((x))(()) ((1)(2)) )

REPEAT( ((h_))((s)(sqrt_s))(()(_zoom1)(_zoom2)) ( (12)(34)(45)(35) ))
REPEAT( ((h_))((t)(sqrt_t))(()(_zoom1)(_zoom2)) ( (23)(15)(34)(35)(45) ))

REPEAT( ((h_))((x_3_))(()(_zoom1)(_zoom2)) ( \
  (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10) ))

REPEAT( ((h_))((x_4_))(()(_zoom1)(_zoom2)) ( \
  (0)(1)(2)(3)(4)(5)(6)(7)(8)(9) \
  (10)(11)(12)(13)(14)(15)(16)(17)(18)(19) \
  (20)(21)(22)(23)(24)(25)(26)(27)(28)(29) \
  (30) ))

const double mh2 = sq(125.), mt2 = sq(172.3);

const double rS = runcards["/analysis/rootS"_jp];
const double beam_E = rS*500;

#elif defined(ANALYSIS_LOOP) // =====================================

if (njets < 2) continue; // -----------------------------------------

const auto x1 = *_x1;
const auto x2 = *_x2;

FILL(x1);
FILL(x2);

const double E1 = beam_E * x1;
const double E2 = beam_E * x2;

const auto& k1 = TLorentzVector(0,0,E1,E1);
const auto& k2 = TLorentzVector(0,0,-E2,E2);
const auto& k3 = higgs;
const auto& j1 = jets[0];
const auto& j2 = jets[1];

#define def_s(i,j) \
  const auto s##i##j = norm2(k##i + k##j); \
  const auto sqrt_s##i##j = ssrt(s##i##j);
#define def_t(i,j) \
  const auto t##i##j = norm2(k##i - k##j); \
  const auto sqrt_t##i##j = ssrt(t##i##j);

def_s(1,2)
def_t(2,3)

#define FILL_zoom1_(X1) h_##X1##_zoom1(X1);
#define FILL_zoom1(X1) FILL_zoom1_(X1)

#define FILL_zoom2_(X1) h_##X1##_zoom2(X1);
#define FILL_zoom2(X1) FILL_zoom2_(X1)

REPEAT( ( (FILL)(FILL_zoom1)(FILL_zoom2) )((s)(sqrt_s))(())( (12) ))
REPEAT( ( (FILL)(FILL_zoom1)(FILL_zoom2) )((t)(sqrt_t))(())( (23) ))

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

  REPEAT( ( (FILL)(FILL_zoom1)(FILL_zoom2) )((s)(sqrt_s))(())( \
        (34)(35)(45) ))
  REPEAT( ( (FILL)(FILL_zoom1)(FILL_zoom2) )((t)(sqrt_t))(())( \
    (15)(34)(35)(45) ))

  // Triangles
  const auto x_3_0 = mh2 - 4*mt2 - s12 + s45 - t23;
  const auto x_3_1 = s34 - 4*mt2;
  const auto x_3_2 = mh2 - 4*mt2 + s12 - s34 - s45;
  const auto x_3_3 = s45 - 4*mt2;
  const auto x_3_4 = s12 - 4*mt2;
  const auto x_3_5 = 4*mt2 + s12 - s34 + t15;
  const auto x_3_6 = 4*mt2 + s45 + t15 - t23;
  const auto x_3_7 = mh2 - 4*mt2 - s34 + t15 - t23;
  const auto x_3_8 =
    mh2*(s12*s45 - 2*(s12+s45)*mt2) + mh2*mh2*mt2 + sq(s12-s45)*mt2;
  const auto x_3_9 =
    mh2*(2*mt2*(s12-s34+s45+2*t15-t23)+(s12-s34+t15)*(s45+t15-t23))
    + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23);
  const auto x_3_10 =
    t15*mh2*mh2 + mt2*sq(s34+t23) - t15*mh2*(4*mt2+s34-t15+t23);

  REPEAT( ((FILL)(FILL_zoom1)(FILL_zoom2))((x_3_))(()) ( \
    (0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10) ))

  // Boxes Double cut
  const auto x_4_0 = s12 - 4*mt2;
  const auto x_4_1 = s34 - 4*mt2;
  const auto x_4_2 = mh2 - 4*mt2 + s12 - s34 - s45;
  const auto x_4_3 = s45 - 4*mt2;
  const auto x_4_4 = 4*mt2 + s45 + t15 - t23;
  const auto x_4_5 = 4*mt2 - t23;
  const auto x_4_6 = 4*mt2 + s12 - s34 + t15;
  const auto x_4_7 = mh2 - 4*mt2 - s12 + s45 - t23;
  const auto x_4_8 = mh2 - 4*mt2 - s34 + t15 - t23;
  const auto x_4_9 = -mh2 + 4*mt2 + s12 - s45 + t23;

  // Boxes Triple cut
  const auto x_4_10 = mh2*(s12*s45 - 2*(s12+s45)*mt2)
                    + mh2*mh2*mt2 + sq(s12-s45)*mt2;
  const auto x_4_11 = mh2*( 2*mt2*(s12-s34+s45+2*t15-t23)
                          + (s12-s34+t15)*(s45+t15-t23) )
                    + mh2*mh2*mt2 + mt2*sq(s12-s34-s45+t23);
  const auto x_4_12 = t15*mh2*mh2 + mt2*sq(s34+t23)
                    - t15*mh2*(4*mt2+s34-t15+t23);

  // Boxes Quadruple cut
  const auto x_4_13 = (s12-s34)*(mh2-s34) + s45*(s34-4*mh2);
  const auto x_4_14 = 4*mt2*(s34-s12)*(mh2-s34) + s34*s45*(s34-4*mt2);
  const auto x_4_15 = -s45*(mh2*mh2 + sq(-s12+s34+s45))
                    + 2*mh2*(2*(s12-s34)*mt2 + s45*(-s12+s34+s45))
                    + 4*s34*(-s12+s34+s45)*mt2;
  const auto x_4_16 = 4*mh2*mt2*(s45+t15)
                    - t23*(4*(-s12+s34+s45)*mt2 + t23*(s12-s34+t15));
  const auto x_4_17 = mh2*(s45+t15) - 4*mt2*(s12-s34+t15) + t23*(s12-s34-s45);
  const auto x_4_18 = mh2*mh2*(s12-s34+t15)
                    + 2*mh2*((s12-s34-s45)*(s12-s34+t15)-2*mt2*(s45+t15))
                    + (s12-s34-s45)*((s12-s34-s45)*(s12-s34+t15)-4*mt2*t23);
  const auto x_4_19 = 4*mt2*(mh2*t15+(s12-s34-s45)*(s12-s45+t23))
                    - t15*sq(mh2-s12+s45-t23);
  const auto x_4_20 = t15*(mh2-4*mt2) + (s12-s34-s45)*(s12-s45+t23);
  const auto x_4_21 = 4*mt2*(mh2*t15 + (s12-s34-s45)*(s12-s45+t23))
                    - t15*sq(mh2+s12-s34-s45);
  const auto x_4_22 = mh2*(4*mt2*t15-t23*t23)
                    + t23*(t23*(s34-t15+t23)-4*mt2*s34);
  const auto x_4_23 = mh2*(t15-4*mt2) + 4*mt2*(s34-t15+t23) - s34*t23;
  const auto x_4_24 = -sq(s34)*(mh2+t15-t23) + 4*mt2*(mh2*t15-s34*t23)
                    + s34*s34*s34;
  const auto x_4_25 = sq(mh2)*(s45+t15-t23)
                    - 2*mh2*(2*mt2*(s12+t15) + (s45+t15-t23)*(s12-s45+t23))
                    - (s12-s45+t23)*(-4*mt2*s34-(s45+t15-t23)*(s12-s45+t23));
  const auto x_4_26 = mh2*(s12+t15) - 4*mt2*(s45+t15-t23) - s34*(s12-s45+t23);
  const auto x_4_27 = 4*mh2*mt2*(s12+t15)
                    + s34*(-4*mt2*(s12-s45+t23)-s34*(s45+t15-t23));
  const auto x_4_28 = mh2*mh2*s12 + (s12-s45+t23)*(s12*(s12-s45+t23)-4*mt2*t23)
                    - 2*mh2*(2*mt2*(s45-t23)+s12*(s12-s45+t23));
  const auto x_4_29 = 4*mh2*mt2*(t23-s45) + t23*(s12*t23-4*mt2*(s12-s45+t23));
  const auto x_4_30 = mh2*(s45-t23) - 4*mt2*s12 + t23*(s12-s45+t23);

  REPEAT( ((FILL)(FILL_zoom1)(FILL_zoom2))((x_4_))(()) ( \
    (0)(1)(2)(3)(4)(5)(6)(7)(8)(9) \
    (10)(11)(12)(13)(14)(15)(16)(17)(18)(19) \
    (20)(21)(22)(23)(24)(25)(26)(27)(28)(29) \
    (30) ))
}

bin_t::id<swap_45>(int( swap_45::all ));

#endif
