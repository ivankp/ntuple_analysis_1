#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

#define CATEGORIES (isp)

TLorentzVector operator>>(TLorentzVector v, const TVector3& b) {
  v.Boost(b);
  return v;
}
TLorentzVector operator>>(const fastjet::PseudoJet& j, const TVector3& b) {
  return TLorentzVector(j[0],j[1],j[2],j[3]) >> b;
}
double cos(const TVector3& a, const TVector3& b) noexcept {
  const double ptot2 = a.Mag2()*b.Mag2();
  if (ptot2 <= 0) {
    return 0;
  } else {
    const double c = a.Dot(b)/std::sqrt(ptot2);
    if (c >  1) return  1;
    if (c < -1) return -1;
    return c;
  }
}
double cos(const TLorentzVector& a, const TLorentzVector& b) noexcept {
  return cos(a.Vect(),b.Vect());
}

bool approx(double a, double b, double t) noexcept {
  a /= b;
  return (1-t) < a && a < (1+t);
}

#elif defined(ANALYSIS_INIT) // =====================================

branch_reader<Double_t> _x1(reader,"x1"), _x2(reader,"x2");

h_(x1) h_(x2)

const double /*mh = 125.,*/ mt = 172.3;

const double rS = runcards["/analysis/rootS"_jp];
const double beam_E = rS*500;

h_(cos_Hj1_j2)
  h_(cos_Hj1_j2_cuts1)
  h_(cos_Hj1_j2_cuts1y)
  h_(cos_Hj1_j2_cuts2)
  h_(cos_Hj1_j2_cuts4)
h_(cos_Hj2_j1)
h_(cos_Hj1_k1)

#elif defined(ANALYSIS_LOOP) // =====================================

const auto x1 = *_x1;
const auto x2 = *_x2;

FILL(x1)
FILL(x2)

const double E1 = beam_E * x1;
// const double E2 = beam_E * x2;

const auto& k1 = TLorentzVector(0,0,E1,E1);
// const auto& k2 = TLorentzVector(0,0,-E2,E2);

if (njets < 1) continue;

const auto   Hj1 = higgs + jets[0];
const double Hj1_mass = Hj1.M(); // √s34
const double j1_pT = jets[0].pt();

if (njets == 1 && jets[0].constituents().size() == 2) { // ----------

  const auto boost_cm = -Hj1.BoostVector();
  const auto j1_cm = jets[0] >> boost_cm;

  const auto k1_cm = k1 >> boost_cm;
  const auto cos_Hj1_k1 = cos(j1_cm,k1_cm);
  FILL(cos_Hj1_k1)

} else if (njets == 2) { // -----------------------------------------

  const double H_pT = higgs.Pt();
  const double j2_pT = jets[1].pt();

  // const auto jj = jets[0]+jets[1];
  // const double jj_mass = jj.m();

  const auto boost_cm = -(Hj1+jets[1]).BoostVector();

  const auto H_cm = higgs >> boost_cm;
  const auto j1_cm = jets[0] >> boost_cm;
  const auto j2_cm = jets[1] >> boost_cm;

  const auto Hj1_cm = H_cm + j1_cm;
  const auto boost_cm_Hj1 = -Hj1_cm.BoostVector();

  const auto j1_cm_Hj1 = j1_cm >> boost_cm_Hj1;
  const auto j2_cm_Hj1 = j2_cm >> boost_cm_Hj1;

  const auto Hj2_cm = H_cm + j2_cm;
  const auto boost_cm_Hj2 = -Hj2_cm.BoostVector();

  const auto j2_cm_Hj2 = j2_cm >> boost_cm_Hj2;
  const auto j1_cm_Hj2 = j1_cm >> boost_cm_Hj2;

  const auto cos_Hj1_j2 = cos(j1_cm_Hj1,j2_cm_Hj1);
  FILL(cos_Hj1_j2)

  const auto cos_Hj2_j1 = cos(j2_cm_Hj2,j1_cm_Hj2);
  FILL(cos_Hj2_j1)

  const double  H_y = higgs.Rapidity();
  const double j1_y = jets[0].rap();
  const double j2_y = jets[1].rap();
  const double  H_absy = std::abs( H_y);
  const double j1_absy = std::abs(j1_y);
  const double j2_absy = std::abs(j2_y);

  if ((2*mt < Hj1_mass && Hj1_mass < (2*mt+100)) &&
      // approx(2*mt,Hj1_mass,0.2) &&
      approx(H_pT,j1_pT,0.2) &&
      j1_pT < 2*j2_pT
  ) {
    h_cos_Hj1_j2_cuts1(cos_Hj1_j2);

    if (approx(j1_absy,H_absy,0.2) && H_absy < j2_absy)
      h_cos_Hj1_j2_cuts1y(cos_Hj1_j2);
  }

  if ((2*mt < Hj1_mass && Hj1_mass < (2*mt+100)) &&
      // Hj1_mass > 2*mt &&
      j1_pT > H_pT &&
      approx(j1_absy,j2_absy,0.2) &&
      std::abs(j1_y-H_y) < std::abs(j2_y-H_y) &&
      j2_pT > 0.5*j1_pT
  ) {
    h_cos_Hj1_j2_cuts2(cos_Hj1_j2);
  }

  if ((2*mt < Hj1_mass && Hj1_mass < (2*mt+100)) &&
      approx(j2_pT,H_pT,0.2) &&
      approx(j1_absy,j2_absy,0.2)
  ) {
    h_cos_Hj1_j2_cuts4(cos_Hj1_j2);
  }
}

#endif
