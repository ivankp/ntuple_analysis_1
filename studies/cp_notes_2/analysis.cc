#ifndef ANALYSIS

#define OUTPUT_ROOT

#define HIST_HJ "analysis.cc"
#include "hist_Hjets.hh"

#elif defined(ANALYSIS_GLOBAL) // ===================================

MAKE_ENUM(swap_45,(all)(no)(yes))

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

#elif defined(ANALYSIS_INIT) // =====================================

branch_reader<Double_t> _x1(reader,"x1"), _x2(reader,"x2");

h_(x1) h_(x2)

// const double mh = 125., mh2 = sq(mh), mt = 172.3, mt2 = sq(mt);
// const double mt_a = mt*0.9, mt_b = mt*1.1;

const double rS = runcards["/analysis/rootS"_jp];
const double beam_E = rS*500;

h_(Hj1_mass)
h_(Hj1_mass,j1_pT)
h_(Hj1_mass,j2_pT)
h_(Hj1_mass,j1_pT,j2_pT)

h_(s12)
h_(jj_mass)

h_(cos_Hj1_j2)
h_(cos_Hj2_j1)
h_(cos_Hj1_k1)

#elif defined(ANALYSIS_LOOP) // =====================================

const auto x1 = *_x1;
const auto x2 = *_x2;

FILL(x1)
FILL(x2)

const double E1 = beam_E * x1;
const double E2 = beam_E * x2;

const auto& k1 = TLorentzVector(0,0,E1,E1);
// const auto& k2 = TLorentzVector(0,0,-E2,E2);

const double s12 = sq(E1+E2)-sq(E1-E2);
FILL(s12)

if (njets < 1) continue;

const auto   Hj1 = higgs + jets[0];
const double Hj1_mass = Hj1.M(); // √s34
const double j1_pT = jets[2].pt();

FILL(Hj1_mass)
FILL(Hj1_mass,j1_pT)

const auto boost_cm_Hj1 = -Hj1.BoostVector();
const auto j1_cm_Hj1 = jets[0] >> boost_cm_Hj1;

if (njets == 1) { // ------------------------------------------------

  const auto k1_cm_Hj1 = k1 >> boost_cm_Hj1;
  const auto cos_Hj1_k1 = cos(j1_cm_Hj1,k1_cm_Hj1);
  FILL(cos_Hj1_k1)

} else if (njets == 2) { // -----------------------------------------

  const double j2_pT = jets[1].pt();

  FILL(Hj1_mass,j2_pT)
  FILL(Hj1_mass,j1_pT,j2_pT)

  const auto jj = jets[0]+jets[1];
  const double jj_mass = jj.m();

  FILL(jj_mass)

  const auto Hj2 = higgs + jets[1];
  const auto boost_cm_Hj2 = -Hj2.BoostVector();
  const auto j2_cm_Hj2 = jets[1] >> boost_cm_Hj2;

  const auto j2_cm_Hj1 = jets[1] >> boost_cm_Hj1;
  const auto j1_cm_Hj2 = jets[0] >> boost_cm_Hj2;

  const auto cos_Hj1_j2 = cos(j1_cm_Hj1,j2_cm_Hj1);
  FILL(cos_Hj1_j2)

  const auto cos_Hj2_j1 = cos(j2_cm_Hj2,j1_cm_Hj2);
  FILL(cos_Hj2_j1)

  // if (2*mt_a < Hj1_mass && Hj1_mass < 2*mt_b)
}

#endif
