#ifndef ANALYSIS
#define ANALYSIS hist_Hjets.hh
#include "analysis.hh"
#endif

#ifdef ANALYSIS_GLOBAL // ===========================================

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "ivanp/math/math.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/binner.hh"
#include "ivanp/binner/category_bin.hh"
#include "ivanp/binner/re_axes.hh"
#include "bin_defs.hh"
#include "ivanp/container.hh"

#include "json/JetAlgorithm.hh"
#include "json/binner.hh"
#include "json/bins.hh"

#include "Higgs2diphoton.hh"

using ivanp::reserve;

namespace fj = fastjet;
using namespace ivanp::math;

MAKE_ENUM(isp,(all)(gg)(gq)(qq))

isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

MAKE_ENUM(photon_cuts,(all)(with_photon_cuts))

#define HIST_HJ_GLOBAL
#include STR(HIST_HJ)
#undef HIST_HJ_GLOBAL

#ifndef CATEGORIES
#define CATEGORIES (photon_cuts)(isp)
#endif

using cat_bin = ivanp::category_bin<nlo_bin,BOOST_PP_SEQ_ENUM(CATEGORIES)>;

using bin_t = multiweight_bin<cat_bin>;
template <bool... OF>
using hist = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<typename re_axes::axis_type,OF,OF>...> >;

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

TLorentzVector operator+(const TLorentzVector& a, const fastjet::PseudoJet& b) {
  return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] };
}
TLorentzVector& operator+=(TLorentzVector& a, const fastjet::PseudoJet& b) {
  a[0] += b[0]; a[1] += b[1]; a[2] += b[2]; a[3] += b[3]; return a;
}

#endif
#ifdef ANALYSIS_INIT // =============================================

branch_reader<Int_t> _id(reader,"id");
branch_reader<Int_t> _nparticle(reader,"nparticle");
branch_reader<Int_t[]> _kf(reader,"kf");

floats_reader
  _px(reader,"px"), _py(reader,"py"), _pz(reader,"pz"), _E(reader,"E");

boost::optional<branch_reader<Int_t>> _ncount;
for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
  if (!strcmp(bo->GetName(),"ncount")) {
    _ncount.emplace(reader,"ncount"); break;
  }
}
branch_reader<Int_t> _id1(reader,"id1"), _id2(reader,"id2");

const auto& conf = runcards.at("analysis");
const unsigned njets_required = conf.value("/jets/njets"_jp,0);

const std::string bfname = conf.at("binning");
cout << "\033[36mBinning\033[0m: " << bfname << '\n' << endl;
re_axes ra(bfname);

// Define histograms ==============================================
bin_t::weights.resize(_weights.size());

ivanp::binner<bin_t, std::tuple<ivanp::axis_spec<
    ivanp::uniform_axis<int>,0,1
  >>> h_Njets({njets_required+2u,0,(int)njets_required+2});

#define H_MACRO(_1,_2,_3,NAME,...) NAME
#define h_(...) H_MACRO(__VA_ARGS__, h3_, h2_, h1_)(__VA_ARGS__)

#define h1_(X) hist<1> h_##X(#X,ra[#X]);
#define h2_(X1,X2) hist<1,0> \
  h_##X1##_##X2(#X1"-"#X2, \
    ra[#X1":"#X1"-"#X2], \
    ra[#X2":"#X1"-"#X2]);
#define h3_(X1,X2,X3) hist<1,0,0> \
  h_##X1##_##X2##_##X3(#X1"-"#X2"-"#X3, \
    ra[#X1":"#X1"-"#X2"-"#X3], \
    ra[#X2":"#X1"-"#X2"-"#X3], \
    ra[#X3":"#X1"-"#X2"-"#X3]);

const double jet_pt_cut  = conf.at("/jets/cuts/pT"_jp);
const double jet_eta_cut = conf.at("/jets/cuts/eta"_jp);
const fj::JetDefinition jet_def = conf.at("/jets/alg"_jp);

Int_t prev_id = -1;
size_t num_entries=0, num_events=0, ncount_total=0;

std::vector<fj::PseudoJet> partons, jets;
Higgs2diphoton Hdecay;
Higgs2diphoton::photons_type photons;
TLorentzVector higgs;

#define HIST_HJ_INIT
#include STR(HIST_HJ)
#undef HIST_HJ_INIT

fj::ClusterSequence::print_banner(); // get it out of the way
cout << jet_def.description() << endl;
cout << "\033[36mNjets\033[0m >= " << njets_required << endl;

#endif
#ifdef ANALYSIS_LOOP // =============================================

// Reset ----------------------------------------------------------
const size_t np = *_nparticle;
partons.clear();

for (unsigned i=_weights.size(); i--; ) // set weights
  bin_t::weights[i] = *_weights[i];

// Keep track of multi-entry events -------------------------------
++num_entries;
nlo_bin::current_id = *_id;
const bool new_id = (prev_id != nlo_bin::current_id);
if (new_id) {
  prev_id = nlo_bin::current_id;
  ncount_total += (_ncount ? **_ncount : 1);
  ++num_events;
}

// Read particles -------------------------------------------------
unsigned n22 = 0; // number of photons
unsigned n25 = 0; // number of Higgs
for (size_t i=0; i<np; ++i) {
  if (_kf[i] == 25) {
    if (n25) throw std::runtime_error("more than one Higgs");
    higgs.SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
    ++n25;
  } else if (_kf[i] == 22) {
    if (n22>1) throw std::runtime_error("more than two photons");
    photons[n22].SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
    ++n22;
  } else {
    partons.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
  }
}
if (!n25 && n22!=2) throw std::runtime_error("missing Higgs or photons");

cat_bin::id<isp>() = (unsigned)get_isp(*_id1,*_id2);
// ----------------------------------------------------------------

// Jets -----------------------------------------------------------
auto fj_seq = fj::ClusterSequence(partons,jet_def);
jets = fj_seq.inclusive_jets(jet_pt_cut); // apply pT cut
// apply eta cut
for (auto it=jets.begin(); it!=jets.end(); ) {
  if (std::abs(it->eta()) > jet_eta_cut) it = jets.erase(it);
  else ++it;
}
// sort by pT
std::sort( jets.begin(), jets.end(),
  [](const fj::PseudoJet& a, const fj::PseudoJet& b){
    return ( a.pt() > b.pt() );
  });
// resulting number of jets
const unsigned njets = jets.size();

// Higgs decay or add photons -------------------------------------
if (n25) photons = Hdecay(higgs,new_id);
else higgs = photons[0] + photons[1];

double A1_pT = photons[0].Pt(), A2_pT = photons[1].Pt();
if (A1_pT < A2_pT) {
  std::swap(photons[0],photons[1]);
  std::swap(A1_pT,A2_pT);
}

// Photon cuts ----------------------------------------------------
const double A1_eta = photons[0].Eta(), A2_eta = photons[1].Eta();

cat_bin::id<photon_cuts>() = !(
  (A1_pT < 0.35*125.) or
  (A2_pT < 0.25*125.) or
  photon_eta_cut(std::abs(A1_eta)) or
  photon_eta_cut(std::abs(A2_eta))
);

// Fill Histograms ------------------------------------------------
h_Njets.fill_bin(njets);

#define HIST_HJ_LOOP
#include STR(HIST_HJ)
#undef HIST_HJ_LOOP

#endif
#ifdef ANALYSIS_END // ==============================================

nlohmann::json out;
auto& ann = out["annotation"];

ann["runcard"] = runcards;

auto& cnt = ann["count"];
cnt["entries"] = num_entries;
cnt["events"] = num_events;
cnt["ncount"] = ncount_total;
cnt["norm"] = ncount_total;

#define CATEGORY_ANN(r, data, elem) \
  ann_bins.push_back({STR(elem),enum_traits<elem>::all_str()});

auto& ann_bins = ann["bins"];
ann_bins.push_back({"weight",weights_names});
BOOST_PP_SEQ_FOR_EACH(CATEGORY_ANN,,CATEGORIES)
ann_bins.push_back({"bin",{"w","w2","n"}});

auto& hists = out["histograms"];

hists["Njets_excl"] = h_Njets;
auto h_Njets_incl = h_Njets;
h_Njets_incl.integrate_left();
hists["Njets_incl"] = h_Njets_incl;

#ifndef HIST_MAX_D
#define HIST_MAX_D 1
#endif

#define WITH_COMMAS(z, n, text) ,text
#define SAVE_HISTS(z, n, text) \
  for (const auto& h : hist<1 BOOST_PP_REPEAT(n,WITH_COMMAS,0)>::all) \
    hists[h.name] = *h;
BOOST_PP_REPEAT(HIST_MAX_D,SAVE_HISTS,)

const string ofname = runcards["output"];
cout << "\033[36mWriting output\033[0m: " << ofname << endl;

try { // write output file
  std::ofstream file(ofname);
  file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
  if (ivanp::ends_with(ofname,".xz")) {
    namespace bio = boost::iostreams;
    bio::filtering_streambuf<bio::output> buf;
    buf.push(bio::lzma_compressor(bio::lzma::best_compression));
    buf.push(file);
    std::ostream(&buf) << out;
  } else file << out;
} catch (const std::exception& e) {
  cerr << "\033[31mError writing file\033[0m \"" << ofname << "\": "
       << e.what() << endl;
  return 1;
}

#endif

