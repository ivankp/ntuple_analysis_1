#ifndef ANALYSIS
#ifndef LOOPSIM
#define ANALYSIS "hist_Hjets.hh"
#include "analysis.hh"
#else
#define ANALYSIS_LS "hist_Hjets.hh"
#include "loopsim_analysis.hh"
#endif

#elif defined(ANALYSIS_GLOBAL) // ===================================

#include <boost/preprocessor/repetition/repeat.hpp>
#include <TLorentzVector.h>
#include <fastjet/ClusterSequence.hh>

#include "ivanp/math/math.hh"
#include "ivanp/timed_counter.hh"
#include "ivanp/binner.hh"
#include "ivanp/binner/category_bin.hh"
#include "ivanp/binner/re_axes.hh"
#include "ivanp/container.hh"
#include "ivanp/scope.hh"
#include "bin_defs.hh"

#include <TDirectory.h>
#include <TH1.h>
#include "ivanp/root/binner.hh"
#include <boost/preprocessor/seq/reverse.hpp>

#include "ivanp/sqlite.hh"
#include <cstdio>

#include "json/JetAlgorithm.hh"
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

#include HIST_HJ

#ifndef CATEGORIES
#define CATEGORIES (photon_cuts)(isp)
#endif

template <typename Cat>
constexpr bool used_cat
= ivanp::is_one_of<Cat,BOOST_PP_SEQ_ENUM(CATEGORIES)>::value;

using nlo_bin_t = nlo_bin<double[]>;
using bin_t = ivanp::category_bin<nlo_bin_t,BOOST_PP_SEQ_ENUM(CATEGORIES)>;
template <bool... OF>
using hist = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<typename re_axes::axis_type,OF,OF>...> >;

void excl_labels(TH1* h, bool excl) {
  auto* ax = h->GetXaxis();
  for (int i=1, n=h->GetNbinsX(); i<=n; ++i)
    ax->SetBinLabel(i,cat(excl ? "=" : ">=", i-1).c_str());
}

namespace ivanp { namespace root {
template <> struct bin_converter<bin_t> {
  auto val (const bin_t& b) const noexcept { return (**b).w; }
  auto num (const bin_t& b) const noexcept { return b->n; }
  auto err2(const bin_t& b) const noexcept { return (**b).w2+sq((**b).wtmp); }
};
}}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

TLorentzVector operator+(const TLorentzVector& a, const fastjet::PseudoJet& b) {
  return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] };
}
TLorentzVector& operator+=(TLorentzVector& a, const fastjet::PseudoJet& b) {
  a[0] += b[0]; a[1] += b[1]; a[2] += b[2]; a[3] += b[3]; return a;
}

#elif defined(ANALYSIS_INIT) // =====================================
// this is pasted inside main()

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
const unsigned njets_min = conf.value("/jets/njets_min"_jp,0);

const std::string bfname = conf.at("binning");
cout << "\033[36mBinning\033[0m: " << bfname << '\n' << endl;
re_axes ra(bfname);

// Define histograms ================================================
nlo_bin_t::weights.resize(weights.size());

ivanp::binner<bin_t, std::tuple<ivanp::axis_spec<
    ivanp::uniform_axis<int>,0,1
  >>> h_Njets({njets_min+2u,0,(int)njets_min+2});

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

#define FILL_MACRO(_1,_2,_3,NAME,...) NAME
#define FILL(...) FILL_MACRO(__VA_ARGS__, FILL3_, FILL2_, FILL1_)(__VA_ARGS__)

#define FILL1_(X1) h_##X1(X1);
#define FILL2_(X1,X2) h_##X1##_##X2(X1,X2);
#define FILL3_(X1,X2,X3) h_##X1##_##X2##_##X3(X1,X2,X3);

const double jet_pt_cut  = conf.at("/jets/cuts/pT"_jp);
const double jet_eta_cut = conf.at("/jets/cuts/eta"_jp);
const fj::JetDefinition jet_def = conf.at("/jets/alg"_jp);

Int_t prev_id = -1;
size_t num_entries=0, num_events=0, ncount_total=0;

std::vector<fj::PseudoJet> particles, jets;
Higgs2diphoton Hdecay;
Higgs2diphoton::photons_type photons;
TLorentzVector higgs;

#include HIST_HJ

fj::ClusterSequence::print_banner(); // get it out of the way
cout << jet_def.description() << endl;
cout << "\033[36mNjets\033[0m >= " << njets_min << endl;

#elif defined(ANALYSIS_LOOP) // =====================================

#ifndef LOOPSIM
// Reset ------------------------------------------------------------
const unsigned np = *_nparticle;
particles.clear();

for (unsigned i=weights.size(); i--; ) // set weights
  nlo_bin_t::weights[i] = weights[i];

// Keep track of multi-entry events ---------------------------------
++num_entries;
nlo_bin_t::current_id = *_id;
const bool new_id = (prev_id != nlo_bin_t::current_id);
if (new_id) {
  prev_id = nlo_bin_t::current_id;
  ncount_total += (_ncount ? **_ncount : 1);
  ++num_events;
}

// Read particles ---------------------------------------------------
unsigned n22 = 0; // number of photons
unsigned n25 = 0; // number of Higgs
for (unsigned i=0; i<np; ++i) {
  if (_kf[i] == 25) {
    if (n25) throw std::runtime_error("more than one Higgs");
    higgs.SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
    ++n25;
  } else if (_kf[i] == 22) {
    if (n22>1) throw std::runtime_error("more than two photons");
    photons[n22].SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
    ++n22;
  } else {
    particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
    particles.back().set_user_index(i);
  }
}
if (!n25 && n22!=2) throw std::runtime_error("missing Higgs or photons");
// ------------------------------------------------------------------

// Jets -------------------------------------------------------------
auto fj_seq = fj::ClusterSequence(particles,jet_def);
jets = fj_seq.inclusive_jets(); // get clustered jets

#endif // not loopsim

// apply jet cuts
jets.erase( std::remove_if( jets.begin(), jets.end(), [=](const auto& jet){
  return (jet.pt() < jet_pt_cut)
  or (std::abs(jet.eta()) > jet_eta_cut);
}), jets.end() );
// sort by pT
std::sort( jets.begin(), jets.end(),
  [](const auto& a, const auto& b){ return ( a.pt() > b.pt() ); });
// resulting number of jets
const unsigned njets = jets.size();

// Higgs decay or add photons ---------------------------------------
if (n25) photons = Hdecay(higgs,new_id);
else higgs = photons[0] + photons[1];

double A1_pT = photons[0].Pt(), A2_pT = photons[1].Pt();
if (A1_pT < A2_pT) {
  std::swap(photons[0],photons[1]);
  std::swap(A1_pT,A2_pT);
}

// Photon cuts ------------------------------------------------------
const double A1_eta = photons[0].Eta(), A2_eta = photons[1].Eta();

bin_t::id<photon_cuts>(used_cat<photon_cuts> ? !(
  (A1_pT < 0.35*125.) or
  (A2_pT < 0.25*125.) or
  photon_eta_cut(std::abs(A1_eta)) or
  photon_eta_cut(std::abs(A2_eta))
) : 0);

bin_t::id<isp>(used_cat<isp> ? (unsigned)get_isp(*_id1,*_id2) : 0);

// Fill Histograms --------------------------------------------------
SCOPE_EXIT { h_Njets.fill_bin(njets); };

#include HIST_HJ

#elif defined(ANALYSIS_END) // ======================================

#include HIST_HJ

nlohmann::json info;

info["runcard"] = runcards;

auto& cnt = info["count"];
cnt["entries"] = num_entries;
cnt["events"] = num_events;
cnt["ncount"] = ncount_total;
cnt["norm"] = ncount_total;

cout << "\033[36mPreparing output\033[0m" << endl;

// ==================================================================

#ifndef HIST_MAX_D
#define HIST_MAX_D 1
#endif

auto h_Njets_incl = h_Njets;
h_Njets_incl.integrate_left();

const string& ofname = runcards["output"];
if (ivanp::ends_with(ofname,".root")) {

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(ofname.c_str(),"recreate");
  if (fout->IsZombie()) return 1;
  TDirectory *dir = fout.get();

  // write root historgrams
  nlo_bin_t::wi = 0;
  for (const auto& w : weights_names) {
    dir = dir->mkdir(w.c_str());
    cout << dir->GetName() << endl;

#define CATEGORY_TOP(r, data, elem) \
  bin_t::id<elem>() = 0; \
  for (const char* dir_name : enum_traits<elem>::all_str()) { \
    dir = dir->mkdir(dir_name);
    BOOST_PP_SEQ_FOR_EACH(CATEGORY_TOP,,CATEGORIES)
#undef CATEGORY_TOP

      dir->cd();

      using ivanp::root::to_root;
      using ivanp::root::slice_to_root;

      auto* _h_Njets_excl = to_root(h_Njets,"Njets_excl");
      excl_labels(_h_Njets_excl,true);
      auto* _h_Njets_incl = to_root(h_Njets_incl,"Njets_incl");
      _h_Njets_incl->SetEntries( _h_Njets_excl->GetEntries() );
      excl_labels(_h_Njets_incl,false);

      for (auto& h : hist<1>::all) to_root(*h,h.name);
      for (auto& h : hist<1,0>::all) {
        const auto vars = ivanp::rsplit<1>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1]);
      }
      for (auto& h : hist<1,0,0>::all) {
        const auto vars = ivanp::rsplit<2>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1],vars[2]);
      }

#define CATEGORY_BOT(r, data, elem) \
    dir = dir->GetMotherDir(); \
    ++bin_t::id<elem>(); \
  }
    BOOST_PP_SEQ_FOR_EACH(CATEGORY_BOT,,BOOST_PP_SEQ_REVERSE(CATEGORIES))
#undef CATEGORY_BOT

    dir = dir->GetMotherDir();
    ++nlo_bin_t::wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount_total);
  h_N->SetEntries(num_events);

  fout->Write();

} else if (ivanp::ends_with(ofname,".db")) {

  std::remove(ofname.c_str());
  ivanp::sqlite db(ofname);
  db("PRAGMA page_size=4096")
    ("PRAGMA cache_size=8000");

  db("CREATE TABLE runcard(runcard TEXT)");
  { auto stmt = db.prepare("INSERT INTO runcard VALUES (?)");
    stmt.bind_row(runcards.dump());
  }

  db("CREATE TABLE num(value TEXT, n INT)");
  { auto stmt = db.prepare("INSERT INTO num VALUES (?,?)");
    stmt.bind_row("entries",num_entries)
        .bind_row("events",num_events)
        .bind_row("ncount",ncount_total)
        .bind_row("norm",ncount_total);
  }

  using axis_ptr = const typename hist<1>::axis_type<0>*;
  auto axis_cmp = [](axis_ptr a, axis_ptr b) noexcept { return (*a) < (*b); };
  std::set<axis_ptr,decltype(axis_cmp)> axes_set(axis_cmp);
  std::map<axis_ptr,decltype(axes_set)::iterator> axes_map;
  for (auto& h : hist<1>::all) {
    const axis_ptr a = &h->axis();
    axes_map.emplace(a,axes_set.emplace(a).first);
  }
  std::map<axis_ptr,int> axes_index;
  for (const auto& a : axes_map)
    axes_index.emplace(
      a.first,
      std::distance(axes_set.begin(),a.second) + 2
    );

  db("CREATE TABLE axis_dict(uniform INT, edges TEXT)");
  { auto stmt = db.prepare("INSERT INTO axis_dict VALUES (?,?)");

    auto write = [&](const auto& axis){
      if (axis.is_uniform()) {
        stmt.bind(1,1);
        stmt.bind(2,nlohmann::json{axis.nbins(),axis.min(),axis.max()}.dump());
      } else {
        stmt.bind(1,0);
        nlohmann::json edges;
        auto nbins = axis.nbins();
        for (decltype(nbins) i=0; i<nbins; ++i)
          edges.push_back(axis.edge(i));
        stmt.bind(2,edges.dump());
      }

      stmt.step();
      stmt.reset();
    };

    db("BEGIN");
    write(h_Njets.axis());
    for (auto a : axes_set) write(*a);
    db("END");
  }

#define LAMBDA(r, data, elem) ",\n " BOOST_PP_STRINGIZE(elem TEXT)

  db("CREATE TABLE hist(\n weight TEXT"
    BOOST_PP_SEQ_FOR_EACH(LAMBDA,,CATEGORIES)
    ",\n var1 TEXT,\n _axis INT,\n _head TEXT,\n _data BLOB\n)");

#undef LAMBDA
#define LAMBDA(r, data, elem) ",?"

{ auto stmt = db.prepare("INSERT INTO hist VALUES (?"
    BOOST_PP_SEQ_FOR_EACH(LAMBDA,,CATEGORIES) ",?,?,?,?)");

  static constexpr int n = BOOST_PP_SEQ_SIZE(CATEGORIES) + 2;
  stmt.bind(n+2,"f8,f8,u4",-1,false); // data header

  auto write = [&](const auto& h, const auto& name){
    stmt.bind(n,name,-1,false); // variable name

    const auto& hbins = h.bins();
    const unsigned nbins = hbins.size();

    size_t data_len = nbins*(
      sizeof(double)+sizeof(double)+sizeof(uint32_t)
    );
    void* data = malloc(data_len);
    auto push = [p=data](auto x) mutable {
      (*(reinterpret_cast<decltype(x)*&>(p)++)) = x;
    };
    for (unsigned i=0; i<nbins; ++i) {
      const auto& b = *hbins[i];
      const auto& w = *b;
      push(w.w);
      push(w.w2);
      push((uint32_t)b.n);
    }

    stmt.bind(n+3,data,data_len);
    free(data);

    stmt.step();
    stmt.reset();
  };

  db("BEGIN");

#undef LAMBDA

  nlo_bin_t::wi = 0;
  for (const auto& col : weights_names) {
    stmt.bind(1,col,-1,false);

#define CATEGORY_TOP(r, data, I, elem) \
  bin_t::id<elem>() = 0; \
  for (const char* col : enum_traits<elem>::all_str()) { \
    stmt.bind(I+2,col,-1,false);
    BOOST_PP_SEQ_FOR_EACH_I(CATEGORY_TOP,,CATEGORIES)
#undef CATEGORY_TOP

      stmt.bind(n+1,1); // axis
      write(h_Njets,"Njets_excl");
      write(h_Njets_incl,"Njets_incl");

      for (auto& h : hist<1>::all) {
        stmt.bind(n+1,axes_index[&h->axis()]); // axis
        write(*h,h.name);
      }

#define CATEGORY_BOT(r, data, elem) \
    ++bin_t::id<elem>(); \
  }
    BOOST_PP_SEQ_FOR_EACH(CATEGORY_BOT,,BOOST_PP_SEQ_REVERSE(CATEGORIES))
#undef CATEGORY_BOT

    ++nlo_bin_t::wi;
  }

  db("END");
}

}

#endif // ANALYSIS

