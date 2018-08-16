#include "analysis.hh"
#include <TLorentzVector.h>
#include <fastjet/ClusterSequence.hh>
#include "ivanp/math.hh"
#include "ivanp/enum_traits.hh"
#include "ivanp/binner.hh"
#include "ivanp/binner/re_axis.hh"
#include "ivanp/binner/category_bin.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
// using ivanp::cat;
namespace fj = fastjet;

#define CATEGORIES (photon_cuts)(isp)

MAKE_ENUM(isp,(all)(gg)(gq)(qq))

isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

MAKE_ENUM(photon_cuts,(all)(with_photon_cuts))

using cat_bin = ivanp::category_bin<nlo_bin,BOOST_PP_SEQ_ENUM(CATEGORIES)>;

using bin_t = multiweight_bin<cat_bin>;
template <bool... OF>
using hist = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<typename re_axes::axis_type,OF,OF>...> >;

#define H_MACRO(_1,_2,_3,NAME,...) NAME
#define h_(...) H_MACRO(__VA_ARGS__, h3_, h2_, h1_)(__VA_ARGS__)

#define h1_(X) hist<1> h_##X(#X,ra[#X]);
#define h2_(X1,X2) hist<1,0> h_##X1##_##X2(#X1"-"#X2,ra[#X1"_2"],ra[#X2"_2"]);
#define h3_(X1,X2,X3) hist<1,0,0> \
  h_##X1##_##X2##_##X3(#X1"-"#X2"-"#X3,ra[#X1"_2"],ra[#X2"_2"],ra[#X3"_2"]);

#define hj_(X) auto h_jet_##X = reserve<hist<1>>(njets_expected+1); \
  for (unsigned i=0; i<=njets_expected; ++i) { \
    const auto name = cat("jet",i+1,"_"#X); \
    h_jet_##X.emplace_back(name,ra[name]); \
  }

class diphoton_analysis_base: analysis_base {
  TTreeReaderValue<Int_t> _id1(reader,"id1"), _id2(reader,"id2");

  double jet_pt_cut, jet_eta_cut;
  unsigned njets_required;

  size_t num_events = 0, ncount_total = 0;

  vector<fj::PseudoJet> partons;
  Higgs2diphoton Hdecay;
  std::pair<TLorentzVector,TLorentzVector> diphoton;
  boost::optional<TLorentzVector> higgs;

  ivanp::binner<bin_t, std::tuple<ivanp::axis_spec<
    ivanp::uniform_axis<int>,0,1>>> h_Njets;

public:
  diphoton_analysis_base(analysis_args&& args)
  : analysis_base(std::move(args)),
    _id1(reader,"id1"), _id2(reader,"id2")
  {
    const auto& jet_cuts = conf["jets"]["cuts"];
    jet_pt_cut  = jet_cuts["pT"];
    jet_eta_cut = jet_cuts["eta"];
    njets_required = conf["jets"]["Nrequired"];
    h_Njets({njets_required+2u,0,int(njets_required+2)});
  }

  void event_loop();
  virtual void fill_hists() = 0;
};

void excl_labels(TH1* h, bool excl) {
  auto* ax = h->GetXaxis();
  for (int i=1, n=h->GetNbinsX(); i<=n; ++i)
    ax->SetBinLabel(i,cat(excl ? "=" : ">=", i-1).c_str());
}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

TLorentzVector operator+(const TLorentzVector& a, const fj::PseudoJet& b) {
  return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] };
}
TLorentzVector& operator+=(TLorentzVector& a, const fj::PseudoJet& b) {
  a[0] += b[0]; a[1] += b[1]; a[2] += b[2]; a[3] += b[3]; return a;
}

void diphoton_analysis_base::event_loop() {
  // Reset ----------------------------------------------------------
  const size_t np = *_nparticle;
  partons.clear();
  higgs = boost::none;

  bin_t::weights = weights; // set weights

  // Keep track of multi-entry events -------------------------------
  nlo_bin::current_id = *_id;
  const bool new_id = (prev_id != nlo_bin::current_id);
  if (new_id) {
    prev_id = nlo_bin::current_id;
    ncount_total += ( _ncount ? **_ncount : 1);
    ++num_events;
  }

  // Read particles -------------------------------------------------
  unsigned n22 = 0; // number of photons
  for (size_t i=0; i<np; ++i) {
    if (_kf[i] == 25) {
      if (higgs) throw std::runtime_error("more than one Higgs");
      higgs.emplace(_px[i],_py[i],_pz[i],_E[i]);
    } else if (_kf[i] == 22) {
      if (n22>=2) throw std::runtime_error("more than two photons");
      (n22 ? diphoton.second : diphoton.first)
        .SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
      ++n22;
    } else {
      partons.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
    }
  }
  if (!higgs && n22!=2) throw std::runtime_error("missing Higgs or photons");

  cat_bin::id<isp>() = (unsigned)get_isp(*_id1,*_id2);
  // ----------------------------------------------------------------

  // Jets -----------------------------------------------------------
  auto fj_seq = fj::ClusterSequence(partons,jet_def);
  auto fj_jets = fj_seq.inclusive_jets(jet_pt_cut); // apply pT cut
  // apply eta cut
  for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
    if (std::abs(it->eta()) > jet_eta_cut) it = fj_jets.erase(it);
    else ++it;
  }
  // sort by pT
  std::sort( fj_jets.begin(), fj_jets.end(),
    [](const fj::PseudoJet& a, const fj::PseudoJet& b){
      return ( a.pt() > b.pt() );
    });
  // resulting number of jets
  const unsigned njets = fj_jets.size();

  // Higgs decay or add photons -------------------------------------
  if (higgs) diphoton = Hdecay(*higgs,new_id);
  else higgs = diphoton.first + diphoton.second;

  TLorentzVector *A1 = &diphoton.first, *A2 = &diphoton.second;
  double A1_pT = A1->Pt(), A2_pT = A2->Pt();
  if (A1_pT < A2_pT) {
    std::swap(A1,A2);
    std::swap(A1_pT,A2_pT);
  }

  // Photon cuts ----------------------------------------------------
  const double A1_eta = A1->Eta();
  const double A2_eta = A2->Eta();

  cat_bin::id<photon_cuts>()
    =  (A1_pT < 0.35*125.)
    && (A2_pT < 0.25*125.)
    && photon_eta_cut(std::abs(A1_eta))
    && photon_eta_cut(std::abs(A2_eta));

  // Fill Histograms ------------------------------------------------
  h_Njets.fill_bin(njets);

  // kinematic distributions must have at least N jets
  if (njets < njets_required) continue;

  this->fill_hists();
}

void diphoton_analysis_base::write_output() {
  cout << "Processed events: " << num_events << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(ofname,"recreate");
  if (fout->IsZombie()) return 1;
  TDirectory *dir = fout.get();

  auto h_Njets_integrated = h_Njets;
  h_Njets_integrated.integrate_left();

#define CATEGORY_TOP(r, data, elem) \
  cat_bin::id<elem>() = 0; \
  for (const char* dir_name : enum_traits<elem>::all_str()) { \
    dir = dir->mkdir(dir_name);

#define CATEGORY_BOT(r, data, elem) \
    dir = dir->GetMotherDir(); \
    ++cat_bin::id<elem>(); \
  }

  // write root historgrams
  wi = 0;
  for (const auto& _w : _weights) {
    dir = dir->mkdir(_w.GetBranchName());
    cout << dir->GetName() << endl;

    BOOST_PP_SEQ_FOR_EACH(CATEGORY_TOP,,CATEGORIES)

      dir->cd();

      using ivanp::root::to_root;
      using ivanp::root::slice_to_root;

      auto* h_Njets_excl = to_root(h_Njets,"Njets_excl");
      excl_labels(h_Njets_excl,true);
      auto* h_Njets_incl = to_root(h_Njets_integrated,"Njets_incl");
      h_Njets_incl->SetEntries( h_Njets_excl->GetEntries() );
      excl_labels(h_Njets_incl,false);

      for (auto& h : hist<1>::all) to_root(*h,h.name);
      for (auto& h : hist<1,0>::all) {
        const auto vars = ivanp::rsplit<1>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1]);
      }
      for (auto& h : hist<1,0,0>::all) {
        const auto vars = ivanp::rsplit<2>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1],vars[2]);
      }

    BOOST_PP_SEQ_FOR_EACH(CATEGORY_BOT,,BOOST_PP_SEQ_REVERSE(CATEGORIES))

    dir = dir->GetMotherDir();
    ++wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount_total);
  h_N->SetEntries(num_events);

  fout->Write();

  cout << "\n\033[32mWritten output file: " << fout->GetName() << endl;
}

