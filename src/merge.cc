#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <functional>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/regex.hpp>

#include <nlohmann/json.hpp>

#include "ivanp/error.hh"
#include "ivanp/tuple.hh"
#include "ivanp/container.hh"
#include "ivanp/functional.hh"
#include "ivanp/program_options.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::get;
using std::tie;
using nlohmann::json;
using ivanp::cat;
using ivanp::y_combinator;

inline json::json_pointer operator "" _jp(const char* s, size_t n) {
  return json::json_pointer(std::string(s,n));
}

template <typename T>
void compat(const T& key, const json& a, const json& b) {
  if (a.at(key) != b.at(key)) throw ivanp::error("Incompatible ",key);
}

template <typename T>
void add(const std::tuple<json*,json*>& j) {
  *get<0>(j) = ( get<0>(j)->get<T>() + get<1>(j)->get<T>() );
}
template <typename T> T get(const json& j) { return j.get<T>(); }

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname;
  bool merge_xsec = false,
       merge_variations = false,
       remove_blank = false;
  int verbose = 0;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input JSON files",req(),pos())
      (ofname,'o',"output JSON file",req())
      (merge_xsec,{"-x","--xsec","--nlo"},
       "merge cross sections (e.g. NLO parts)")
      (merge_variations,{"-u","--unc"},"merge scale & pdf variations")
      (remove_blank,{"-b","--rm-blank"},"remove blank histograms")
      (verbose,'v',"verbose [0,1,2]",switch_init(1))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  json out;
  for (const char* ifname : ifnames) {
    if (verbose>=1) cout << ifname << endl;
    static bool first = true;
    json in;
    try {
      std::ifstream file(ifname);
      file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
      if (ivanp::ends_with(ifname,".xz")) {
        namespace bio = boost::iostreams;
        bio::filtering_streambuf<bio::input> buf;
        buf.push(bio::lzma_decompressor());
        buf.push(file);
        std::istream(&buf) >> in;
      } else file >> in;
    } catch (const std::exception& e) {
      cerr << "\033[31mError reading file\033[0m \"" << ifname << "\": "
           << e.what() << endl;
      return 1;
    }

    if (merge_xsec) { // scale input to scross section
      const double norm = in.at("/annotation/count/norm"_jp);
      if (norm!=1) {
        const auto bins_depth = in.at("/annotation/bins"_jp).size();
        for (auto& h : in.at("histograms")) {
          y_combinator([norm](auto rec, auto& bin, unsigned depth)->void {
            if (bin.is_null()) return;
            if (depth) for (auto& b : bin) rec(b,depth-1);
            else for (auto& b : bin[0]) {
              b[0] = get<double>(b.at(0))/norm;
              b[1] = get<double>(b.at(1))/(norm*norm);
            }
          })(h.at("bins"),bins_depth-1+h.at("axes").size());
        }
      }
    }

    if (first) { // first file
      first = false;
      out = std::move(in);
      auto& runcard = out.at("/annotation/runcard"_jp);
      runcard.erase("entry_range");
      runcard.erase("output");
      if (verbose>=2) {
        auto& hists = out.at("histograms");
        for (auto h=hists.begin(), end=hists.end(); h!=end; ++h)
          cout << "  " << h.key() << endl;
      }
    } else {
      try {
        compat("/annotation/bins"_jp,out,in);
        compat("/annotation/runcard/analysis"_jp,out,in);
        compat("/annotation/runcard/reweighting"_jp,out,in);
      } catch (const std::exception& e) {
        cerr << "In file \"" << ifname << "\": " << e << endl;
        return 1;
      }

      auto& input_in  = in .at("/annotation/runcard/input"_jp);
      auto& input_out = out.at("/annotation/runcard/input"_jp);
      for (unsigned i=0, n=input_in.size(); i<n; ++i)
        for (auto&& f : input_in[i].at("files"))
          input_out[i].at("files").emplace_back(std::move(f));

      auto count = tie(out,in) | [](auto& x){
        return &x.at("/annotation/count"_jp);
      };
      for (auto it=get<0>(count)->begin(),
               end=get<0>(count)->end(); it!=end; ++it)
      {
        add<long unsigned>(
          std::make_tuple(&it.value(),&get<1>(count)->at(it.key())));
      }

      auto hists = tie(out,in) | [](auto& x){ return &x.at("histograms"); };
      for (auto it=get<0>(hists)->begin(),
               end=get<0>(hists)->end(); it!=end; ++it)
      {
        if (verbose>=2) cout << "  " << it.key() << endl;
        auto& h_out = it.value();
        auto& h_in  = get<1>(hists)->at(it.key());
        compat("axes",h_out,h_in);

        y_combinator([](auto rec, const auto& bins) -> void {
          switch (get<0>(bins)->type()) {
            case (json::value_t::array):
              for (unsigned i=0, n=get<0>(bins)->size(); i<n; ++i)
                rec(bins | [i](auto* x){ return &x->at(i); });
              break;
            case (json::value_t::number_float): add<double>(bins); break;
            case (json::value_t::number_unsigned):
            case (json::value_t::number_integer): add<long>(bins); break;
            case (json::value_t::null): break;
            default: throw ivanp::error(
              "non-numeric bin value \"",*get<0>(bins),"\"");
          }
        })(
          tie(h_out,h_in) | [](auto& x){ return &x.at("bins"); }
        );
      }
    }
  }

  if (merge_xsec) out.at("/annotation/count/norm"_jp) = 1;

  if (remove_blank) {
    auto& hists = out.at("histograms");
    const auto bins_depth = out.at("/annotation/bins"_jp).size();
    for (auto h=hists.begin(); h!=hists.end(); ) {
      if (
        !y_combinator([](auto f, const auto& bin, unsigned depth) -> bool {
          if (depth) {
            for (const auto& b : bin) if (f(b,depth-1)) return true;
          } else {
            if (get<long>(bin[1])!=0) return true;
          }
          return false;
        })( h->at("bins"), bins_depth-1+(h->at("axes").size()) )
      ){
        if (verbose>=1) cout << "\033[31m✘\033[0m " << h.key() << endl;
        h = hists.erase(h);
      } else ++h;
    }
  }

  if (merge_variations) {
    if (verbose>=1) cout << "\033[36mmerging variations\033[0m" << endl;
    auto& weights = out.at("/annotation/weights"_jp);
    struct w_struct { string pdf,ren,fac; vector<unsigned> scale_i,pdf_i; };
    std::map<string,w_struct> ws;
    vector<unsigned> other_weights;
    using namespace boost;
    const regex re("([^:]+):(\\d+)(?: ren:([\\d.]+))?(?: fac:([\\d.]+))?");
    smatch m;
    for (unsigned i=0, n=weights.size(); i<n; ++i) {
      if (regex_match(weights[i].get<string>(),m,re)) {
        auto& w = ws[m[1]];
        if (!w.scale_i.size()) {
          w.pdf = m[2];
          w.ren = m[3];
          w.fac = m[4];
          w.scale_i.push_back(i);
          w.  pdf_i.push_back(i);
        } else {
          if (m[2]==w.pdf && (m[3]!=w.ren || m[4]!=w.fac))
            w.scale_i.push_back(i);
          else if (m[2]!=w.pdf && m[3]==w.ren && m[4]==w.fac)
            w.pdf_i.push_back(i);
          else throw ivanp::error("unexpected weight \"",m[0],'\"');
        }
      } else other_weights.push_back(i);
    }
    auto weights2 = json::array();
    for (auto i : other_weights) weights2.push_back(weights[i]);
    for (const auto& w : ws) weights2.push_back(w.first);
    weights = weights2;

    auto& ann_bins = out.at("/annotation/bins"_jp);
    auto& bin_vars = std::find_if(ann_bins.begin(),ann_bins.end(),
      [](const auto& x){ return x[0]=="bin"; })->at(1).at(0);
    bin_vars.push_back({"scale",{"min","max"}});
    bin_vars.push_back({"pdf"  ,{"min","max"}});
    const auto bins_depth = ann_bins.size();
    auto& hists = out.at("histograms");
    for (auto h=hists.begin(), end=hists.end(); h!=end; ++h) {
      if (verbose>=2) cout << "  " << h.key() << endl;
      y_combinator([&](auto f, auto& bin, unsigned depth) -> void {
        if (depth) for (auto& b : bin) f(b,depth-1);
        else {
          auto& b = bin[0];
          auto b2 = json::array();
          for (auto i : other_weights) b2.push_back(b[i]);
          for (const auto& w : ws) {
            auto _b = b[w.second.scale_i[0]];
            std::tie(w.second.scale_i,w.second.pdf_i) | [&](const auto& is) {
              if (is.size()<2) _b.push_back(nullptr);
              else {
                _b.push_back(ivanp::minmax(
                  is | [&](auto i) -> double { return b[i][0]; }
                ));
              }
            };
            b2.push_back(_b);
          }
          b = b2;
        }
      })( h->at("bins"), bins_depth-1+(h->at("axes").size()) );
    }
  }

  if (verbose>=1) cout << "Writing: " << ofname << endl;
  try {
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
}
