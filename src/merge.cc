#include <iostream>
#include <fstream>
#include <functional>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>

#include "json/nlohmann.hpp"

#include "ivanp/error.hh"
#include "ivanp/tuple.hh"
#include "ivanp/container.hh"
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

// https://stackoverflow.com/a/40873657/2640636
template <class F> struct Ycombinator {
  F f;
  template <class... Args>
  decltype(auto) operator()(Args&&... args) const {
    return f(std::ref(*this), std::forward<Args>(args)...);
  }
};
template <class F>
Ycombinator<std::decay_t<F>> make_Ycombinator(F&& f) {
  return { std::forward<F>(f) };
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname;
  bool merge_xsec = false;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input JSON files",req(),pos())
      (ofname,'o',"output JSON file",req())
      (merge_xsec,{"-x","--xsec","--nlo"},
       "merge cross sections (e.g. NLO parts)")
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  json out;
  for (const char* ifname : ifnames) {
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
      for (auto& h : in.at("histograms")) {
        make_Ycombinator([norm](auto rec, auto& bin, unsigned depth) -> void {
          if (bin.is_null()) return;
          if (depth) for (auto& b : bin) rec(b,depth-1);
          else {
            bin[0] = get<double>(bin.at(0))/norm;
            bin[1] = get<double>(bin.at(1))/(norm*norm);
          }
        })(h.at("bins"),in.at("/annotation/bins"_jp).size());
      }
    }

    if (first) { // first file
      first = false;
      out = std::move(in);
      auto& output = out.at("/annotation/runcard/output"_jp);
      output = { output };
      out.at("/annotation/runcard"_jp).erase("entry_range");
    } else {
      try {
        compat("/annotation/bins"_jp,out,in);
        compat("/annotation/runcard/analysis"_jp,out,in);
      } catch (const std::exception& e) {
        cerr << "In file \"" << ifname << "\": " << e << endl;
        return 1;
      }

      auto& input_in  = in .at("/annotation/runcard/input"_jp);
      auto& input_out = out.at("/annotation/runcard/input"_jp);
      for (unsigned i=0, n=input_in.size(); i<n; ++i)
        for (auto&& f : input_in[i].at("files"))
          input_out[i].at("files").emplace_back(std::move(f));

      out.at("/annotation/runcard/output"_jp).emplace_back(std::move(
        in.at("/annotation/runcard/output"_jp)));

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
        auto& h_out = it.value();
        auto& h_in  = get<1>(hists)->at(it.key());
        compat("axes",h_out,h_in);

        make_Ycombinator([](auto rec, const auto& bins) -> void {
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
