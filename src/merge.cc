#include <iostream>
#include <fstream>
#include <functional>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>

#include "json/nlohmann.hpp"

#include "ivanp/error.hh"
#include "ivanp/tuple.hh"
#include "ivanp/container.hh"

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
  if (argc<2 || std::any_of(argv+1,argv+argc,[](const char* arg){
    return !strcmp(arg,"-h") || !strcmp(arg,"--help");
  })) {
    cout << "usage: " << argv[0] << " output.json input.json ..." << endl;
    return 1;
  }

  json out;
  for (int i=2; i<argc; ++i) {
    json in;
    try {
      std::ifstream file(argv[i]);
      file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
      if (ivanp::ends_with(argv[i],".xz")) {
        namespace bio = boost::iostreams;
        bio::filtering_streambuf<bio::input> buf;
        buf.push(bio::lzma_decompressor());
        buf.push(file);
        std::istream(&buf) >> in;
      } else file >> in;
    } catch (const std::exception& e) {
      cerr << "\033[31mError reading file\033[0m \"" << argv[i] << "\": "
           << e.what() << endl;
      return 1;
    }
    if (i==2) { // first file
      auto& output = in.at("/annotation/runcard/output"_jp);
      output = { output };
      in.at("/annotation/runcard"_jp).erase("entry_range");
      out = in;
    } else {
      try {
        compat("/annotation/bins"_jp,out,in);
        compat("/annotation/runcard/analysis"_jp,out,in);
      } catch (const std::exception& e) {
        cerr << "In file \"" << argv[i] << "\": " << e << endl;
        return 1;
      }

      auto& input_in  = in .at("/annotation/runcard/input"_jp);
      auto& input_out = out.at("/annotation/runcard/input"_jp);
      for (unsigned i=0, n=input_in.size(); i<n; ++i)
        for (auto&& f : input_in[i].at("files"))
          input_out[i].at("files").emplace_back(std::move(f));

      out.at("/annotation/runcard/output"_jp).emplace_back(std::move(
        in.at("/annotation/runcard/output"_jp)));

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

  try {
    std::ofstream file(argv[1]);
    file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    if (ivanp::ends_with(argv[1],".xz")) {
      namespace bio = boost::iostreams;
      bio::filtering_streambuf<bio::output> buf;
      buf.push(bio::lzma_compressor(bio::lzma::best_compression));
      buf.push(file);
      std::ostream(&buf) << out;
    } else file << out;
  } catch (const std::exception& e) {
    cerr << "\033[31mError writing file\033[0m \"" << argv[1] << "\": "
         << e.what() << endl;
    return 1;
  }
}
