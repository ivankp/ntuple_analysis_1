#include <iostream>
#include <fstream>
#include <functional>

#include "json/nlohmann.hpp"

#include "ivanp/error.hh"
#include "ivanp/tuple.hh"
#include "ivanp/container.hh"
#include "lzma_compress.hh"
#include "ivanp/debug/type_str.hh"

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
  if (a.at(key) != b.at(key)) throw ivanp::error("Incompatible \"",key,"\"");
}

template <typename T>
void add(const std::tuple<json*,json*>& j) {
  *get<0>(j) = get<0>(j)->get<T>() + get<1>(j)->get<T>();
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
      file >> in;
    } catch (const std::exception& e) {
      cerr << "\033[31mError reading file\033[0m \"" << argv[i] << "\": "
           << e.what() << endl;
      return 1;
    }
    if (i==2) {
      auto& output = in.at("/annotation/runcard/output"_jp);
      output = { output };
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

      // const auto bins_sizes = out.at("/annotation/bins"_jp)
      //   | [](const auto& x){ return x.at(1).size(); };

      auto hists = tie(in,out) | [](auto& x){ return x.at("histograms"); };
      for (auto it=get<0>(hists).begin(),
               end=get<0>(hists).end(); it!=end; ++it)
      {
        auto& h_out = get<1>(hists).at(it.key());
        auto& h_in  = it.value();
        compat("axes",h_out,h_in);

        // auto bins = tie(h_out,h_in) | [](auto& x){ return &x.at("bins"); };
        /*{
          auto impl = [&](auto& f) mutable {
            const auto prev = bins;
            for (unsigned i=0, n=get<0>(bins)->size(); i<n; ++i) {
              bins = prev | [i](auto* x){ return &x->at(i); };
              if (get<0>(bins)->is_array()) f(f);
              else get<0>(bins) += get<1>(bins);
            }
            bins = prev;
          };
          impl(impl);
        }*/
        make_Ycombinator([](auto rec, const auto& bins) -> void {
          switch (get<0>(bins)->type()) {
            case (json::value_t::array):
              for (unsigned i=0, n=get<0>(bins)->size(); i<n; ++i)
                rec(bins | [i](auto* x){ return &x->at(i); });
              break;
            case (json::value_t::number_float): add<double>(bins); break;
            case (json::value_t::number_integer): add<long int>(bins); break;
            case (json::value_t::number_unsigned): add<long unsigned>(bins); break;
            case (json::value_t::null): break;
            default: throw ivanp::error(
              "non-numeric type of bin value \"",*get<0>(bins),"\"");
          }
        })(
          tie(h_out,h_in) | [](auto& x){ return &x.at("bins"); }
        );
      }
    }
  }

  if (ivanp::ends_with(argv[1],".xz")) {
    std::ofstream(argv[1]) << lzma_compress(out.dump(),6);
  } else std::ofstream(argv[1]) << out;
}
