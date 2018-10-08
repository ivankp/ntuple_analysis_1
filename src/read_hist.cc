#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <nlohmann/json.hpp>
#include "ivanp/scribe.hh"
#include "ivanp/functional.hh"
#include "ivanp/error.hh"
#include "ivanp/string.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;

using ivanp::scribe::size_type;
using namespace ivanp;

int main(int argc, char* argv[]) {
  ivanp::scribe::reader sr(argv[1]);
  TEST(sr.head())
  // TEST(((const void*)sr.head().data()))
  auto head = nlohmann::json::parse(sr.head());
  nlohmann::json sel;
  std::ifstream(argv[2]) >> sel;

  std::vector<std::string> weights;
  struct bin_t {
    #define WEIGHT_N_VAL 2
    std::array<double,WEIGHT_N_VAL> w;
    long unsigned n;
    bin_t(long unsigned n): w(), n(n) { }
  };
  std::vector<bin_t> hist;

  y_combinator([&](auto f, const auto& bins) -> void {
    for (auto bin : bins) {
      const auto ui = bin.union_index();
      if (ui==0) continue;
      bin = *bin;
      if (ui==2) f(bin);
      for (;;) {
        const char* type_name = bin.type_name();
        // TEST(type_name)
        if (starts_with(type_name,"nlo_bin<")) {
          if (weights.empty()) {
            const auto& ws = head["types"][type_name][0];
            auto it = ++ws.begin();
            for (auto end=ws.end(); it!=end; ++it) {
              weights.emplace_back();
              it->get_to(weights.back());
            }
          }
          hist.emplace_back(bin["n"].template cast<long unsigned>());
          auto it = bin.begin();
          const decltype(it) end = weights.size();
          auto& w = hist.back().w;
          for (; it!=end; ++it) { // loop over weights
            const auto val = *it;
            if (val.size()!=WEIGHT_N_VAL) throw error(
              "more values per weight than expected");
            const auto* _w = &val.template cast<double>();
            for (unsigned i=WEIGHT_N_VAL; i; ) --i, w[i] = _w[i];
          }
          break;
        }
        bin = bin[sel.at(type_name)];
      }
    }
  })( sr[sel.at("hist")]["bins"] );

  for (const auto& w : weights)
    cout << w << endl;

  for (const auto& bin : hist)
    cout << bin.w[0] << ' ' << bin.w[1] << ' ' << bin.n << endl;
}
