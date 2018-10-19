#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include "ivanp/scribe.hh"
#include "ivanp/scribe/json.hh"
#include "ivanp/functional.hh"
#include "ivanp/error.hh"
#include "ivanp/string.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;

using ivanp::scribe::size_type;
using namespace ivanp;
using nlohmann::json;

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " file.dat[.xz] selection.json\n";
    return 1;
  }
  ivanp::scribe::reader sr(argv[1]);
  TEST(sr.head_str())
  // TEST(((const void*)sr.head().data()))
  const json& head = sr.head();
  json sel;
  std::ifstream(argv[2]) >> sel;
  json out;
  auto& hist = out["hist"];
  auto& obins = hist["bins"] = json::array();

  std::vector<std::string> weights;
  /*
  struct bin_t {
    std::vector<std::vector<double>> ws;
    long unsigned n;
    bin_t(const double* w, size_type nw, long unsigned n): ws(w,w+nw), n(n) { }
  };
  std::vector<bin_t> hist;
  */

  y_combinator([&](auto f, const auto& bins) -> void {
    for (auto bin : bins) {
      const auto ui = bin.union_index();
      if (ui==0) {
        obins.push_back(nullptr);
        continue;
      }
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
          obins.push_back(bin);
          /*
          TEST(bin.type_name())
          TEST(weights.size())
          TEST(weights.size())
          hist.emplace_back(
            &bin[0].template cast<double>(), weights.size(),
            bin["n"].template cast<long unsigned>()
          );

          auto it = bin.begin();
          const decltype(it) end = weights.size();
          auto& ws = hist.back().w;
          for (; it!=end; ++it) { // loop over weights
            const auto val = *it;
            const auto* _w = &val.template cast<double>();
            for (unsigned i=WEIGHT_N_VAL; i; ) --i, w[i] = _w[i];
          }
          */
          break;
        }
        bin = bin[sel.value(type_name,"all")];
      }
    }
  })( sr[sel.at("hist")]["bins"] );

  for (const auto& w : weights)
    cout << w << endl;

  // for (const auto& bin : hist)
  //   cout << bin.ws[0] << ' ' << bin.ws[1] << ' ' << bin.n << endl;
  cout << out << endl;
}
