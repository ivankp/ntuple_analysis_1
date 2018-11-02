#include <iostream>
#include <fstream>
#include <vector>
#include <array>

#include "ivanp/io/mem_file.hh"
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

  const mem_file file = (
    ends_with(argv[1],".xz") || ends_with(argv[1],".lzma")
    ? mem_file::pipe(cat("unxz -c ",argv[1]).c_str())
    : mem_file::mmap(argv[1])
  );

  ivanp::scribe::reader sr(file.mem(),file.size());
  TEST(sr.head_str())
  const json& head = sr.head();
  json sel;
  std::ifstream(argv[2]) >> sel;
  json out;
  auto& hist = out["hist"];
  auto& obins = hist["bins"] = json::array();

  std::vector<std::string> weights;

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
          break;
        }
        bin = bin[sel.value(type_name,"all")];
      }
    }
  })( sr[sel.at("hist")]["bins"] );

  for (const auto& w : weights)
    cout << w << endl;

  cout << out << endl;
}
