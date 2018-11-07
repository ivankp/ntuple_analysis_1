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

  // std::vector<std::string> weights;

  y_combinator([&](auto f, const auto& bins) -> void {
    for (auto bin : bins) {
      y_combinator([&](auto g, const auto& bin) -> void {
        const auto type = bin.get_type();
        TEST(type.name())
        if (type.is_null()) obins.push_back(nullptr);
        else if (type.is_union()) {
          TEST((unsigned)bin.union_index())
          g(*bin);
        }
        else {
          const char* type_name = type.name();
          if (starts_with(type_name,"nlo_bin<")) {
            obins.push_back(bin);
          } else {
            auto it = sel.find(type_name);
            if (it==sel.end()) f(bin);
            else g(bin[it->get<std::string>()]);
          }
        }
      })(bin);
    }
  })( sr[sel.at("hist")]["bins"] );

  // for (const auto& w : weights)
  //   cout << w << endl;

  cout << out << endl;
}
