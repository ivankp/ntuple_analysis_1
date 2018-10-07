#include <iostream>
#include <fstream>
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

  /*
  auto hist = sr[sel["hist"]];
  TEST(hist.ptr())
  TEST(hist.type_name())
  auto bins = hist["bins"];
  TEST(bins.ptr())
  TEST(bins.type_name())
  TEST(bins.size())
  // for (const auto& bin : bins) {
  //   TEST(bin.type_name())
  //   TEST((unsigned)bin.union_index())
  // }
  TEST(bins[1].type_name())
  TEST((unsigned)bins[1].union_index())
  TEST(bins[1][0].type_name())
  TEST(bins[1][0][0].type_name())
  TEST(bins[1][0][0]["n"].type_name())
  TEST(bins[1]["all"]["all"]["w"].type_name())
  TEST(bins[4]["all"]["all"]["w"][0][0].cast<double>())
  */

  /*
  for (auto bin : sr[sel.at("hist")]["bins"]) {
    if (bin.union_index()==0) continue;
    bin = *bin;
    for (;;) {
      const char* type_name = bin.type_name();
      TEST(type_name)
      if (starts_with(type_name,"nlo_bin<")) break;
      bin = bin[sel.at(type_name)];
    }
  }
  */

  y_combinator([&](auto f, const auto& bins) -> void {
    for (auto bin : bins) {
      const auto ui = bin.union_index();
      if (ui==0) continue;
      bin = *bin;
      if (ui==2) f(bin);
      for (;;) {
        const char* type_name = bin.type_name();
        TEST(type_name)
        if (starts_with(type_name,"nlo_bin<")) break;
        bin = bin[sel.at(type_name)];
      }
    }
  })( sr[sel.at("hist")]["bins"] );


  /*
  TEST(hist.ptr())
  TEST(type.num_children())
  TEST(hist["axes"].ptr())
  TEST(hist["axes"].type_name())

  TEST(type.begin()->name)
  TEST(type.begin()->type.name())
  TEST((type.begin()+1)->name)
  TEST((type.begin()+1)->type.name())
  TEST(type[1].name())

  TEST(type[0].is_array())
  TEST(type[0].is_union())
  TEST(type[0].memlen())
  TEST(type[0].size())
  TEST(type[1].is_array())
  TEST(type[1].is_union())
  TEST(type[1].memlen())
  TEST(type[1].size())

  TEST(hist[1].ptr())
  */
  // for (auto t : type) TEST(t->name())
  // for (auto t : type[1]) TEST(t->name())
  // for (auto t : type[1][0]) TEST(t->name())
  // for (auto t : type[1][0][1]) TEST(t->name())
  // TEST(hist["bins"].type_name())
  // TEST(hist["bins"].type_name(0))
}
