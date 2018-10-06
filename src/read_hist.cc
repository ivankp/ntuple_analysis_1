#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include "ivanp/scribe.hh"
#include "ivanp/functional.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;

using ivanp::scribe::size_type;

int main(int argc, char* argv[]) {
  ivanp::scribe::reader sr(argv[1]);
  TEST(sr.head())
  TEST(((const void*)sr.head().data()))
  auto head = nlohmann::json::parse(sr.head());
  nlohmann::json sel;
  std::ifstream(argv[2]) >> sel;

  TEST(sel["hist"])
  TEST(sr.ptr())
  auto hist = sr[sel["hist"]];
  TEST(hist.type_name())
  const auto type = hist.get_type();
  TEST(type.name())
  /*
  ivanp::y_combinator([](auto f, const auto t, auto name, int indent=0) -> void {
    for (int i=indent; i; --i) cout << "  ";
    cout << t.name();
    if (!name.empty()) cout << " : " << name;
    cout << endl;
    for (auto t : t) f(t.type,t.name,indent+1);
  })(type,std::string{});
  */

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
