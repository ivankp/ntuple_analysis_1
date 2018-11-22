#include <iostream>
#include <fstream>
#include <vector>

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

template <typename T>
inline T& get(json& x) { return x.get_ref<T&>(); }
template <typename T>
inline const T& get(const json& x) { return x.get_ref<const T&>(); }

int main(int argc, char* argv[]) {
  if (argc!=3 && argc!=2) {
    cout << "usage: " << argv[0] << " file.dat[.xz] [selection.json]\n";
    return 1;
  }

  const mem_file file = (
    ends_with(argv[1],".xz") || ends_with(argv[1],".lzma")
    ? mem_file::pipe(cat("unxz -c ",argv[1]).c_str())
    : mem_file::mmap(argv[1])
  );

  ivanp::scribe::reader sr(file.mem(),file.size());
  const json& head = sr.head();
  // TEST(sr.head_str())

  json sel, out;
  if (argc==3) {
    std::ifstream(argv[2]) >> sel;
  } else {
    std::ios::sync_with_stdio(false);
    std::cin >> sel;
  }
  const auto& hist_in = sr[sel.at("hist")];
  auto& hist = out["hists"][hist_in.get_name()];

  auto& axes = hist["axes"];
  for (const auto& axis : hist_in["axes"]) {
    json axis_obj;
    for (const auto& x : *axis)
      axis_obj[x.get_name()] = x;
    axes.push_back(axis_obj);
  }

  auto& obins = hist["bins"];
  auto& hvals = hist["values"];
  y_combinator([&](auto f, const auto& bins) -> void {
    for (auto bin : bins) {
      y_combinator([&](auto g, const auto& bin) -> void {
        const auto type = bin.get_type();
        // TEST(type.name())
        if (type.is_null()) obins.push_back(nullptr);
        else if (type.is_union()) g(*bin);
        else {
          const char* type_name = type.name();
          auto& name = sel[type_name];
          if (!strcmp(type_name,"weights")) {
            if (!name.is_null()) {
              obins.push_back(bin[get<std::string>(name)]);
            } else {
              const auto& w = bin[0];
              name = w.get_name();
              obins.push_back(w);
            }
            if (hvals.empty()) {
              const auto& w_type = bin[get<std::string>(name)].get_type();
              // TEST(w_type.name())
              for (const auto& weight_type : w_type)
                hvals.push_back(weight_type.name);
            }
          } else {
            if (!name.is_null()) g(bin[get<std::string>(name)]);
            // else f(bin); // loop over all
            else {
              const auto& b = bin[0];
              name = b.get_name();
              g(b);
            }
          }
        }
      })(bin);
    }
  })( hist_in["bins"] );

  auto& cats = hist["categories"];
  auto& cats_sel = hist["selection"];
  y_combinator([&](auto f, const auto& type, bool last=false) -> void {
    const char* type_name = type.name();
    auto& cat = cats[type_name];
    for (const auto& child : type) {
      if (sel[type_name]==child.name)
        cats_sel[type_name] = cat.size();
      // cats_sel.push_back(cat.size());
      cat.push_back(child.name);
    }
    auto next_type = type[0];
    if (!last) f(next_type,starts_with(next_type.name(),"weights"));
  })(hist_in.get_type().find("bins")[0][1]);

  auto& info = out["info"] = head.at("info");
  try {
    info.at("runcard").at("analysis").erase("binning");
    info.at("runcard").erase("reweighting");
    info.at("runcard").erase("entry_range");
  } catch (...) { }

  cout << out << endl;
}
