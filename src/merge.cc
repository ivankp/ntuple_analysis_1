#include <iostream>
#include <fstream>
#include <vector>
#include <map>
// #include <functional>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/regex.hpp>

#include "ivanp/error.hh"
#include "ivanp/tuple.hh"
#include "ivanp/container.hh"
#include "ivanp/functional.hh"
#include "ivanp/program_options.hh"
#include "ivanp/io/mem_file.hh"
#include "ivanp/scribe.hh"

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
using namespace ivanp;

inline json::json_pointer operator "" _jp(const char* s, size_t n) {
  return json::json_pointer(std::string(s,n));
}

template <typename T>
void compat(const T& key, const json& a, const json& b) {
  if (a.at(key) != b.at(key)) throw ivanp::error("Incompatible ",key);
}
template <typename T>
void compat(std::initializer_list<T> keys, const json& a, const json& b) {
  for (const auto& key : keys) compat(key,a,b);
}

template <typename T>
inline void add(scribe::value_node& a, const scribe::value_node& b) noexcept {
  // TEST(a.ptr())
  // TEST(reinterpret_cast<void*>(&a.cast<T&>()))
  // TEST(a.cast<T&>())
  a.cast<T&>() += b.cast<T>();
}

template <typename A, typename B, typename F>
inline void together(A& a, B& b, F&& f) {
  auto it_a = begin(a);
  auto it_b = begin(b);
  const auto end_a = end(a);
  for (; it_a!=end_a; ++it_a, ++it_b) f(*it_a,*it_b);
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname;
  bool merge_xsec = false,
       merge_variations = false,
       remove_blank = false;
  int verbose = 0;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input files",req(),pos())
      (ofname,'o',"output file",req())
      (merge_xsec,{"-x","--xsec","--nlo"},
       "merge cross sections (e.g. NLO parts)")
      (merge_variations,{"-u","--unc"},"merge scale & pdf variations")
      (remove_blank,{"-b","--rm-blank"},"remove blank histograms")
      (verbose,'v',"verbose [0,1,2]",switch_init(1))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  mem_file first_file;
  scribe::reader first, in;

  for (const char* ifname : ifnames) {
    // TEST(ifname)
    mem_file file = (
      ends_with(ifname,".xz") || ends_with(ifname,".lzma")
      ? mem_file::pipe(cat("unxz -c ",ifname).c_str())
      : mem_file::read(ifname)
    );

    in = { file.mem(), file.size() };
    auto& head = in.head();
    auto& info = head.at("info");
    auto& runcard = info.at("runcard");
    runcard.erase("entry_range");
    runcard.erase("input");
    runcard.erase("output");

    const auto offset = [m = file.mem()](void* p) -> auto& {
      return cout << "\033[32m" << (void*)((char*)p-m) << "\033[0m";
    };
    // offset(in.ptr());
    // TEST(in.memlen())
    // TEST(file.size())
    // TEST(in.head_str().size())
    // TEST(in.memlen() + in.head_str().size())

    if (!first) {
      first_file = std::move(file);
      first = std::move(in);
    } else {
      try {
        compat({"/info/runcard"_jp,"/root"_jp,"/types"_jp},first.head(),head);

        // TODO: merge counts

        // cout << "***********" << endl;
        together(first,in,[&](auto first, auto in){
          if (first["axes"] != in["axes"]) throw error(
            "histograms \"","!!TODO!!","\" have different axes"
          );
          // TEST(in.type_name())
          // offset(in.ptr());
          // TEST(first.memlen())
          y_combinator([&](auto f, auto first, auto in) -> void {
            const auto type = first.get_type();
            // TEST(type.name())
            offset(in.ptr()) << " " << type.name() << endl;
            together(first,in,[&](auto first, auto in) -> void {
              const auto type = first.get_type();
              offset(in.ptr()) << " " << type.name()
                << " " << type.is_union();
              if (type.is_union()) { cout
                << " " << (unsigned)in.union_index()
                << " " << (*in).get_type().is_null()
              ;}
              cout << endl;
              const char* type_name = type.name();
              // TEST(type.name())
              // offset(in.ptr());
              if (type.is_null()) return;
              if (type.is_union()) {
                if (first.union_index() != in.union_index())
                  throw error("union index mismatch");
                first = *first;
                if (first.get_type().is_null()) return;
                f(first,*in);
              }
              else if (!type.is_fundamental()) f(first,in);
              else {
                if (strlen(type_name)==2) { // add values
                  const char t = type_name[0], s = type_name[1];
                  if (t=='f') {
                    if (s=='8') { add<double  >(first,in); return; }
                    if (s=='4') { add<float   >(first,in); return; }
                  } else if (t=='u') {
                    if (s=='8') { add<uint64_t>(first,in); return; }
                    if (s=='4') { add<uint32_t>(first,in); return; }
                    if (s=='2') { add<uint16_t>(first,in); return; }
                    if (s=='1') { add<uint8_t >(first,in); return; }
                  } else if (t=='i') {
                    if (s=='8') { add<int64_t >(first,in); return; }
                    if (s=='4') { add<int32_t >(first,in); return; }
                    if (s=='2') { add<int16_t >(first,in); return; }
                    if (s=='1') { add<int8_t  >(first,in); return; }
                  }
                }
                throw error("unknown type \"",type_name,"\"");
              }
            });
          })(first["bins"],in["bins"]);
        });
      } catch (const std::exception& e) {
        cerr << "while reading file \"" << ifname << "\":\n" << e << endl;
        return 1;
      }
    }
  }

}
