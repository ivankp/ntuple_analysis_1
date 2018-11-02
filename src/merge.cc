#include <iostream>
#include <fstream>
#include <vector>
#include <map>
// #include <functional>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>
// #include <boost/regex.hpp>

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
  scribe::reader first;

  for (const char* ifname : ifnames) {
    cout << "> " << ifname << endl;
    mem_file file = (
      ends_with(ifname,".xz") || ends_with(ifname,".lzma")
      ? mem_file::pipe(cat("unxz -c ",ifname).c_str())
      : mem_file::read(ifname)
      // : (!first ? mem_file::read(ifname) : mem_file::mmap(ifname))
    );

    scribe::reader in(file.mem(), file.size());
    auto& head = in.head();
    auto& info = head.at("info");
    auto& runcard = info.at("runcard");
    runcard.erase("entry_range");
    runcard.erase("input");
    runcard.erase("output");

    // const auto offset = [m = file.mem()](void* p) -> auto& {
    //   return cout << "\033[32m" << (void*)((char*)p-m) << "\033[0m";
    // };

    if (!first) {
      first_file = std::move(file);
      first = std::move(in);
    } else {
      try {
        compat({"/info/runcard"_jp,"/root"_jp,"/types"_jp},first.head(),head);

        // TODO: merge counts

        together(first,in,[](auto first, auto in){
          if (first["axes"] != in["axes"]) throw error(
            "histograms \"","!!TODO!!","\" have different axes"
          );
          y_combinator([](auto f, auto first, auto in) -> void {
            together(first,in,
            y_combinator([&f](auto g, auto first, auto in) -> void {
              const auto type = first.get_type();
              if (type.is_null()) return;
              else if (type.is_union()) {
                if (first.union_index() != in.union_index())
                  throw error("union index mismatch");
                g(*first,*in);
                // do {
                //   first = *first;
                //   in = *in;
                //   type = first.get_type();
                // } while (type.is_union());
                // if (!type.is_null()) f(first,in);
              }
              else if (type.is_fundamental()) { // add values
                const char* type_name = type.name();
                const char t = type_name[0], s = type_name[1];
                if (t=='f') {
                  if (s=='8') add<double  >(first,in); else
                  if (s=='4') add<float   >(first,in);
                } else if (t=='u') {
                  if (s=='8') add<uint64_t>(first,in); else
                  if (s=='4') add<uint32_t>(first,in); else
                  if (s=='2') add<uint16_t>(first,in); else
                  if (s=='1') add<uint8_t >(first,in);
                } else if (t=='i') {
                  if (s=='8') add<int64_t >(first,in); else
                  if (s=='4') add<int32_t >(first,in); else
                  if (s=='2') add<int16_t >(first,in); else
                  if (s=='1') add<int8_t  >(first,in);
                }
              }
              else f(first,in);
            }));
          })(first["bins"],in["bins"]);
        });
      } catch (const std::exception& e) {
        cerr << "while reading file \"" << ifname << "\":\n" << e << endl;
        return 1;
      }
    }
  } // end file loop

  cout << "< " << ofname << std::flush;
  try { // write output file
    std::ofstream file(ofname);
    file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    if (ivanp::ends_with(ofname,".xz")) {
      namespace bio = boost::iostreams;
      bio::filtering_streambuf<bio::output> buf;
      buf.push(bio::lzma_compressor(bio::lzma::best_compression));
      buf.push(file);
      (std::ostream(&buf) << first.head())
        .write(first.data_ptr(),first.data_len());
    } else {
      (file << first.head()).write(first.data_ptr(),first.data_len());
    }
  } catch (const std::exception& e) {
    cout << " \033[31;1m✘\033[0m" << endl;
    cerr << "\033[31mError writing file\033[0m \"" << ofname << "\": "
         << e.what() << endl;
    return 1;
  }
  cout << " \033[32;1m✔\033[0m" << endl;
}
