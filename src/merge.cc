#include <iostream>
#include <fstream>
#include <vector>
#include <map>

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

bool merge_xsec = false,
     merge_variations = false,
     remove_blank = false;
int verbose = 0;

using count_t = long unsigned;

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
  a.cast<T&>() += b.cast<T>();
}

template <typename T>
inline T& as(const scribe::value_node& a) noexcept { return a.cast<T&>(); }

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
  if (verbose && merge_xsec)
    cout << "Will normalize to cross section" << endl;

  mem_file first_file;
  scribe::reader first;

  std::map<std::string,count_t> total_count;

  for (const char* ifname : ifnames) {
    if (verbose) cout << "> " << ifname << endl;
    if (!std::ifstream(ifname)) {
      cerr << "\033[31mError reading file\033[0m" << endl;
      return 1;
    }
    mem_file file = (
      ends_with(ifname,".xz") || ends_with(ifname,".lzma")
      ? mem_file::pipe(cat("unxz -c ",ifname).c_str())
      : mem_file::read(ifname)
    );

    scribe::reader in(file.mem(), file.size());
    auto& head = in.head();
    auto& info = head.at("info");
    auto& runcard = info.at("runcard");
    runcard.erase("entry_range");
    runcard.erase("input");
    runcard.erase("output");

    // add counts
    auto& count = info.at("count");
    for (auto it=count.begin(), end=count.end(); it!=end; ++it)
      total_count[it.key()] += it->get<count_t>();
    const auto norm = count.at("norm").get<count_t>();
    const auto norm2 = norm*norm;

    if (merge_xsec) { // normalize to cross section
      for (auto hist : in) {
        y_combinator([=](auto f, auto bins) -> void {
          for (auto bin : bins) {
            y_combinator([=](auto g, auto bin) -> void {
              const auto type = bin.get_type();
              if (type.is_null()) return;
              else if (type.is_union()) g(*bin);
              else if (starts_with(type.name(),"nlo_bin<")) {
                /*
                for (const auto& w : bin) { // loop over weights
                  if (strcmp(w.type_name(),"f8#2")) continue;
                  // w[0].cast<double&>() /= norm;
                  as<double>(w[0]) /= norm;
                  as<double>(w[1]) /= norm2;
                }
                */
                const unsigned n = bin.size()-1;
                double* w = &as<double&>(bin);
                for (unsigned i=0; i<n; ) {
                  w[i++] /= norm;
                  w[i++] /= norm2;
                }
              }
              else f(bin);
            })(bin);
          }
        })(hist["bins"]);
      }
    }

    if (!first) { // first file
      first_file = std::move(file);
      first = std::move(in);
    } else { // further files
      try {
        compat({"/info/runcard"_jp,"/root"_jp,"/types"_jp},first.head(),head);

        together(first,in,[](auto first, auto in){
          if (verbose>=2) cout << "  " << first.get_name() << endl;
          if (first["axes"] != in["axes"]) throw error(
            "histograms \"",first.get_name(),"\" have different axes"
          );
          y_combinator([](auto f, auto first, auto in) -> void {
            together(first,in,
            y_combinator([f](auto g, auto first, auto in) -> void {
              const auto type = first.get_type();
              if (type.is_null()) return;
              else if (type.is_union()) {
                if (first.union_index() != in.union_index())
                  throw error("union type mismatch");
                g(*first,*in);
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

  auto& head = first.head();

  if (merge_xsec) total_count["norm"] = 1;
  head["/info/count"_jp] = total_count;

  char* out = first.data_ptr();
  size_t out_len = first.data_len();
  bool new_out = false;

  if (merge_variations) { // merge variations
    if (verbose) cout << "Merging variations" << endl;

    struct w_struct { string pdf,ren,fac; vector<unsigned> scale_i,pdf_i; };
    std::map<string,w_struct> ws;
    vector<unsigned> other_weights;
    vector<string> other_names;
    auto& types = head.at("types");
    auto& nlo_bin_type = types.at("nlo_bin<f8#>");
    const auto& w_names = nlo_bin_type.at(0);
    const unsigned n_w_names = w_names.size()-1;
    {
      using namespace boost;
      const regex re("([^:]+):(\\d+)(?: ren:([\\d.]+))?(?: fac:([\\d.]+))?");
      smatch m;
      int i = -1;
      for (string&& w_name : w_names) {
        if (i < 0) { ++i; continue; }
        if (regex_match(w_name,m,re)) {
          auto& w = ws[m[1]];
          if (!w.scale_i.size()) {
            w.pdf = m[2];
            w.ren = m[3];
            w.fac = m[4];
            w.scale_i.push_back(i);
            w.  pdf_i.push_back(i);
          } else {
            if (m[2]==w.pdf && (m[3]!=w.ren || m[4]!=w.fac))
              w.scale_i.push_back(i);
            else if (m[2]!=w.pdf && m[3]==w.ren && m[4]==w.fac)
              w.pdf_i.push_back(i);
            else throw ivanp::error("unexpected weight \"",m[0],'\"');
          }
        } else {
          other_weights.push_back(i);
          other_names.emplace_back(std::move(w_name));
        }
        ++i;
      }
    }

    if (verbose>=2) {
      for (auto i : other_weights) cout << "  " << other_names[i] << '\n';
      for (const auto& w : ws) cout << "  " << w.first << '\n';
    }
    cout << std::flush;

    // Note: assumes that all input weights are represented by pairs of doubles
    // (8 byte floats) [ weight, sumw2 ]
    // followed by an 8 byte bin count

    types["envelope"] = R"(
      [ [ "f8#2", "central" ], [ "[null,f8#2]", "scale_unc", "pdf_unc" ] ]
    )"_json;
    // auto new_nlo_bin_type = R"(
    //   [ [ "f8#2" ], [ "envelope" ], [ "u8", "n" ] ]
    // )"_json;
    auto new_nlo_bin_type = R"(
      [ [ "f8#2" ], [ "u8", "n" ] ]
    )"_json;
    for (auto i : other_weights) new_nlo_bin_type[0].push_back(other_names[i]);
    // for (const auto& w : ws) new_nlo_bin_type[1].push_back(w.first);
    nlo_bin_type = std::move(new_nlo_bin_type);

    const unsigned old_size = n_w_names*8;
    const unsigned new_size = (other_weights.size()*8 + ws.size()*(8+(8+1)*2));
    const unsigned size_rat = (new_size + old_size - 1) / old_size;
    char* cur_old = out;
    if (size_rat > 1) {
      if (verbose>=2)
        cout << "\033[33moutput will be larger than input\033[0m" << endl;
      new_out = true;
      out = new char[out_len*size_rat];
    }
    char* cur_new = out;

    for (auto hist : first) {
      TEST(hist.get_name())
      y_combinator([&](auto f, auto bins) -> void {
        for (auto _bin : bins) {
          y_combinator([&](auto g, auto bin) -> void {
            const auto type = bin.get_type();
            if (type.is_null()) return;
            else if (type.is_union()) g(*bin);
            else if (starts_with(type.name(),"nlo_bin<")) {
              { // data between bins
                const auto len = bin.ptr()-cur_old;
                memmove(cur_new,cur_old,len);
                cur_old += len + bin.memlen();
                cur_new += len;
              }
              const unsigned nw = bin.size()-1;
              auto* w = reinterpret_cast<double(*)[2]>(bin.ptr());
              for (auto i : other_weights) { // not merged weights
                const auto len = sizeof(*w);
                memmove(cur_new, w+i, len);
                cur_new += len;
              }
              /*
              { // merged weights : TODO

              }
              */
              { // n entries
                const auto len = 8;
                memmove(cur_new, w+nw, len);
                cur_new += len;
              }
            }
            else f(bin);
          })(_bin);
        }
      })(hist["bins"]);
    }
    const auto len = (first.data_ptr()+first.data_len())-cur_old;
    TEST(len)
    if (len) { // data after bins
      memmove(cur_new,cur_old,len);
      cur_old += len;
      cur_new += len;
    }
    out_len = cur_new - out;
  }

  if (verbose) cout << "< " << ofname << std::flush;
  try { // write output file
    std::ofstream file(ofname);
    file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    if (ivanp::ends_with(ofname,".xz")) {
      namespace bio = boost::iostreams;
      bio::filtering_streambuf<bio::output> buf;
      buf.push(bio::lzma_compressor(bio::lzma::best_compression));
      buf.push(file);
      (std::ostream(&buf) << first.head()).write(out,out_len);
    } else {
      (file << first.head()).write(out,out_len);
    }
  } catch (const std::exception& e) {
    if (verbose) cout << " \033[31;1m✘\033[0m" << endl;
    cerr << "\033[31mError writing file\033[0m \"" << ofname << "\": "
         << e.what() << endl;
    return 1;
  }
  if (verbose) cout << " \033[32;1m✔\033[0m" << endl;

  if (new_out) delete[] out;
}
