#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstring>
#include <cstdio>

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << (var) << std::endl;

#include "ivanp/sqlite.hh"
#include "ivanp/program_options.hh"
#include "ivanp/string.hh"

using std::cout;
using std::endl;
using std::cerr;
using std::get;
using ivanp::cat;

template <typename C, typename T>
bool has(const C& cont, T&& x) {
  return std::find_if(cont.begin(),cont.end(),std::forward<T>(x))!=cont.end();
}

#define CAST(T,X) reinterpret_cast<T*&>(X)
#define ADDER_FCN(S,T) \
  else if (!strncmp(str,S,i)) { \
    fs.emplace_back([](void*& a, const void*& b) noexcept { \
      (*CAST(T,a)) += (*CAST(const T,b)); ++CAST(T,a); ++CAST(const T,b); \
    }); \
  }

static_assert(sizeof(float)==4);
static_assert(sizeof(double)==8);

class adder {
  typedef void (*f_t)(void*&, const void*&);
  std::vector<f_t> fs;
  struct error: std::runtime_error {
    using std::runtime_error::runtime_error;
  };
public:
  adder(const char* str) {
    unsigned i = 0;
    for (;;) {
      char c = str[i];
      if (c==',' || c=='\0') {
        if (!i) throw error("blank type");
        ADDER_FCN("f8",double)
        ADDER_FCN("f4",float)
        ADDER_FCN("u8",uint64_t)
        ADDER_FCN("i8",int64_t)
        ADDER_FCN("u4",uint32_t)
        ADDER_FCN("i4",int32_t)
        ADDER_FCN("u2",uint16_t)
        ADDER_FCN("i2",int16_t)
        ADDER_FCN("u1",uint8_t)
        ADDER_FCN("i1",int8_t)
        else throw error("unknown type");
        if (c=='\0') break;
        str += i+1;
        i = 0;
      } else ++i;
    }
  }
  void operator()(void* a, const void* b, size_t len) noexcept {
    const void* last = reinterpret_cast<const char*>(a) + len;
    while (a < last) for (auto& f : fs) f(a,b);
  }
};

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname;
  bool merge_xsec = false,
       normalize = false,
       verbose = false;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input files to merge",req(),pos())
      (ofname,'o',"output file",req())
      (normalize,{"-n","--norm"},"normalize after adding")
      (merge_xsec,{"-x","--xsec","--nlo"},
       "merge cross sections (e.g. NLO parts)")
      (verbose,'v',"verbose")
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  std::vector<std::array<std::string,2>> cols;
  std::vector<std::set<ivanp::sqlite::value>> values;

  auto hists_cmp = [](const auto& a, const auto& b){
    return std::lexicographical_compare(
      a.begin(), a.end(),
      b.begin(), b.end(),
      [](const auto& a, const auto& b){ return (*a) < (*b); }
    );
  };
  std::map<
    std::vector<decltype(values)::value_type::iterator>,
    std::vector<ivanp::sqlite::value>,
    decltype(hists_cmp)
  > hists(hists_cmp);

  decltype(hists)::key_type key;
  decltype(hists)::mapped_type data;
  unsigned ncols=0, ndata=0;
  int i_head=-1, i_data=-1;

  std::map<std::string,long unsigned int> counts;

  std::map<
    std::string,
    std::vector<std::tuple<
      std::array<std::string,2>,
      std::vector<ivanp::sqlite::value>
    >>
  > dicts;

  for (unsigned f=0; f<ifnames.size(); ++f) {
    if (verbose) cout << '\n';
    cout << ifnames[f] << endl;
    ivanp::sqlite db(ifnames[f]);

    // TODO: compare axis definitions

    // read counts ==================================================
    { auto stmt = db.prepare("SELECT * from num");
      while (stmt.step()) {
        const std::string name = stmt.column_text(0);
        auto cnt = stmt.column_int64(1);
        if (merge_xsec && name=="norm") {
          if (cnt!=1) {
            cerr << "\033[31m" "input files must be normalized before merging"
              " cross sections" "\033[0m" << endl;
            return 1;
          }
        } else counts[name] += cnt;
      }
    }

    // read and verify column names =================================
    { auto stmt = db.prepare("PRAGMA table_info(hist)");
      if (!f) {
        while (stmt.step()) {
          cols.push_back({
            stmt.column_text(1),
            stmt.column_text(2)
          });
          const auto& col = cols.back()[0];
          if (col[0]=='_') {
            if (col=="_data") i_data = ndata;
            else if (col=="_head") i_head = ndata;
            ++ndata;
          }
        }
        for (const auto& col : cols)
          cout << col[0] <<' '<< col[1] << endl;
        if (i_data==-1 || i_head==-1) {
          cerr << "\033[31m" "missing column \""
               << (i_data==-1 ? "_data" : "_head")
               << "\" in file \"" << ifnames[f] << "\"" "\033[0m" << endl;
          return 1;
        }
        ncols = cols.size();
        values.resize(ncols);
        key.reserve(ncols-ndata);
        data.reserve(ndata);

      } else {
        for (unsigned i=0; stmt.step(); ++i) {
          const decltype(cols)::value_type col = {
            stmt.column_text(1),
            stmt.column_text(2)
          };
          if (cols[i]!=col) {
            cerr << "\033[31m" "column " << (i+1) << ":" "\033[0m" " got \""
                 << col[0] <<' '<< col[1] << "\", expected \""
                 << cols[i][0] <<' '<< cols[i][1] << "\""
                 << endl;
            return 1;
          }
        }
      }
    }

    // read dictionaries ============================================
    { auto stmt = db.prepare(
        "SELECT name FROM sqlite_master WHERE type='table'");
      while (stmt.step()) {
        const std::string name = stmt.column_text(0);
        if (!(ivanp::ends_with(name,"_dict") &&
            has(cols,[
              col_name=("_"+name.substr(0,name.size()-5))
            ](const auto& col){ return col[0] == col_name; })
        )) continue;
        auto& dict0 = dicts[name];
        decltype(dicts)::mapped_type dict;
        for (auto stmt = db.prepare(cat("PRAGMA table_info(",name,")"));
             stmt.step(); )
        {
          dict.emplace_back();
          auto& col = get<0>(dict.back());
          col = { stmt.column_text(1), stmt.column_text(2) };
        }
        if (!f) {
          dict0 = std::move(dict);
        } else if (dict != dict0) {
          cerr << "\033[31m" "unequal dictionaries \"" << name << "\""
                  "\033[0m" << endl;
          return 1;
        }
      }
    }

    // read rows ====================================================
    { auto stmt = db.prepare("SELECT * from hist");
      while (stmt.step()) {
        key.clear();
        data.clear();

        for (unsigned i=0; i<ncols; ++i) {
          if (cols[i][0][0]=='_') {
            data.emplace_back(stmt.column_value(i));
          } else {
            key.emplace_back(
              values[i].insert(stmt.column_value(i)).first
            );
          }
        }

        if (verbose) {
          for (const auto& x : key)
            cout << x->as_text() << ' ';
          cout << endl;
        }

        auto h = hists.find(key);
        if (h == hists.end()) {
          hists.emplace(std::move(key),std::move(data));
        } else {
          if (!f) {
            cerr << "\033[31m" "repeated histogram in first file"
                    "\033[0m" << endl;
            return 1;
          }

          const auto& d = h->second;
          for (unsigned i=0; i<ndata; ++i) {
            if (i==(unsigned)i_data) continue;
            if (d[i]!=data[i]) {
              cerr << "\033[31m" "unequal data values:" "\033[0m" " got \""
                   << data[i].as_text() << "\", expected \""
                   << d[i].as_text() << "\"" << endl;
              return 1;
            }
          }

          const auto len = d[i_data].bytes();
          if (data[i_data].bytes() != len) {
            cerr << "\033[31m" "data of different length" "\033[0m"
                 << len << '\n'
                 << data[i_data].bytes() << endl;
            return 1;
          }

          adder(d[i_head].as_text())(
            const_cast<void*>(d[i_data].as_blob()),
            data[i_data].as_blob(),
            len
          );
        }
      }
    }
  }

  cout <<'\n'<< hists.size() << " histograms" << endl;

  // normalize ======================================================
  auto& norm = counts["norm"];
  if (merge_xsec) norm = 1;

  if (normalize && norm!=1) { // TODO: generalize
    const double coeff = 1./norm, coeff2 = coeff*coeff;
    for (auto& h : hists) {
      auto& data = h.second[i_data];
      union {
        const void* v;
        char* c;
        double* d;
      } p;
      p.v = data.as_blob();
      const char* end = p.c + data.bytes();
      while (p.c < end) {
        (*(p.d++)) *= coeff;
        (*(p.d++)) *= coeff2;
        p.c += 4;
      }
    }
    norm = 1;
  }

  // write output ===================================================
  cout << "writing " << ofname << endl;
  std::remove(ofname);
  ivanp::sqlite db(ofname);
  db("PRAGMA page_size=4096")("PRAGMA cache_size=8000");

  // write counts ---------------------------------------------------
  db("CREATE TABLE num(value TEXT, n INT)");
  { auto stmt = db.prepare("INSERT INTO num VALUES (?,?)");
    for (const auto& x : counts)
      stmt.bind_row(x.first,x.second);
  }

  // write dictionaries ---------------------------------------------
  for (const auto& dict : dicts) {
    std::stringstream sql;
    sql << "CREATE TABLE " << dict.first << '(';
    bool first = true;
    for (const auto& col : dict.second) {
      if (first) first = false;
      else sql << ", ";
      sql << get<0>(col)[0] <<' '<< get<0>(col)[1];
    }
    sql << ')';
    db(sql.str());
  }

  // write histograms -----------------------------------------------
  { std::stringstream sql;
    sql << "CREATE TABLE hist(";
    for (unsigned i=0; i<ncols; ++i) {
      if (i) sql << ',';
      sql << "\n " << cols[i][0] << ' ' << cols[i][1];
    }
    sql << "\n)";
    db(sql.str());
  }

  { std::stringstream sql;
    sql << "INSERT INTO hist VALUES (";
    for (unsigned i=0; i<ncols; ++i) {
      if (i) sql << ',';
      sql << '?';
    }
    sql << ')';
    auto stmt = db.prepare(sql.str());
    db("BEGIN");
    for (const auto& h : hists) {
      for (unsigned i=0, d=0, k=0; i<ncols; ++i) {
        const bool is_data = (cols[i][0][0] == '_');
        stmt.bind(i+1, is_data ? h.second[d++] : *h.first[k++]);
      }
      stmt.step();
      stmt.reset();
    }
    db("END");
  }
}
