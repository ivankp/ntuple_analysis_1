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
using ivanp::cat;
using ivanp::y_combinator;

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

  // ivanp::scribe::reader sr(ifnames);
}
