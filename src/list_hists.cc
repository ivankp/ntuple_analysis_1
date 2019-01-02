#include <iostream>

#include "ivanp/io/mem_file.hh"
#include "ivanp/scribe.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;

using ivanp::scribe::size_type;
using namespace ivanp;
using nlohmann::json;

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cout << "usage: " << argv[0] << " file.dat[.xz]\n";
    return 1;
  }

  const mem_file file = (
    ends_with(argv[1],".xz") || ends_with(argv[1],".lzma")
    ? mem_file::pipe(cat("unxz -c ",argv[1]).c_str())
    : mem_file::mmap(argv[1])
  );

  ivanp::scribe::reader sr(file.mem(),file.size());
  cout << "[";
  bool first = true;
  for (const auto& x : sr.get_type()) {
    if (first) first = false;
    else cout << ',';
    cout << '\"' << x.name << '\"';
  }
  cout << "]";
}
