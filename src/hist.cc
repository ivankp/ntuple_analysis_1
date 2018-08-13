#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <regex>

// #include <sqlite3.h>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "nlohmann_json.hpp"

#include "ivanp/string.hh"
#include "ivanp/error.hh"
#include "ivanp/root/tkey.hh"

#include "float_or_double_reader.hh"
#include "glob.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::get;

using ivanp::cat;

namespace nlohmann {
template <> struct adl_serializer<std::regex> {
  static void from_json(const json& j, std::regex& re) {
    re.assign(j.get<std::string>());
  }
};
}

template <typename Key>
std::string opt_str(const nlohmann::json& json, const Key& key){
  auto it = json.find(key);
  if (it!=json.end()) return *it;
  return { };
}

int main(int argc, char* argv[]) {
  if (argc<3 || std::any_of(argv+1,argv+argc,[](const char* arg){
    return !strcmp(arg,"-h") || !strcmp(arg,"--help");
  })) {
    cout << "usage: " << argv[0] << " analysis.so runcard.json ..." << endl;
    return 1;
  }

  // Read runcards ==================================================
  nlohmann::json runcards;
  std::ifstream(argv[2]) >> runcards;
  for (int i=3; i<argc; ++i) {
    nlohmann::json runcard;
    std::ifstream(argv[i]) >> runcard;
    runcards.merge_patch(runcard);
  }

  // Chain and friend input files ===================================
  std::list<TChain> chains;
  std::vector<std::string> weights_names;
  for (const auto& input : runcards["input"]) {
    const auto info = opt_str(input,"info");
    if (!info.empty()) // print info
      cout << "\033[36mInfo\033[0m: " << info << endl;

    const bool getnentries = input["getnentries"];

    TChain* chain = nullptr;
    auto tree_name = opt_str(input,"tree");
    for (const std::string& file_glob : input["files"]) { // Chain input files
      for (const auto& file_name : ivanp::glob(file_glob)) {
        if (!chain) {
          if (tree_name.empty()) {
            std::vector<std::string> trees;
            TFile file(file_name.c_str());
            for (const TKey& key : list_cast<TKey>(file.GetListOfKeys())) {
              const auto* key_class = TClass::GetClass(key.GetClassName(),true);
              if (!key_class) continue;
              if (key_class->InheritsFrom(TTree::Class()))
                trees.emplace_back(key.GetName());
            }
            if (trees.size()>1) throw ivanp::error(
              "Multiple TTrees in file ",file_name,": ",ivanp::lcat(trees));
            tree_name = trees[0];
          }
          cout << "\033[36mTTree\033[0m: " << tree_name << '\n';
          cout << "\033[36mChaining input files\033[0m:" << endl;
          chains.emplace_back(tree_name.c_str());
          chain = &chains.back();
        }
        if (!chain->Add(file_name.c_str(),getnentries?0:TTree::kMaxEntries))
          return 1;
        cout << "  " << file_name << endl;
      }
    }

    auto weights = input.find("weights"); // Add weights
    if (weights!=input.end()) {
      std::vector<std::regex> res = *weights;
      for (const auto* b : *chain->GetListOfBranches()) {
        const char* bname = b->GetName();
        for (const auto& re : res)
          if (std::regex_match(bname,re)) {
            weights_names.emplace_back(bname);
            break;
          }
      }
    }

    if (chains.size()>1) // Friend the chains
      chains.front().AddFriend(chain);
  }
  TChain& chain = chains.front();

  if (weights_names.empty()) { // Find default weights
    const char* def[] = { "weight2", "weight" };
    for (unsigned i=0, n=sizeof(def)/sizeof(*def); i<n; ++i)
      for (const auto* b : *chain.GetListOfBranches())
        if (!strcmp(b->GetName(),def[i])) {
          weights_names.emplace_back(def[i]);
          goto endloop;
        }
    cout << "\033[31mCannot find weight branches\033[0m" << endl;
    return 1;
    endloop: ;
  }

  TTreeReader reader(&chain);

  std::vector<float_or_double_value_reader> _weights;
  _weights.reserve(weights_names.size());
  cout << "\033[36mWeights\033[0m:\n";
  for (const auto& name : weights_names) { // Make weight readers
    cout << "  " << name << endl;
    _weights.emplace_back(reader,name.c_str());
  }

  // Analysis =======================================================
  // analysis ana(chain);
}
