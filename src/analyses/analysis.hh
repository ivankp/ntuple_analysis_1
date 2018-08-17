#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <regex>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "float_or_double_reader.hh"

#include "json/nlohmann.hpp"
#include "json/print_value_t.hh"
#include "json/std_regex.hh"

#include "ivanp/string.hh"
#include "ivanp/root/tkey.hh"
#include "ivanp/error.hh"
#include "ivanp/timed_counter.hh"
#include "glob.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using ivanp::cat;

template <typename T, typename Key>
T opt(const nlohmann::json& json, const Key& key, T def={}){
  auto it = json.find(key);
  if (it!=json.end()) return *it;
  return def;
}

#ifndef ANALYSIS
#error "ANALYSIS is not defined"
#endif

#define ANALYSIS_GLOBAL
#include STR(ANALYSIS)
#undef ANALYSIS_GLOBAL

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
  vector<string> weights_names;
  for (const auto& input : runcards["input"]) {
    const string info = opt<string>(input,"info");
    if (!info.empty()) cout << "\033[36mInfo\033[0m: " << info << endl;

    const bool getnentries = opt<bool>(input,"getnentries",true);

    TChain* chain = nullptr;
    string tree_name = opt<string>(input,"tree");
    for (const string& file_glob : input["files"]) { // Chain input files
      for (const string& file_name : ivanp::glob(file_glob)) {
        if (!chain) {
          if (tree_name.empty()) {
            vector<string> trees;
            TFile file(file_name.c_str());
            for (const TKey& key : list_cast<TKey>(file.GetListOfKeys())) {
              const auto* key_class = TClass::GetClass(key.GetClassName(),true);
              if (!key_class) continue;
              if (key_class->InheritsFrom(TTree::Class()))
                trees.emplace_back(key.GetName());
            }
            if (trees.size()>1) {
              cerr << "\033[31mMultiple TTrees in file " << file_name
                   << ": " << ivanp::lcat(trees) << "\033[0m" << endl;
              return 1;
            }
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
      vector<std::regex> res;
      if (weights->is_array()) res = weights->get<decltype(res)>();
      else if (weights->is_string()) res.emplace_back(weights->get<string>());
      else {
        cerr << "\033[31mWrong \"weights\" value type ("
             << weights->type() << ")\033[0m" << endl;
        return 1;
      }
      for (const auto* b : *chain->GetListOfBranches()) {
        const char* bname = b->GetName();
        for (const auto& re : res)
          if (std::regex_match(bname,re)) {
            weights_names.emplace_back(bname);
            break;
          }
      }
    }

    if (chains.size() > 1) // Friend the chains
      chains.front().AddFriend(chain);
  }
  TChain& chain = chains.front();

  if (weights_names.empty()) { // Find default weights
    const char* names[] = { "weight2", "weight" };
    for (unsigned i=0, n=sizeof(names)/sizeof(*names); i<n; ++i)
      for (const auto* b : *chain.GetListOfBranches())
        if (!strcmp(b->GetName(),names[i])) {
          weights_names.emplace_back(names[i]);
          goto endloop;
        }
    cerr << "\033[31mCannot find weight branches\033[0m" << endl;
    return 1;
    endloop: ;
  }

  TTreeReader reader(&chain);

  vector<float_or_double_value_reader> _weights;
  _weights.reserve(weights_names.size());
  cout << "\033[36mWeights\033[0m:\n";
  for (const auto& name : weights_names) { // Make weight readers
    cout << "  " << name << endl;
    _weights.emplace_back(reader,name.c_str());
  }

#define ANALYSIS_INIT
#include STR(ANALYSIS)
#undef ANALYSIS_INIT

  // LOOP ===========================================================
  using counter = ivanp::timed_counter<Long64_t>;
  // TODO: fix counter
  for (counter ent(reader.GetEntries(true)); reader.Next(); ++ent) {

    // TODO: add reweighting

#define ANALYSIS_LOOP
#include STR(ANALYSIS)
#undef ANALYSIS_LOOP

  }

#define ANALYSIS_END
#include STR(ANALYSIS)
#undef ANALYSIS_END
}