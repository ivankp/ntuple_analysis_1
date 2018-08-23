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
#include "ivanp/program_options.hh"
#include "glob.hh"
#include "copy_file.hh"

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

#include <cstdio>
#include <boost/core/noncopyable.hpp>
class tmp_file_wrap: private boost::noncopyable {
  string fname;
public:
  tmp_file_wrap(tmp_file_wrap&& x): fname(std::move(x.fname)) { }
  template <typename... T>
  tmp_file_wrap(T&&... x): fname(std::forward<T>(x)...) { }
  ~tmp_file_wrap() {
    if (!fname.empty())
      if (remove(fname.c_str()))
        perror(("\033[31mError deleting file\033[31m \""+fname+"\"").c_str());
  }
};

inline nlohmann::json::json_pointer operator "" _jp(const char* s, size_t n) {
  return nlohmann::json::json_pointer(std::string(s,n));
}

#ifndef ANALYSIS
#error "ANALYSIS is not defined"
#endif

#define ANALYSIS_GLOBAL
#include STR(ANALYSIS)
#undef ANALYSIS_GLOBAL

int main(int argc, char* argv[]) {
  std::vector<const char*> card_names;
  std::string tmp_dir;

  try {
    using namespace ivanp::po;
    if (program_options()
      (card_names,'i',"json runcards",req(),pos())
      (tmp_dir,"--tmp-dir","copy input ntuples to temporary directory")
      .parse(argc,argv,true)) return 0;
    if (!tmp_dir.empty() && tmp_dir.back()!='/') tmp_dir += '/';
  } catch (const std::exception& e) {
    cerr << "\033[31m" << e.what() << "\033[0m" << endl;
    return 1;
  }
  // ================================================================

  // Read runcards ==================================================
  nlohmann::json runcards;
  std::ifstream(card_names[0]) >> runcards;
  for (const auto& card_name : card_names) {
    nlohmann::json runcard;
    std::ifstream(card_name) >> runcard;
    runcards.merge_patch(runcard);
  }

  // Chain and friend input files ===================================
  std::list<tmp_file_wrap> tmp_files;
  std::list<TChain> chains;
  vector<string> weights_names;
  for (const auto& input : runcards.at("input")) {
    const string info = input.value("info","");
    if (!info.empty()) cout << "\033[36mInfo\033[0m: " << info << endl;

    // const bool getnentries = input.value("getnentries",true);

    TChain* chain = nullptr;
    string tree_name = input.value("tree","");
    for (const string& file_glob : input.at("files")) { // Chain input files
      for (string& file_name : ivanp::glob(file_glob)) {
        if (file_name.empty() || file_name.back()=='/') {
          cerr << "\033[31mBad file_name\033[0m: \""
               << file_name << "\"" << endl;
          return 1;
        }
        if (!tmp_dir.empty()) { // copy input files to tmp_dir
          string tmp_file_name = tmp_dir + string(
            std::find(file_name.rbegin(), file_name.rend(), '/').base(),
            file_name.end());

          cout << "\033[36mCopying\033[0m: "
               << file_name << " -> " << tmp_file_name;
          using namespace std::chrono;
          using clock = high_resolution_clock;
          clock::time_point t1 = clock::now();
          copy_file(file_name.c_str(), tmp_file_name.c_str());
          // system(cat("cp -v ",file_name,' ',tmp_file_name).c_str());
          clock::time_point t2 = clock::now();
          cout << cat(std::setprecision(2),std::fixed,
              " (", duration_cast<duration<float>>(t2-t1).count(), " sec)"
            ) << endl;

          tmp_files.emplace_back(tmp_file_name);
          file_name = std::move(tmp_file_name);
        }
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
        if (!chain->Add(file_name.c_str(),0)) return 1;
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
  const std::array<Long64_t,2> entry_range =
    runcards.value<decltype(entry_range)>("entry_range",{0,0});
  reader.SetEntriesRange(entry_range[0],entry_range[1]);
  for ( ivanp::timed_counter<Long64_t> ent( entry_range[0],
          entry_range[1] > entry_range[0]
          ? entry_range[1]
          : reader.GetEntries(true)
        ); reader.Next(); ++ent)
  {

    // TODO: add reweighting

#define ANALYSIS_LOOP
#include STR(ANALYSIS)
#undef ANALYSIS_LOOP

  }

#define ANALYSIS_END
#include STR(ANALYSIS)
#undef ANALYSIS_END
}
