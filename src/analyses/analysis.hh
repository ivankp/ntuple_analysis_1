#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <regex>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TChain.h>
#include "branch_reader.hh"
#include "reweighter.hh"

#include "json/nlohmann.hpp"
#include "json/print_value_t.hh"
#include "json/std_regex.hh"
#include "json/reweighter.hh"

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
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  // Read runcards ==================================================
  nlohmann::json runcards;
  for (const auto& filename : card_names) {
    nlohmann::json runcard;
    std::ifstream file(filename);
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      file >> runcard;
    } catch (const std::exception& e) {
      cerr << "\033[31mError reading file\033[0m \"" << filename << "\": "
           << e.what() << endl;
      return 1;
    }
    runcards.merge_patch(runcard);
  }

  // Chain and friend input files ===================================
  vector<string> tmp_files;
  std::list<TChain> chains;
  vector<string> weights_names;
  for (auto& input : runcards.at("input")) {
    const string info = input.value("info","");
    if (!info.empty()) cout << "\033[36mInfo\033[0m: " << info << endl;

    TChain* chain = nullptr;
    string tree_name = input.value("tree","");
    vector<string> file_names;
    for (const string& file_glob : input.at("files")) { // Chain input files
      for (string& file_name : ivanp::glob(file_glob)) {
        if (file_name.empty() || file_name.back()=='/') {
          cerr << "\033[31mBad file_name\033[0m: \""
               << file_name << "\"" << endl;
          return 1;
        }
        file_names.push_back(file_name);
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
    input["files"] = file_names;

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

  // delete the blocking file
  // system("ls -lh /home/ivanp/.condor_cp_mutex");
  // system("rm -fv /home/ivanp/.condor_cp_mutex");
  // if (condor_cp_mutex) condor_cp_mutex = boost::none;

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

  // Make weight readers
  vector<float_reader> _weights;
  _weights.reserve(weights_names.size());
  cout << "\033[36mWeights\033[0m:\n";
  for (const auto& name : weights_names) {
    cout << "  " << name << endl;
    _weights.emplace_back(reader,name.c_str());
  }

  // Make reweighters
  vector<reweighter> reweighters;
  const auto reweighting = runcards.find("reweighting");
  if (reweighting!=runcards.end()) {
    reweighters.reserve(reweighting->size());
    for (const auto& j : *reweighting)
      reweighters.emplace_back(reader,j);
  }

  // Reserve weights vector
  vector<double> weights(std::accumulate(
    reweighters.begin(), reweighters.end(), _weights.size(),
    [](unsigned n, auto& rew){ return n + rew.nweights(); }
  ));

  for (auto& rew : reweighters) // get weights names
    for (unsigned i=0; i<rew.nweights(); ++i)
      weights_names.emplace_back(rew.weight_name(i));

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
    { unsigned wi = 0;
      for (; wi<_weights.size(); ++wi)
        weights[wi] = *_weights[wi];
      for (auto& rew : reweighters) {
        rew(); // reweight this event
        for (unsigned i=0; i<rew.nweights(); ++i, ++wi)
          weights[i] = rew[i];
      }
    }

    // TEST(*_weights[0])
    // TEST(*_w2)
    // TEST((*_weights[0]==*_w2))
    //
    // if (ent < 10) continue; else break;

#define ANALYSIS_LOOP
#include STR(ANALYSIS)
#undef ANALYSIS_LOOP

  }

#define ANALYSIS_END
#include STR(ANALYSIS)
#undef ANALYSIS_END

  for (const string& fname : tmp_files) { // delete temporary files
    if (!fname.empty())
      if (remove(fname.c_str()))
        perror(("\033[31mError deleting file\033[0m \""+fname+"\"").c_str());
  }
}
