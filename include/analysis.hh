#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <cstring>
#include <vector>

#include <boost/optional.hpp>

#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "nlohmann/json.hpp"
#include "nlohmann/print_value_t.hh"
#include "nlohmann/std_regex.hh"

#include "float_or_double_reader.hh"

struct analysis_args {
  TTreeReader& reader;
  const std::vector<double>& weights;
  const nlohmann::json& conf;
};

class analysis_base: public analysis_args {
protected:
  TTreeReaderValue<Int_t> _id, _nparticle;
  TTreeReaderArray<Int_t> _kf;
  float_or_double_array_reader _px, _py, _pz, _E;
  boost::optional<TTreeReaderValue<Int_t>> _ncount;

public:
  analysis_base(analysis_args&& args)
  : analysis_args(std::move(args)),
    _id(reader,"id"), _nparticle(reader,"nparticle"), _kf(reader,"kf"),
    _px(reader,"px"), _py(reader,"py"), _pz(reader,"pz"), _E(reader,"E" )
  {
    for (auto* b : *reader.GetTree()->GetListOfBranches()) {
      if (!strcmp(b->GetName(),"ncount")) {
        _ncount.emplace(reader,"ncount"); break;
      }
    }
  }

  virtual void event_loop() = 0;
  virtual void write_output() = 0;
};

#define ANALYSIS_FACTORY(CLASS_NAME) \
  extern "C" analysis_base* make_analysis(analysis_args&& args) { \
    return new CLASS_NAME(std::move(args)); \
  }

#endif

