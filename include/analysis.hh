#ifndef ANALYSIS_HH
#define ANALYSIS_HH

#include <cstring>

#include <boost/optional.hpp>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "nlohmann/json.hpp"
#include "nlohmann/print_value_t.hh"
#include "nlohmann/std_regex.hh"

#include "float_or_double_reader.hh"

class basic_analysis {
  TTreeReaderValue<Int_t> _id, _nparticle;
  TTreeReaderArray<Int_t> _kf;
  float_or_double_array_reader _px, _py, _pz, _E;
  boost::optional<TTreeReaderValue<Int_t>> _ncount;
public:
  basic_analysis(TTreeReader& r)
  : _id(r,"id"), _nparticle(r,"nparticle"), _kf(r,"kf"),
    _px(r,"px"), _py(r,"py"), _pz(r,"pz"), _E(r,"E" )
  {
    for (auto* b : *r.GetTree()->GetListOfBranches()) {
      if (!strcmp(b->GetName(),"ncount")) {
        _ncount.emplace(r,"ncount"); break;
      }
    }
  }
};

#endif
