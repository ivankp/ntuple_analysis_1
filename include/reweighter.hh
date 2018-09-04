#ifndef NTUPLES_REWEIGHTER_HH
#define NTUPLES_REWEIGHTER_HH

#include <vector>
#include <string>
#include <boost/optional.hpp>
#include <TTreeReader.h>

class reweighter_impl;

class reweighter {
  reweighter_impl * const impl;

public:
  struct args_struct {
    std::string pdf, scale;
    bool pdf_var;
    template <typename T>
    struct ren_fac { boost::optional<T> ren, fac; };
    std::vector<double> K;
    std::vector<ren_fac<unsigned>> Ki;
    void add_scale(const ren_fac<double>& k);
  };

  reweighter(TTreeReader& reader, args_struct args);
  ~reweighter();

  unsigned nweights() const;
  double operator[](unsigned i) const;
  std::string name(unsigned i) const;
};

#endif
