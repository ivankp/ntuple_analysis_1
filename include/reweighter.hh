#ifndef NTUPLES_REWEIGHTER_HH
#define NTUPLES_REWEIGHTER_HH

#include <vector>
#include <string>
#include <boost/optional.hpp>
#include <TTreeReader.h>

class reweighter_impl;


class reweighter {
  reweighter_impl *impl;

public:
  template <typename T>
  struct ren_fac { boost::optional<T> ren, fac; };

  struct args_struct {
    std::string pdf, scale;
    bool pdf_var;
    std::vector<double> Kr, Kf;
    std::vector<ren_fac<unsigned>> Ki;
    void add_scale(const ren_fac<double>& k);
  };

  reweighter(TTreeReader& reader, args_struct args);
  reweighter() = delete;
  reweighter(const reweighter&) = delete;
  reweighter(reweighter&&);
  reweighter& operator=(const reweighter&) = delete;
  reweighter& operator=(reweighter&&) = delete;
  ~reweighter();

  void operator()();
  unsigned nweights() const;
  double operator[](unsigned i) const;
  std::string weight_name(unsigned i) const;
};

#endif
