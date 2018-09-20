#ifndef BIN_DEFS_HH
#define BIN_DEFS_HH

#include <vector>
#include "ivanp/debug/at.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

struct lo_bin {
  double w = 0, w2 = 0;
  long unsigned n = 0;
  void operator()() noexcept {
    ++w; ++w2; ++n;
  }
  void operator()(double weight) noexcept {
    w  += weight;
    w2 += weight*weight;
    ++n;
  }
  lo_bin& operator+=(const lo_bin& b) noexcept {
    w  += b.w;
    w2 += b.w2;
    n  += b.n;
    return *this;
  }
  double err2() const noexcept { return w2; }
};

struct nlo_bin_base {
  int id = 0; // event id
  static int current_id;
};
int nlo_bin_base::current_id;

template <typename T, typename SFINAE=void> struct nlo_bin;

template <typename T>
struct nlo_bin<T,std::enable_if_t<std::is_arithmetic<T>::value>>
: nlo_bin_base {
  using type = T;
  T w=0, wtmp=0, w2=0;
  long unsigned n = 0;
  void operator()() noexcept {
    if (id == current_id) ++wtmp;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      ++wtmp;
    }
    ++w;
    ++n;
  }
  void operator()(const T& weight) noexcept {
    if (id == current_id) wtmp += weight;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      wtmp = weight;
    }
    w += weight;
    ++n;
  }
  nlo_bin& operator+=(const nlo_bin& b) noexcept {
    wtmp += b.wtmp;
    w  += b.w;
    w2 += b.w2;
    n  += b.n;
    return *this;
  }
  void finalize() noexcept { w2 += wtmp*wtmp; }
};

template <typename T>
struct nlo_bin<T[],std::enable_if_t<std::is_arithmetic<T>::value>>
: nlo_bin_base {
  using type = T;
  static std::vector<T> weights;
  static unsigned wi;
  struct w_struct { T w=0, wtmp=0, w2=0; };
  std::vector<w_struct> ws;
  long unsigned n = 0;
  nlo_bin(): nlo_bin_base(), ws(weights.size()) { }
  void operator()() noexcept {
    if (id == current_id) {
      for (unsigned i=ws.size(); i; ) {
        --i;
        auto  w = weights AT(i);
        auto& _ = ws AT(i);
        _.w    += w;
        _.wtmp += w;
      }
    } else {
      id = current_id;
      for (unsigned i=ws.size(); i; ) {
        --i;
        auto  w = weights AT(i);
        auto& _ = ws AT(i);
        _.w   += w;
        _.w2  += _.wtmp*_.wtmp;
        _.wtmp = w;
      }
    }
    ++n;
  }
  nlo_bin& operator+=(const nlo_bin& b) noexcept {
    for (unsigned i=ws.size(); i; ) {
      --i;
      auto& _b = b.ws AT(i);
      auto& _  =   ws AT(i);
      _.w    += _b.w;
      _.wtmp += _b.wtmp;
      _.w2   += _b.w2;
    }
    n += b.n;
    return *this;
  }
  void finalize() noexcept {
    for (unsigned i=ws.size(); i; ) {
      --i;
      auto& _ = ws AT(i);
      _.w2 += _.wtmp*_.wtmp;
    }
  }
};
#define NLO_BIN_STATIC(TYPE,NAME) \
  template <typename T> \
  TYPE nlo_bin<T[],std::enable_if_t<std::is_arithmetic<T>::value>>::NAME

NLO_BIN_STATIC(std::vector<T>,weights);
NLO_BIN_STATIC(unsigned,wi);

struct profile_bin {
  double w = 0, // total weight
         m = 0, // mean
         s = 0; // variance * (n-1)
  long unsigned n = 0;
  void operator()(double x) noexcept {
    ++w;
    if (n == 0) {
      m = x;
    } else {
      const double d = x-m;
      m += d/w;
      s += d*(x-m);
    }
    ++n;
  }
  void operator()(double weight, double x) noexcept {
    w += weight;
    if (n == 0) {
      m = x;
    } else {
      const double wd = weight*(x-m);
      m += wd/w;
      s += wd*(x-m);
    }
    ++n;
  }
  double var() const noexcept { return (n > 1) ? n*s/((n-1)*w) : 0.; }
};

struct multiweight_bin_base {
  static std::vector<double> weights;
  static unsigned wi;
};
std::vector<double> multiweight_bin_base::weights;
unsigned multiweight_bin_base::wi;

template <typename Bin>
struct multiweight_bin: multiweight_bin_base {
  std::vector<Bin> bins;
  multiweight_bin(): bins(weights.size()) { }

  template <typename... T>
  multiweight_bin& operator()(T... x) noexcept {
    for (unsigned i=weights.size(); i; ) {
      --i; bins AT(i)(weights AT(i),x...);
    }
    return *this;
  }
  multiweight_bin& operator+=(const multiweight_bin& rhs) noexcept {
    for (unsigned i=weights.size(); i; ) {
      --i; bins AT(i) += rhs.bins AT(i);
    }
    return *this;
  }
  Bin& operator->() noexcept { return bins AT(wi); }
  const Bin& operator->() const noexcept { return bins AT(wi); }
};

#endif
