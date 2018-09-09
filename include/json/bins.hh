#ifndef NLOHMANN_JSON_BINS_HH
#define NLOHMANN_JSON_BINS_HH

#include "ivanp/detect.hh"

template <typename T>
using member_bins_t = decltype(std::declval<const T&>().bins);

namespace nlohmann {

template <typename T>
struct adl_serializer<
  T, std::enable_if_t<ivanp::is_detected<member_bins_t,T>::value>
> {
  static void to_json(json& j, const T& bin) {
    for (const auto& b : bin.bins) j.emplace_back(b);
  }
};

template <typename T>
struct adl_serializer<
  nlo_bin<T>, std::enable_if_t<!std::is_array<T>::value>
> {
  static void to_json(json& j, const nlo_bin<T>& bin) {
    bin.finalize();
    j = { bin.w, bin.w2 + bin.w*bin.w, bin.n };
  }
};

template <typename T>
struct adl_serializer<
  nlo_bin<T>, std::enable_if_t<std::is_array<T>::value>
> {
  static void to_json(json& j, const nlo_bin<T>& bin) {
    auto j0 = json::array();
    for (const auto& w : bin.ws) j0.push_back({ w.w, w.w2 + w.w*w.w });
    j = { j0, bin.n };
  }
};

}

#endif
