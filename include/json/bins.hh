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

// template <typename T>
// struct adl_serializer<multiweight_bin<T>> {
//   static void to_json(json& j, const multiweight_bin<T>& bin) {
//     for (const auto& b : bin.bins) j.emplace_back(b);
//   }
// };
//
// template <typename Bin, typename E, typename... Es>
// struct adl_serializer<ivanp::category_bin<Bin,E,Es...>> {
//   static void to_json(json& j, const ivanp::category_bin<Bin,E,Es...>& bin) {
//     for (const auto& b : bin.bins) j.emplace_back(b);
//   }
// };

template <>
struct adl_serializer<nlo_bin> {
  static void to_json(json& j, const nlo_bin& bin) {
    j = { bin.w, bin.w2, bin.n };
  }
};

}

#endif
