#ifndef NLOHMANN_JSON_BINS_HH
#define NLOHMANN_JSON_BINS_HH

namespace nlohmann {

template <typename T>
struct adl_serializer<multiweight_bin<T>> {
  static void to_json(json& j, const multiweight_bin<T>& bin) {
    j = json(bin.bins);
  }
};

template <typename Bin, typename E, typename... Es>
struct adl_serializer<ivanp::category_bin<Bin,E,Es...>> {
  static void to_json(json& j, const ivanp::category_bin<Bin,E,Es...>& bin) {
    static constexpr auto str = enum_traits<E>::all_str();
    for (unsigned i=0; i<enum_traits<E>::size; ++i) {
      j.emplace(str[i],bin.bins[i]);
    }
  }
};

template <>
struct adl_serializer<nlo_bin> {
  static void to_json(json& j, const nlo_bin& bin) {
    j = {
      { "w" , bin.w  },
      { "w2", bin.w2 },
      { "n" , bin.n  }
    };
  }
};

}

#endif
