#ifndef NLOHMANN_JSON_BINNER_HH
#define NLOHMANN_JSON_BINNER_HH

namespace nlohmann {

template <typename Bin, typename... Ax, typename Container, typename Filler>
struct adl_serializer<ivanp::binner<Bin,std::tuple<Ax...>,Container,Filler>> {
private:
  using hist = ivanp::binner<Bin,std::tuple<Ax...>,Container,Filler>;
  using i_t = typename hist::size_type;
  using ii_t = typename hist::index_array_type;

  template <unsigned I = sizeof...(Ax)-1>
  static inline std::enable_if_t<I!=0>
  write(json& j, ii_t& ii, const hist& h) {
    using spec = typename hist::template axis_spec<I>;
    const i_t n = std::get<I>(h.axes()).nbins() + spec::under::value;
    for (i_t& i = std::get<I>(ii) = 0; i<=n; ++i) {
      if (i==0 && !spec::under::value) j[i] = nullptr; else
      if (i==n && !spec::over ::value) j[i] = nullptr; else
      write<I>( j[i], ii, h );
    }
  }
  template <unsigned I = sizeof...(Ax)-1>
  static inline std::enable_if_t<I==0>
  write(json& j, ii_t& ii, const hist& h) {
    using spec = typename hist::template axis_spec<I>;
    const i_t n = std::get<I>(h.axes()).nbins() + spec::under::value;
    for (i_t& i = std::get<I>(ii) = 0; i<=n; ++i) {
      if (i==0 && !spec::under::value) j[i] = nullptr; else
      if (i==n && !spec::over ::value) j[i] = nullptr; else
      j[i] = h[ii];
    }
  }
public:
  static void to_json(json& j, const hist& h) {
    j["axes"] = h.axes();
    ii_t ii { };
    write(j["bins"],ii,h);
  }
};

template <typename Container, bool Inherit>
struct adl_serializer<ivanp::container_axis<Container,Inherit>> {
  static void to_json(json& j, const ivanp::container_axis<Container,Inherit>& axis) {
    j = axis.edges();
  }
};

template <typename T, bool Inherit>
struct adl_serializer<ivanp::uniform_axis<T,Inherit>> {
  static void to_json(json& j, const ivanp::uniform_axis<T,Inherit>& axis) {
    j["nbins"] = axis.nbins();
    j["range"] = { axis.min(), axis.max() };
  }
};

template <typename T, typename Ref, bool Inherit>
struct adl_serializer<ivanp::ref_axis<T,Ref,Inherit>> {
  static void to_json(json& j, const ivanp::ref_axis<T,Ref,Inherit>& axis) {
    if (axis.is_uniform()) {
      j["nbins"] = axis.nbins();
      j["range"] = { axis.min(), axis.max() };
    } else {
      j = vector_of_edges(axis);
    }
  }
};

template <typename T>
struct adl_serializer<ivanp::abstract_axis<T>> {
  static void to_json(json& j, const ivanp::abstract_axis<T>& axis) {
    if (axis.is_uniform()) {
      j["nbins"] = axis.nbins();
      j["range"] = { axis.min(), axis.max() };
    } else {
      j = vector_of_edges(axis);
    }
  }
};

}

#endif
