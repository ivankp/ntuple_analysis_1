#ifndef NLOHMANN_JSON_REWEIGHTER_HH
#define NLOHMANN_JSON_REWEIGHTER_HH

#include "json/boost_optional.hh"

namespace nlohmann {

template <>
struct adl_serializer<reweighter::args_struct> {
  static void from_json(const json& j, reweighter::args_struct& args) {
    args.pdf = j.at("pdf");
    args.scale = j.at("scale");
    args.pdf_var = j.value("pdf_var",false);
    for (const reweighter::ren_fac<double>& k : j.at("ren_fac"))
      args.add_scale(k);
  }
};

template <typename T>
struct adl_serializer<reweighter::ren_fac<T>> {
  template <typename U>
  static inline void convert(const json& j, U& x) { x = j.get<U>(); }
  static void from_json(const json& j, reweighter::ren_fac<T>& x) {
    if (j.size()!=2) throw std::runtime_error(
      "reweighter::ren_fac must be represented as an array of size 2");
    convert(j[0],x.ren);
    convert(j[1],x.fac);
  }
};

}

#endif
