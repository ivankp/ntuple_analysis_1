#ifndef NLOHMANN_JSON_JET_ALGORITHM_HH
#define NLOHMANN_JSON_JET_ALGORITHM_HH

#include "ivanp/string.hh"

namespace nlohmann {

template <>
struct adl_serializer<fastjet::JetAlgorithm> {
  static void from_json(const json& j, fastjet::JetAlgorithm& alg) {
    const std::string str = j;
    const char* arg = str.c_str();
    if (!ivanp::strcmpi(arg,"kt"))
      alg = fastjet::kt_algorithm;
    else if (!ivanp::strcmpi(arg,"antikt") || !ivanp::strcmpi(arg,"akt"))
      alg = fastjet::antikt_algorithm;
    else if (!ivanp::strcmpi(arg,"cambridge") || !ivanp::strcmpi(arg,"ca"))
      alg = fastjet::cambridge_algorithm;
    else throw std::runtime_error(cat(
      "cannot interpret \"",arg,"\" as a FastJet algorithm"));
  }
};

template <>
struct adl_serializer<fastjet::JetDefinition> {
  static void from_json(const json& j, fastjet::JetDefinition& def) {
    def = { j[0].get<fastjet::JetAlgorithm>(), j[1].get<double>() };
  }
};

}

#endif
