#ifndef NLOHMANN_JSON_STD_REGEX_HH
#define NLOHMANN_JSON_STD_REGEX_HH

namespace nlohmann {
template <> struct adl_serializer<std::regex> {
  static void from_json(const json& j, std::regex& re) {
    re.assign(j.get<std::string>());
  }
};
}

#endif
