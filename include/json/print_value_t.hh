#ifndef NLOHMANN_JSON_PRINT_VALUE_T_HH
#define NLOHMANN_JSON_PRINT_VALUE_T_HH

std::ostream& operator<<(std::ostream& s, nlohmann::json::value_t x) {
  switch (x) {
    case nlohmann::json::value_t::null: return s << "null";
    case nlohmann::json::value_t::object: return s << "object";
    case nlohmann::json::value_t::array: return s << "array";
    case nlohmann::json::value_t::string: return s << "string";
    case nlohmann::json::value_t::boolean: return s << "boolean";
    case nlohmann::json::value_t::number_integer: return s << "number_integer";
    case nlohmann::json::value_t::number_unsigned: return s << "number_unsigned";
    case nlohmann::json::value_t::number_float: return s << "number_float";
    case nlohmann::json::value_t::discarded: return s << "discarded";
    default: return s;
  }
}

#endif

