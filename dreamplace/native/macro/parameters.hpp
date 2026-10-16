// BSD-3-Clause. See LICENSE.
#pragma once
#include <span>
#include <string>
#include <string_view>
namespace DreamPlace::native::macro {
struct Parameters {
#define H2D_INT(n, d, h) int n = d;
#define H2D_REAL(n, d, h) double n = d;
#define H2D_BOOL(n, d, h) bool n = d;
#define H2D_STRING(n, d, h) std::string n = d;
#include "dreamplace/native/macro/params.def"
#undef H2D_INT
#undef H2D_REAL
#undef H2D_BOOL
#undef H2D_STRING
  void validate() const;
};
enum class ParameterKind { integer, real, boolean, string };
struct ParameterSpec {
  std::string_view name, default_value, description;
  ParameterKind kind;
};
std::span<const ParameterSpec> parameter_specs();
void set_parameter(Parameters &parameters, std::string_view name,
                   std::string_view value);
std::string parameter_value(const Parameters &parameters,
                            std::string_view name);
std::string parameter_key(const Parameters &parameters);
} // namespace DreamPlace::native::macro
