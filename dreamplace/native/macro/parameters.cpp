// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/parameters.hpp"
#include <charconv>
#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <type_traits>
namespace DreamPlace::native::macro {
namespace {
constexpr ParameterSpec specs[] = {
#define H2D_INT(n, d, h) {#n, #d, h, ParameterKind::integer},
#define H2D_REAL(n, d, h) {#n, #d, h, ParameterKind::real},
#define H2D_BOOL(n, d, h) {#n, #d, h, ParameterKind::boolean},
#define H2D_STRING(n, d, h) {#n, d, h, ParameterKind::string},
#include "dreamplace/native/macro/params.def"
#undef H2D_INT
#undef H2D_REAL
#undef H2D_BOOL
#undef H2D_STRING
};
template <typename T> T number(std::string_view v) {
  T out{};
  if constexpr (std::is_floating_point_v<T>) {
    std::istringstream input{std::string(v)};
    input.imbue(std::locale::classic());
    input >> std::noskipws >> out;
    if (!input || input.peek() != std::char_traits<char>::eof())
      throw std::invalid_argument("invalid H2D parameter number");
  } else {
    auto [end, ec] = std::from_chars(v.data(), v.data() + v.size(), out);
    if (ec != std::errc{} || end != v.data() + v.size())
      throw std::invalid_argument("invalid H2D parameter number");
  }
  return out;
}
} // namespace
std::span<const ParameterSpec> parameter_specs() { return specs; }
void set_parameter(Parameters &p, std::string_view name,
                   std::string_view value) {
#define H2D_INT(n, d, h)                                                       \
  if (name == #n) {                                                            \
    p.n = number<int>(value);                                                  \
    return;                                                                    \
  }
#define H2D_REAL(n, d, h)                                                      \
  if (name == #n) {                                                            \
    p.n = number<double>(value);                                               \
    return;                                                                    \
  }
#define H2D_BOOL(n, d, h)                                                      \
  if (name == #n) {                                                            \
    if (value != "true" && value != "false")                                   \
      throw std::invalid_argument("expected H2D boolean");                     \
    p.n = value == "true";                                                     \
    return;                                                                    \
  }
#define H2D_STRING(n, d, h)                                                    \
  if (name == #n) {                                                            \
    p.n = value;                                                               \
    return;                                                                    \
  }
#include "dreamplace/native/macro/params.def"
#undef H2D_INT
#undef H2D_REAL
#undef H2D_BOOL
#undef H2D_STRING
  throw std::invalid_argument("unknown H2D parameter: " + std::string(name));
}
std::string parameter_value(const Parameters &p, std::string_view name) {
  std::ostringstream out;
  out << std::setprecision(17) << std::boolalpha;
#define H2D_INT(n, d, h)                                                       \
  if (name == #n) {                                                            \
    out << p.n;                                                                \
    return out.str();                                                          \
  }
#define H2D_REAL(n, d, h) H2D_INT(n, d, h)
#define H2D_BOOL(n, d, h) H2D_INT(n, d, h)
#define H2D_STRING(n, d, h) H2D_INT(n, d, h)
#include "dreamplace/native/macro/params.def"
#undef H2D_INT
#undef H2D_REAL
#undef H2D_BOOL
#undef H2D_STRING
  throw std::invalid_argument("unknown H2D parameter");
}
std::string parameter_key(const Parameters &p) {
  std::string key;
  for (const auto &s : specs)
    key += std::string(s.name) + "=" + parameter_value(p, s.name) + ";";
  return key;
}
void Parameters::validate() const {
  auto require = [](bool ok) {
    if (!ok)
      throw std::invalid_argument("invalid H2D parameters");
  };
  for (const auto &s : specs)
    if (s.kind == ParameterKind::real)
      require(std::isfinite(number<double>(parameter_value(*this, s.name))));
  require(K >= 1 && K <= 32 && seed >= 0 && iters_per_outer >= 0 &&
          iters_per_outer <= 100000 && bins >= 2 && bins <= 1024 &&
          threads >= 1 && threads <= 128);
  require(freeze_after >= 1 && freeze_after <= 32 && hops_max >= 0 &&
          hops_max <= 16 && fanout_cap >= 1);
  require(chains >= 1 && chains <= 128 && evals >= 0 && evals <= 100000 &&
          seed_sa_evals >= 0 && seed_sa_evals <= 100000);
  require(lambda0 >= 0 && gamma > 0 && gamma <= 2 && delta >= 0 && tau_q >= 0 &&
          tau_q <= 1 && theta_q >= 0 && theta_q <= 1);
  require(beta >= 0 && alpha >= 0 && mu >= 0 && rho >= 0 && kappa >= 0 &&
          M > 0 && M <= 1e9 && min_channel >= 0);
  require(eps >= 0 && hysteresis >= 0 && T0 > 0 && cooling > 0 && cooling <= 1);
  require(move_plus >= 0 && move_minus >= 0 && move_both >= 0 &&
          move_flip >= 0 && move_plus + move_minus + move_both + move_flip > 0);
  require(target_density >= 0 && target_density <= 1 && route_cap_h >= 0 &&
          route_cap_v >= 0 && pitch_x >= 0 && pitch_y >= 0);
  require(search == "none" || search == "lite" || search == "sa");
  require(anchor_mode == "exact" || anchor_mode == "soft");
}
} // namespace DreamPlace::native::macro
