// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace DreamPlace::native {
bool PinTiming::enabled() const {
  return coefficient > 0 && std::any_of(weight.begin(), weight.end(),
                                        [](double w) { return w > 0; });
}
namespace detail {
namespace {
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
} // namespace
void validate_pin_timing(const Design &d) {
  if (!d.pin_timing)
    return;
  const auto &timing = *d.pin_timing;
  require(std::isfinite(timing.coefficient) && timing.coefficient >= 0,
          "invalid pin timing coefficient");
  require(std::isfinite(timing.maximum_gap_fraction) &&
              timing.maximum_gap_fraction >= 0 &&
              timing.maximum_gap_fraction <= 1,
          "invalid maximum gap-pull fraction");
  if (timing.driver.empty() && timing.weight.empty()) {
    require(timing.gap_pulls.empty(), "gap pulls lack physical pin weights");
    return;
  }
  const size_t pins = d.pin2node.size();
  require(timing.driver.size() == pins && timing.weight.size() == pins &&
              d.pin2net.size() == pins,
          "pin timing array size mismatch");
  require(d.pin_offset_x.size() == pins && d.pin_offset_y.size() == pins &&
              d.node_y.size() == d.node_x.size() &&
              d.net_mask.size() == d.net_weights.size(),
          "invalid pin timing geometry arrays");
  std::vector<int> drivers(d.net_weights.size(), -1);
  for (size_t p = 0; p < pins; ++p) {
    require(d.pin2node[p] >= 0 &&
                static_cast<size_t>(d.pin2node[p]) < d.node_x.size() &&
                d.pin2net[p] >= 0 &&
                static_cast<size_t>(d.pin2net[p]) < drivers.size(),
            "invalid pin timing mapping");
    const double weight = timing.weight[p];
    const int driver = timing.driver[p];
    require(std::isfinite(weight) && weight >= 0, "invalid pin timing weight");
    require(driver >= -1 && (driver < 0 || static_cast<size_t>(driver) < pins),
            "invalid timing driver pin");
    if (driver >= 0) {
      require(static_cast<size_t>(driver) != p &&
                  d.pin2net[p] == d.pin2net[driver],
              "timing driver outside sink net or self driven");
    }
    if (weight > 0) {
      require(driver >= 0, "positive timing weight without a driver");
      const int net = d.pin2net[p];
      require(net >= 0 && static_cast<size_t>(net) < drivers.size() &&
                  d.net_mask[net],
              "timing pin on an inactive net");
      require(drivers[net] < 0 || drivers[net] == driver,
              "conflicting timing drivers on one net");
      require(timing.weight[driver] == 0,
              "timing driver cannot also be a weighted sink");
      drivers[net] = driver;
    }
  }
  if (!timing.gap_pulls.empty()) {
    std::vector<long double> mass(pins), fractions(pins);
    for (size_t p = 0; p < pins; ++p)
      if (timing.weight[p] > 0) {
        mass[timing.driver[p]] += timing.weight[p];
        require(mass[timing.driver[p]] <= std::numeric_limits<double>::max(),
                "gap-pull source weight overflow");
      }
    for (const auto &pull : timing.gap_pulls) {
      require(pull.driver >= 0 && static_cast<size_t>(pull.driver) < pins &&
                  mass[pull.driver] > 0,
              "gap pull lacks an active weighted physical driver");
      require(std::isfinite(pull.x) && std::isfinite(pull.y) &&
                  std::isfinite(pull.fraction) && pull.fraction > 0,
              "invalid gap-pull target/fraction");
      fractions[pull.driver] += pull.fraction;
      require(fractions[pull.driver] <=
                  timing.maximum_gap_fraction *
                      (1 + 16 * std::numeric_limits<double>::epsilon()),
              "gap pulls exceed the per-driver fraction bound");
      require(pull.critical || pull.suppressed_sinks.empty(),
              "outside-window gap cannot suppress straight timing pairs");
      std::set<int> unique;
      for (int p : pull.suppressed_sinks)
        require(p >= 0 && static_cast<size_t>(p) < pins &&
                    timing.weight[p] > 0 && timing.driver[p] == pull.driver &&
                    unique.insert(p).second,
                "gap suppression has stale/duplicate/foreign sink identity");
    }
  }
}
Gradient pin_timing_force(const Design &d, double gamma) {
  require(std::isfinite(gamma) && gamma > 0, "invalid pin timing smoothing");
  validate_pin_timing(d);
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  if (!d.pin_timing || !d.pin_timing->enabled())
    return g;
  const auto &timing = *d.pin_timing;
  const auto coordinate = [&](size_t p, bool horizontal) {
    const int node = d.pin2node[p];
    const double value = horizontal ? d.node_x[node] + d.pin_offset_x[p]
                                    : d.node_y[node] + d.pin_offset_y[p];
    require(std::isfinite(value), "nonfinite timing pin coordinate");
    return value;
  };
  std::vector<bool> suppressed;
  std::vector<double> mass;
  if (!timing.gap_pulls.empty()) {
    suppressed.resize(timing.weight.size());
    mass.resize(timing.weight.size());
    for (size_t p = 0; p < timing.weight.size(); ++p)
      if (timing.weight[p] > 0)
        mass[timing.driver[p]] += timing.weight[p];
    for (const auto &pull : timing.gap_pulls)
      for (int p : pull.suppressed_sinks)
        suppressed[p] = true;
  }
  for (size_t p = 0; p < timing.weight.size(); ++p) {
    if (timing.weight[p] == 0 || (!suppressed.empty() && suppressed[p]))
      continue;
    const size_t driver = timing.driver[p];
    const double scale = timing.coefficient * timing.weight[p];
    require(std::isfinite(scale), "pin timing scale overflow");
    for (bool horizontal : {true, false}) {
      const double delta =
          coordinate(p, horizontal) - coordinate(driver, horizontal);
      require(std::isfinite(delta), "pin timing distance overflow");
      const double z = (delta / gamma) * .5;
      const double t = std::tanh(z);
      // WA(max)-WA(min) for two pins, expressed in translation-invariant
      // differences. Saturated z is the finite L1 limit, never inf*0.
      const double cost = delta * t;
      const double force =
          scale * (t + (std::isfinite(z) ? z * (1 - t * t) : 0));
      g.value += scale * cost;
      auto &values = horizontal ? g.x : g.y;
      if (d.pin2node[p] != d.pin2node[driver]) {
        values[d.pin2node[p]] += force;
        values[d.pin2node[driver]] -= force;
      }
      require(std::isfinite(g.value) && std::isfinite(values[d.pin2node[p]]) &&
                  std::isfinite(values[d.pin2node[driver]]),
              "pin timing force overflow");
    }
  }
  for (const auto &pull : timing.gap_pulls) {
    const double scale =
        timing.coefficient * (mass[pull.driver] * pull.fraction);
    require(std::isfinite(scale), "gap-pull scale overflow");
    for (bool horizontal : {true, false}) {
      const double delta =
          coordinate(pull.driver, horizontal) - (horizontal ? pull.x : pull.y);
      require(std::isfinite(delta), "gap-pull distance overflow");
      const double z = (delta / gamma) * .5, t = std::tanh(z);
      const double force =
          scale * (t + (std::isfinite(z) ? z * (1 - t * t) : 0));
      g.value += scale * delta * t;
      auto &values = horizontal ? g.x : g.y;
      values[d.pin2node[pull.driver]] += force;
      require(std::isfinite(g.value) &&
                  std::isfinite(values[d.pin2node[pull.driver]]),
              "gap-pull force overflow");
    }
  }
  return g;
}
} // namespace detail
} // namespace DreamPlace::native
