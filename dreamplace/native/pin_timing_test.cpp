// BSD-3-Clause. See LICENSE. Assertions remain active under -c opt.
#include "dreamplace/native/placement.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace DreamPlace::native;
namespace {
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
void near(double a, double b, double eps, const char *message) {
  check(std::abs(a - b) <= eps, message);
}
template <typename F> void rejects(F fn) {
  bool failed = false;
  try {
    fn();
  } catch (const std::exception &) {
    failed = true;
  }
  check(failed, "expected invalid timing data rejection");
}
Design network() {
  Design d;
  d.layout = {0, 0, 100, 8, 1, 1};
  d.node_x = {30, 0, 90};
  d.node_y = {2, 2, 2};
  d.node_size_x = {1, 0, 0};
  d.node_size_y = {1, 0, 0};
  d.num_movable = 1;
  d.pin2node = {1, 0, 2};
  d.pin_offset_x = {0, .2, 0};
  d.pin_offset_y = {0, .3, 0};
  d.flat_net2pin = {0, 1, 2};
  d.flat_net2pin_start = {0, 3};
  d.build_maps();
  d.validate();
  return d;
}
PinTiming selective(double coefficient = .25) {
  return {{-1, 0, -1}, {0, 1, 0}, coefficient};
}
Options controls(int iterations = 20) {
  Options o;
  o.max_iterations = o.min_iterations = iterations;
  o.bins_x = 8;
  o.bins_y = 4;
  o.target_density = .5;
  o.add_fillers = false;
  o.legalize = false;
  o.detailed_place = false;
  o.num_threads = 1;
  return o;
}
void gradients_and_selectivity() {
  auto d = network();
  d.pin_timing = selective();
  d.validate();
  const double gamma = .7;
  const auto g = detail::pin_timing_force(d, gamma);
  check(g.x[0] > 0 && g.x[1] < 0, "critical sink/driver did not attract");
  near(g.x[2], 0, 0, "noncritical sink received timing force");
  near(g.y[2], 0, 0, "noncritical sink received timing y force");
  near(g.x[0] + g.x[1] + g.x[2], 0, 1e-12, "x reaction is not conserved");
  near(g.y[0] + g.y[1] + g.y[2], 0, 1e-12, "y reaction is not conserved");
  for (size_t i = 0; i < d.node_x.size(); ++i)
    for (bool horizontal : {true, false}) {
      auto high = d, low = d;
      (horizontal ? high.node_x : high.node_y)[i] += 1e-5;
      (horizontal ? low.node_x : low.node_y)[i] -= 1e-5;
      const double numeric = (detail::pin_timing_force(high, gamma).value -
                              detail::pin_timing_force(low, gamma).value) /
                             2e-5;
      near(numeric, (horizontal ? g.x : g.y)[i], 1e-8,
           "finite difference timing gradient");
    }
  d.pin_timing = PinTiming{};
  const auto cleared = detail::pin_timing_force(d, gamma);
  near(cleared.value, 0, 0, "cleared pin retained timing cost");
  for (double v : cleared.x)
    near(v, 0, 0, "removed pin retained force");
}
void pair_wirelength_oracle() {
  auto d = network();
  d.node_x[0] = .8;
  d.node_y[0] = 2;
  d.pin2node.resize(2);
  d.pin_offset_x = {0, 0};
  d.pin_offset_y = {0, 0};
  d.flat_net2pin = {0, 1};
  d.flat_net2pin_start = {0, 2};
  d.build_maps();
  d.pin_timing = PinTiming{{-1, 0}, {0, 1}, 1};
  d.validate();
  for (double gamma : {.1, .7, 5.}) {
    const auto timing = detail::pin_timing_force(d, gamma);
    const auto ordinary = detail::wirelength(d, gamma, 1);
    near(timing.value, ordinary.value, 1e-12,
         "pin smoothing differs from native two-pin WA");
    for (size_t i = 0; i < d.node_x.size(); ++i) {
      near(timing.x[i], ordinary.x[i], 1e-12,
           "pin gradient differs from native two-pin WA");
      near(timing.y[i], ordinary.y[i], 1e-12,
           "pin y gradient differs from native two-pin WA");
    }
  }
}
void translation_and_units() {
  auto d = network();
  d.pin_timing = selective();
  const auto g = detail::pin_timing_force(d, .7);
  auto moved = d;
  for (double &x : moved.node_x)
    x += 123;
  for (double &y : moved.node_y)
    y -= 73;
  const auto translated = detail::pin_timing_force(moved, .7);
  near(translated.value, g.value, 1e-12, "translation changed timing cost");
  auto scaled = d;
  for (auto *values :
       {&scaled.node_x, &scaled.node_y, &scaled.node_size_x,
        &scaled.node_size_y, &scaled.pin_offset_x, &scaled.pin_offset_y})
    for (double &v : *values)
      v *= 1000;
  const auto units = detail::pin_timing_force(scaled, 700);
  near(units.value, g.value * 1000, 1e-8,
       "length units changed timing objective scale");
  for (size_t i = 0; i < g.x.size(); ++i) {
    near(units.x[i], g.x[i], 1e-12, "length units changed x multiplier");
    near(units.y[i], g.y[i], 1e-12, "length units changed y multiplier");
    near(translated.x[i], g.x[i], 1e-12, "translation changed x force");
  }
}
void numeric_limits() {
  auto d = network();
  d.pin_timing = selective();
  d.pin_offset_x[1] = d.pin_offset_y[1] = 0;
  d.node_x[0] = d.node_x[1];
  d.node_y[0] = d.node_y[1];
  auto g = detail::pin_timing_force(d, 1);
  near(g.value, 0, 0, "coincident pins have timing cost");
  near(g.x[0], 0, 0, "coincident pins have force");
  d.node_x[0] = 1e-10;
  d.node_x[1] = 0;
  d.node_y[0] = d.node_y[1] = 0;
  g = detail::pin_timing_force(d, 1);
  check(g.value > 0, "small timing distance lost to cancellation");
  d.node_x[0] = 1e150;
  g = detail::pin_timing_force(d, 1e150);
  check(std::isfinite(g.value) && std::isfinite(g.x[0]),
        "large finite timing distance overflow");
  auto same = network();
  same.pin2node[1] = same.pin2node[0];
  same.build_maps();
  same.pin_timing = selective();
  same.validate();
  g = detail::pin_timing_force(same, .7);
  for (double v : g.x)
    near(v, 0, 0, "same-node pins exert a translational force");
}
void validation() {
  auto base = network();
  base.pin_timing = selective();
  for (auto mutate : {0, 1, 2, 3, 4, 5, 6, 7}) {
    auto d = base;
    switch (mutate) {
    case 0:
      d.pin_timing->weight.pop_back();
      break;
    case 1:
      d.pin_timing->weight[1] = -1;
      break;
    case 2:
      d.pin_timing->weight[1] = std::numeric_limits<double>::quiet_NaN();
      break;
    case 3:
      d.pin_timing->driver[1] = 1;
      break;
    case 4:
      d.pin_timing->driver[1] = 99;
      break;
    case 5:
      d.pin_timing->coefficient = -1;
      break;
    case 6:
      d.pin_timing->driver[1] = -1;
      break;
    case 7:
      d.net_mask[0] = 0;
      break;
    }
    rejects([&] { d.validate(); });
  }
  auto d = base;
  d.flat_net2pin_start = {0, 1, 3};
  d.net_weights = {1, 1};
  d.net_mask = {0, 1};
  d.build_maps();
  rejects([&] { d.validate(); });
  d = base;
  d.pin_timing->driver[2] = 1;
  d.pin_timing->weight[2] = 1;
  rejects([&] { d.validate(); });
  rejects([&] { detail::pin_timing_force(base, 0); });
}
void gp_integration_and_disabled_replay() {
  auto baseline = network(), zero = baseline, weighted = baseline;
  zero.pin_timing = selective(0);
  weighted.pin_timing = selective(4);
  const auto o = controls();
  const auto control = place(baseline, o), off = place(zero, o),
             on = place(weighted, o);
  check(baseline.node_x == zero.node_x && baseline.node_y == zero.node_y,
        "disabled pin objective changed GP coordinates");
  check(control.history.size() == off.history.size(),
        "disabled pin objective changed convergence history");
  for (size_t i = 0; i < control.history.size(); ++i) {
    near(control.history[i].hpwl, off.history[i].hpwl, 0,
         "disabled objective changed HPWL history");
    near(off.history[i].pin_timing_cost, 0, 0,
         "disabled objective has timing cost");
  }
  check(weighted.node_x[0] < baseline.node_x[0],
        "native GP did not apply selective timing pull");
  bool observed = false;
  for (const auto &sample : on.history)
    observed |= sample.pin_timing_cost > 0;
  check(observed, "GP did not report active timing cost");
  for (size_t i = 1; i < weighted.node_x.size(); ++i)
    near(weighted.node_x[i], network().node_x[i], 0,
         "pin timing moved a fixed node");
}
void checkpoint_and_atomicity() {
  auto d = network();
  auto o = controls(6);
  int calls = 0;
  o.checkpoint_interval = 1;
  o.checkpoint = [&](const Design &snapshot, int iteration) {
    ++calls;
    GpUpdate u;
    if (iteration == 1)
      u.pin_timing = selective(2);
    if (iteration == 2)
      check(snapshot.pin_timing && snapshot.pin_timing->enabled(),
            "checkpoint did not retain pin timing");
    if (iteration == 3)
      u.pin_timing = PinTiming{};
    if (iteration >= 4)
      check(snapshot.pin_timing && !snapshot.pin_timing->enabled(),
            "clear retained stale timing weights");
    return u;
  };
  const auto r = place(d, o);
  check(calls == 6, "checkpoint iteration count");
  check(d.pin_timing && !d.pin_timing->enabled(),
        "cleared objective not written back");
  check(d.net_weights == network().net_weights,
        "pin update changed net criticality weights");
  auto bad = network();
  const auto before = bad;
  o.checkpoint = [](const Design &, int) {
    GpUpdate u;
    u.pin_timing = selective();
    u.pin_timing->driver[1] = 99;
    return u;
  };
  rejects([&] { place(bad, o); });
  check(bad.node_x == before.node_x && bad.node_y == before.node_y &&
            bad.pin_timing == before.pin_timing,
        "invalid update partially committed caller state");
}
void continuation_and_topology() {
  auto d = network();
  d.pin_timing = selective();
  auto o = controls(3);
  GpState state;
  place(d, o, &state);
  check(state.pin_timing == d.pin_timing, "GP state lost pin timing");
  place(d, o, &state);
  auto changed = d;
  changed.pin_timing->weight[1] = .5;
  const auto saved = state;
  rejects([&] { place(changed, o, &state); });
  check(state.x == saved.x && state.pin_timing == saved.pin_timing,
        "failed continuation changed saved objective");
  auto replacement = d;
  o.checkpoint_interval = 1;
  o.checkpoint = [](const Design &snapshot, int iteration) {
    GpUpdate u;
    if (iteration == 1) {
      u.replacement = snapshot;
      u.previous_node = {0, 1, 2};
    }
    if (iteration >= 2)
      check(!snapshot.pin_timing,
            "topology replacement retained obsolete pin identities");
    return u;
  };
  place(replacement, o);
  check(!replacement.pin_timing,
        "replacement did not clear old timing objective");
}
PinTiming gap_timing(bool critical = false) {
  auto timing = selective(2);
  timing.gap_pulls = {{0, 70, 5, .1, critical, {}}, {0, 60, 1, .15, false, {}}};
  if (critical)
    timing.gap_pulls[0].suppressed_sinks = {1};
  return timing;
}
void gap_gradient_translation_and_units() {
  for (bool critical : {false, true}) {
    auto d = network();
    d.pin_timing = gap_timing(critical);
    const auto g = detail::pin_timing_force(d, .7);
    for (size_t i = 0; i < d.node_x.size(); ++i)
      for (bool horizontal : {true, false}) {
        auto high = d, low = d;
        (horizontal ? high.node_x : high.node_y)[i] += 1e-5;
        (horizontal ? low.node_x : low.node_y)[i] -= 1e-5;
        const double numeric = (detail::pin_timing_force(high, .7).value -
                                detail::pin_timing_force(low, .7).value) /
                               2e-5;
        near(numeric, (horizontal ? g.x : g.y)[i], 1e-8,
             "gap finite-difference gradient");
      }
    auto moved = d;
    for (double &x : moved.node_x)
      x += 123;
    for (double &y : moved.node_y)
      y -= 73;
    for (auto &pull : moved.pin_timing->gap_pulls) {
      pull.x += 123;
      pull.y -= 73;
    }
    const auto translated = detail::pin_timing_force(moved, .7);
    near(translated.value, g.value, 1e-10, "gap translation changed objective");
    auto scaled = d;
    for (auto *values :
         {&scaled.node_x, &scaled.node_y, &scaled.node_size_x,
          &scaled.node_size_y, &scaled.pin_offset_x, &scaled.pin_offset_y})
      for (double &v : *values)
        v *= 1000;
    for (auto &pull : scaled.pin_timing->gap_pulls) {
      pull.x *= 1000;
      pull.y *= 1000;
    }
    const auto units = detail::pin_timing_force(scaled, 700);
    near(units.value, g.value * 1000, 1e-7,
         "gap unit scaling changed objective");
    for (size_t i = 0; i < g.x.size(); ++i) {
      near(translated.x[i], g.x[i], 1e-12, "gap translation changed force");
      near(units.x[i], g.x[i], 1e-12, "gap length unit changed x force");
      near(units.y[i], g.y[i], 1e-12, "gap length unit changed y force");
    }
    if (critical)
      near(g.x[0], 0, 0, "suppressed sink retained straight timing force");
    d.pin_timing = PinTiming{};
    near(detail::pin_timing_force(d, .7).value, 0, 0,
         "cleared gap retained cost");
  }
}
void gap_two_pin_oracle() {
  auto d = network();
  d.pin_timing = gap_timing(true);
  d.pin_timing->gap_pulls.resize(1);
  auto pair = network();
  pair.node_x = {0, 70};
  pair.node_y = {2, 5};
  pair.node_size_x = pair.node_size_y = {0, 0};
  pair.num_movable = 0;
  pair.pin2node = {0, 1};
  pair.pin_offset_x = pair.pin_offset_y = {0, 0};
  pair.flat_net2pin = {0, 1};
  pair.flat_net2pin_start = {0, 2};
  pair.build_maps();
  pair.validate();
  for (double gamma : {.1, .7, 5.}) {
    const auto g = detail::pin_timing_force(d, gamma);
    const auto oracle = detail::wirelength(pair, gamma, 1);
    near(g.value, .2 * oracle.value, 1e-12,
         "gap target differs from fixed two-pin WA oracle");
    near(g.x[1], .2 * oracle.x[0], 1e-12, "gap x force differs from WA oracle");
    near(g.y[1], .2 * oracle.y[0], 1e-12, "gap y force differs from WA oracle");
    near(g.x[2], 0, 0, "gap applied force to an unrelated sink");
  }
}
void gap_validation_and_limits() {
  for (int which = 0; which < 13; ++which) {
    auto d = network();
    d.pin_timing = gap_timing(true);
    auto &timing = *d.pin_timing;
    switch (which) {
    case 0:
      timing.gap_pulls[0].driver = 99;
      break;
    case 1:
      timing.gap_pulls[0].driver = 1;
      break;
    case 2:
      timing.gap_pulls[0].fraction = .2;
      break;
    case 3:
      timing.gap_pulls[0].x = std::numeric_limits<double>::quiet_NaN();
      break;
    case 4:
      timing.gap_pulls[0].fraction = 0;
      break;
    case 5:
      timing.weight[1] = 0;
      break;
    case 6:
      timing.gap_pulls[0].critical = false;
      break;
    case 7:
      timing.gap_pulls[0].suppressed_sinks = {1, 1};
      break;
    case 8:
      timing.gap_pulls[0].suppressed_sinks = {2};
      break;
    case 9:
      timing.maximum_gap_fraction = 0;
      break;
    case 10:
      timing.maximum_gap_fraction = 2;
      break;
    case 11:
      timing.driver.clear();
      timing.weight.clear();
      break;
    case 12:
      timing.gap_pulls[0].fraction = std::numeric_limits<double>::infinity();
      break;
    }
    rejects([&] { d.validate(); });
  }
  auto d = network();
  d.pin_timing = gap_timing();
  d.pin_timing->coefficient = 0;
  near(detail::pin_timing_force(d, .7).value, 0, 0,
       "disabled gaps retained force");
  d.pin_timing = gap_timing();
  d.node_x[1] = 1e150;
  check(std::isfinite(detail::pin_timing_force(d, 1e150).value),
        "large finite gap target failed");
}
void gap_priority_and_gp_motion() {
  auto straight = network(), gap = network();
  straight.pin2node = gap.pin2node = {0, 1, 2};
  straight.build_maps();
  gap.build_maps();
  straight.pin_timing = selective(4);
  gap.pin_timing = selective(4);
  gap.pin_timing->gap_pulls = {{0, 70, 2, .25, true, {1}}};
  const auto force = detail::pin_timing_force(gap, .7);
  check(force.x[0] < 0,
        "critical target lost priority to conflicting straight pair");
  near(force.x[1], 0, 0, "suppressed straight pair retained force");
  const auto o = controls();
  place(straight, o);
  const auto r = place(gap, o);
  check(gap.node_x[0] > straight.node_x[0],
        "real GP ignored critical gap direction");
  check(gap.node_x[1] == 0 && gap.node_x[2] == 90, "gap moved fixed cells");
  check(!r.history.empty() && r.history.back().pin_timing_cost > 0,
        "gap GP cost was not reported");
}
void gap_checkpoint_clear_and_atomicity() {
  auto d = network();
  auto o = controls(5);
  o.checkpoint_interval = 1;
  o.checkpoint = [](const Design &snapshot, int iteration) {
    GpUpdate u;
    if (iteration == 1)
      u.pin_timing = gap_timing(true);
    if (iteration == 2)
      check(snapshot.pin_timing == gap_timing(true),
            "checkpoint lost gap payload");
    if (iteration == 3)
      u.pin_timing = PinTiming{};
    if (iteration >= 4)
      check(snapshot.pin_timing->gap_pulls.empty(),
            "clear retained gap identities");
    return u;
  };
  place(d, o);
  auto bad = network();
  const auto before = bad;
  o.checkpoint = [](const Design &, int) {
    GpUpdate u;
    u.pin_timing = gap_timing(true);
    u.pin_timing->gap_pulls[0].fraction = .9;
    return u;
  };
  rejects([&] { place(bad, o); });
  check(bad.node_x == before.node_x && bad.node_y == before.node_y &&
            bad.pin_timing == before.pin_timing,
        "invalid gap checkpoint partially changed the caller");
  d = network();
  d.pin_timing = gap_timing();
  GpState state;
  o = controls(2);
  place(d, o, &state);
  check(state.pin_timing == d.pin_timing, "continuation lost gap payload");
  const auto saved = state;
  d.pin_timing->gap_pulls[0].x += 1;
  rejects([&] { place(d, o, &state); });
  check(state.x == saved.x && state.pin_timing == saved.pin_timing,
        "stale gap continuation modified retained state");
}
} // namespace
int main() {
  try {
    gradients_and_selectivity();
    pair_wirelength_oracle();
    translation_and_units();
    numeric_limits();
    validation();
    gp_integration_and_disabled_replay();
    checkpoint_and_atomicity();
    continuation_and_topology();
    gap_gradient_translation_and_units();
    gap_two_pin_oracle();
    gap_validation_and_limits();
    gap_priority_and_gp_motion();
    gap_checkpoint_clear_and_atomicity();
  } catch (const std::exception &e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  std::cout << "13 pin/gap timing numerical, mapping and GP checks passed\n";
}
