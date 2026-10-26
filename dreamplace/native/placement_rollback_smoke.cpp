// BSD-3-Clause. See LICENSE. Assertions remain active under -c opt.
#include "dreamplace/native/placement.hpp"

#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace DreamPlace::native;
namespace {
void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
template <typename F> void rejects(F fn) {
  bool failed = false;
  try {
    fn();
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  check(failed, "expected checkpoint rejection");
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
Options controls(int iterations) {
  Options o;
  o.max_iterations = o.min_iterations = iterations;
  o.bins_x = 8;
  o.bins_y = 4;
  o.target_density = .5;
  o.add_fillers = true;
  o.filler_area_fraction = .25;
  o.legalize = o.detailed_place = false;
  o.num_threads = 1;
  // Hold physical smoothing fixed so independent short/long runs have exactly
  // the same force schedule; rejected work must not contaminate continuation.
  o.wirelength_gamma_start = o.wirelength_gamma_end = 1;
  o.checkpoint_interval = 10;
  o.checkpoint_rollback = true;
  return o;
}
void same(const Design &a, const Design &b) {
  check(a.pin_timing == b.pin_timing, "rollback lost timing objective");
  check(a.node_x.size() == b.node_x.size(),
        "rollback changed physical topology");
  for (size_t i = 0; i < a.node_x.size(); ++i) {
    check(std::abs(a.node_x[i] - b.node_x[i]) < 1e-12,
          "restored x trajectory differs");
    check(std::abs(a.node_y[i] - b.node_y[i]) < 1e-12,
          "restored y trajectory differs");
  }
}
void resumed_trajectory(bool gaps = false) {
  auto control = network(), trial = network();
  std::map<int, Design> expected, actual;
  auto c = controls(20), t = controls(30);
  c.checkpoint_interval = t.checkpoint_interval = 5;
  c.timing_spread_weight = t.timing_spread_weight = .5;
  const auto timing = [&](double coefficient = .25, double target = 70) {
    auto value = selective(coefficient);
    if (gaps)
      value.gap_pulls = {{0, target, 2, .25, true, {1}}};
    return value;
  };
  c.checkpoint = [&](const Design &d, int iteration) {
    expected[iteration] = d;
    GpUpdate u;
    if (iteration == 10) {
      u.pin_timing = timing();
      u.net_weights = {4};
      u.action = GpCheckpointAction::retain;
    }
    return u;
  };
  t.checkpoint = [&](const Design &d, int iteration) {
    actual[iteration] = d;
    GpUpdate u;
    if (iteration == 10) {
      u.pin_timing = timing();
      u.net_weights = {4};
      u.action = GpCheckpointAction::retain;
    }
    if (iteration == 15) {
      u.pin_timing = timing(9, 20);
      u.net_weights = {40};
    }
    if (iteration == 20)
      u.action = GpCheckpointAction::restore;
    return u;
  };
  auto cr = place(control, c), tr = place(trial, t);
  check(cr.filler_count > 0 && tr.filler_count == cr.filler_count,
        "fixture lacks real virtual fillers");
  same(actual.at(10), expected.at(10));
  check(std::abs(actual.at(20).node_x[0] - expected.at(20).node_x[0]) > 1e-6,
        "rejected objective did not perturb the trial");
  same(actual.at(30), expected.at(20));
  same(trial, control);
  check(tr.iterations == 30 && tr.rollback_history.size() == 1,
        "rollback rewound spent budget");
  check(tr.rollback_history[0].iteration == 20 &&
            tr.rollback_history[0].restored_iteration == 10,
        "wrong checkpoint restored");
  check(trial.node_x[1] == 0 && trial.node_x[2] == 90,
        "fixed coordinates changed");
}
void final_restore_and_objective() {
  auto d = network();
  auto o = controls(20);
  o.checkpoint_interval = 5;
  Design retained;
  o.checkpoint = [&](const Design &p, int iteration) {
    GpUpdate u;
    if (iteration == 10) {
      retained = p;
      retained.pin_timing = selective();
      u.pin_timing = selective();
      u.action = GpCheckpointAction::retain;
    }
    if (iteration == 15) {
      u.pin_timing = selective(9);
      u.net_weights = {4};
    }
    if (iteration == 20)
      u.action = GpCheckpointAction::restore;
    return u;
  };
  auto r = place(d, o);
  same(d, retained);
  check(d.net_weights == retained.net_weights, "net weights not restored");
  check(r.iterations == 20 && r.selected_iteration == 10,
        "final restore mislabeled iteration");
}
void repeated_rejection_budget() {
  auto d = network();
  auto o = controls(50);
  Design retained;
  o.checkpoint = [&](const Design &p, int iteration) {
    GpUpdate u;
    u.action = iteration == 10 ? GpCheckpointAction::retain
                               : GpCheckpointAction::restore;
    if (iteration == 10)
      retained = p;
    return u;
  };
  auto r = place(d, o);
  same(d, retained);
  check(r.iterations == 50 && r.selected_iteration == 10 &&
            r.rollback_history.size() == 4,
        "rejections extended budget");
}
void validation_and_atomicity() {
  // Existing positional aggregate clients retain their field order.
  GpUpdate legacy{false, std::nullopt, {},           {}, {},
                  {},    {},           std::nullopt, {}, std::nullopt};
  check(legacy.action == GpCheckpointAction::update,
        "legacy aggregate changed default action");
  const auto original = network();
  auto d = original;
  auto o = controls(20);
  o.checkpoint = [](const Design &, int) {
    GpUpdate u;
    u.action = GpCheckpointAction::restore;
    return u;
  };
  rejects([&] { place(d, o); });
  same(d, original);
  o.checkpoint = [](const Design &, int i) {
    GpUpdate u;
    u.action =
        i == 10 ? GpCheckpointAction::retain : GpCheckpointAction::restore;
    if (i == 20)
      u.pin_timing = selective();
    return u;
  };
  rejects([&] { place(d, o); });
  same(d, original);
  o.checkpoint_rollback = false;
  rejects([&] { place(d, o); });
  same(d, original);
  o.checkpoint_rollback = true;
  o.local_density_weight = .2;
  rejects([&] { place(d, o); });
  same(d, original);
}
} // namespace
int main() {
  try {
    resumed_trajectory();
    resumed_trajectory(true);
    final_restore_and_objective();
    repeated_rejection_budget();
    validation_and_atomicity();
    std::cout << "5 native full-state rollback checks passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
