// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>

using namespace DreamPlace::native;
namespace {
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("continuation check failed at " +
                             std::to_string(at.line()));
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 80, 40, 1, 1};
  d.num_movable = 3;
  d.node_x = {10, 20, 30, 0, 80};
  d.node_y = {10, 20, 10, 20, 20};
  d.node_size_x = {2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 0, 0};
  d.node_macro.assign(5, 0);
  d.pin2node = {3, 0, 0, 1, 1, 2, 2, 4};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.build_maps();
  return d;
}
Options options() {
  Options o;
  o.num_threads = 1;
  o.bins_x = o.bins_y = 8;
  o.max_iterations = o.min_iterations = o.schedule_iterations = 50;
  o.legalize = o.detailed_place = false;
  o.filler_max_count = 64;
  o.checkpoint_interval = 5;
  return o;
}
void same_state(const GpState &a, const GpState &b) {
  check(a.x == b.x && a.y == b.y && a.prev_x == b.prev_x &&
        a.prev_y == b.prev_y);
  check(a.old_gx == b.old_gx && a.old_gy == b.old_gy && a.old_ex == b.old_ex &&
        a.old_ey == b.old_ey);
  check(
      a.iteration == b.iteration && a.wirelength_epoch == b.wirelength_epoch &&
      a.signature == b.signature && a.density_weight == b.density_weight &&
      a.step == b.step && a.momentum == b.momentum && a.restart == b.restart &&
      a.filler_replenish_pending == b.filler_replenish_pending);
  check(a.route_model == b.route_model && a.route_hold == b.route_hold &&
        a.route_external_last_iteration == b.route_external_last_iteration &&
        a.route_external_overflow == b.route_external_overflow &&
        a.route_external_peak == b.route_external_peak &&
        a.route_external_unroutable_nets == b.route_external_unroutable_nets &&
        a.route_external_blocked_pins == b.route_external_blocked_pins);
  check(a.pin_timing == b.pin_timing && a.inflation == b.inflation &&
        a.movable_holds == b.movable_holds &&
        a.movable_updates_retired == b.movable_updates_retired &&
        a.checkpoint_invoked == b.checkpoint_invoked &&
        a.last_checkpoint == b.last_checkpoint);
}
void external_route_retains_cadence_and_rolls_back_failures() {
  auto d = fixture();
  auto o = options();
  o.add_fillers = false;
  o.min_iterations = o.schedule_iterations = 1000;
  o.max_iterations = o.min_slice_iterations = 70;
  o.routability = o.improved_routability = true;
  o.route.model = RouteModel::external;
  o.route_external_interval = 50;
  int calls = 0;
  o.route_external = [&](const Design &physical, int bx, int by) {
    check(physical.node_x.size() == 5);
    ++calls;
    RudyResult result;
    result.horizontal.assign(bx * by, 0);
    result.vertical.assign(bx * by, 0);
    result.pin_utilization.assign(bx * by, 0);
    return result;
  };
  auto full = d;
  GpState reference;
  check(place(full, o, &reference).iterations == 70 && calls == 2);
  auto split = d;
  GpState state;
  calls = 0;
  o.max_iterations = o.min_slice_iterations = 30;
  place(split, o, &state);
  check(calls == 1);
  o.max_iterations = o.min_slice_iterations = 40;
  place(split, o, &state);
  check(calls == 2 && split.node_x == full.node_x &&
        split.node_y == full.node_y);
  same_state(reference, state);
  const auto before = state;
  const auto physical = split;
  o.max_iterations = o.min_slice_iterations = 100;
  o.route_external = [&](const Design &, int, int) -> RudyResult {
    throw std::runtime_error("external router failed after retained GP work");
  };
  bool failed = false;
  try {
    place(split, o, &state);
  } catch (const std::runtime_error &) {
    failed = true;
  }
  check(failed && split.node_x == physical.node_x &&
        split.node_y == physical.node_y);
  same_state(before, state);
  // Switching models through a host rebind must sample the requested model,
  // even if the retained state previously selected another estimator.
  auto native = options();
  native.add_fillers = false;
  native.max_iterations = native.min_slice_iterations = 5;
  native.min_iterations = 1000;
  GpRebind rebind;
  place(split, native, &state, &rebind);
  check(state.route_model != RouteModel::external);
  calls = 0;
  o.max_iterations = o.min_slice_iterations = 5;
  o.route_external = [&](const Design &, int bx, int by) {
    ++calls;
    RudyResult result;
    result.horizontal.assign(bx * by, 0);
    result.vertical.assign(bx * by, 0);
    result.pin_utilization.assign(bx * by, 0);
    return result;
  };
  place(split, o, &state, &rebind);
  check(state.route_model == RouteModel::external && calls == 1);
  const auto rebound = state;
  const auto rebound_physical = split;
  o.route_external = [](const Design &, int, int) { return RudyResult{}; };
  failed = false;
  try {
    place(split, o, &state, &rebind);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  check(failed && split.node_x == rebound_physical.node_x &&
        split.node_y == rebound_physical.node_y);
  same_state(rebound, state);
}
void rebind_selects_native_route_model() {
  for (bool deferred : {false, true}) {
    for (auto first : {RouteModel::rudy, RouteModel::rsmt, RouteModel::off}) {
      auto d = fixture();
      auto o = options();
      o.add_fillers = false;
      o.max_iterations = 1;
      o.min_iterations = 1000;
      o.routability = o.improved_routability = true;
      o.route.model = first;
      GpState state;
      place(d, o, &state);
      o.route.model = first == RouteModel::rudy   ? RouteModel::rsmt
                      : first == RouteModel::rsmt ? RouteModel::rudy
                                                  : RouteModel::automatic;
      const GpRebind binding;
      if (deferred) {
        o.max_iterations = 0;
        place(d, o, &state, &binding);
        o.max_iterations = 1;
        place(d, o, &state);
      } else {
        place(d, o, &state, &binding);
      }
      check(state.route_model == (o.route.model == RouteModel::automatic
                                      ? RouteModel::rudy
                                      : o.route.model));
    }
  }
}
void external_route_controls_convergence() {
  for (int mode = 0; mode < 3; ++mode) {
    auto d = fixture();
    auto o = options();
    o.add_fillers = false;
    o.max_iterations = o.min_iterations = o.min_slice_iterations = 20;
    o.routability = true;
    o.route.model = RouteModel::external;
    o.route.hold = mode == 1 ? 2 : 1;
    o.route_external = [mode](const Design &, int bx, int by) {
      RudyResult result;
      result.horizontal.assign(bx * by, mode == 0 ? 2 : 0);
      result.vertical.assign(bx * by, 0);
      result.pin_utilization.assign(bx * by, 0);
      result.overflow = mode == 0 ? 1 : 0;
      result.max_utilization = mode == 0 ? 2 : 0;
      return result;
    };
    const auto result = place(d, o);
    check(result.overflow <= o.stop_overflow);
    check(result.converged == (mode == 2));
  }
}
void deferred_rebind_refills_atomically() {
  auto d = fixture();
  auto o = options();
  o.max_iterations = o.min_slice_iterations = 1;
  o.min_iterations = 1000;
  o.target_density = .004;
  o.filler_width = 2;
  o.filler_max_count = 32;
  o.max_step = 1e-30;
  o.checkpoint_interval = 1;
  bool fail = false;
  o.checkpoint = [&](const Design &, int) {
    if (fail)
      throw std::runtime_error("failed timing after deferred filler refill");
    return GpUpdate{};
  };
  GpState old;
  place(d, o, &old);
  d.layout.xh *= 2;
  d.layout.yh *= 2;
  const GpRebind binding;
  auto direct = d;
  auto reference = old;
  place(direct, o, &reference, &binding);
  check(reference.x.size() > old.x.size());
  auto deferred = d;
  auto state = old;
  o.max_iterations = o.min_slice_iterations = 0;
  place(deferred, o, &state, &binding);
  const auto pending = state;
  const auto physical = deferred;
  o.max_iterations = o.min_slice_iterations = 1;
  fail = true;
  bool rejected = false;
  try {
    place(deferred, o, &state);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected && deferred.node_x == physical.node_x &&
        deferred.node_y == physical.node_y);
  same_state(state, pending);
  fail = false;
  place(deferred, o, &state);
  check(state.x.size() == reference.x.size());
  check(deferred.node_x == direct.node_x && deferred.node_y == direct.node_y);
  same_state(state, reference);
}
void split_equivalence() {
  for (bool fillers : {false, true}) {
    for (int mode = 0; mode < 6; ++mode) {
      auto o = options();
      o.add_fillers = fillers;
      std::vector<int> calls;
      int finished_iteration = 0;
      o.observer = [&](const Design &, PlacementStage stage, int iteration) {
        if (stage == PlacementStage::finished)
          finished_iteration = iteration;
      };
      o.checkpoint = [&](const Design &d, int iteration) {
        calls.push_back(iteration);
        check(d.node_x.size() == 5); // private fillers never reach STA
        GpUpdate u;
        if (mode >= 1) {
          u.net_weights = d.net_weights;
          u.net_weights[1] = iteration % 10 ? 2 : 1;
          u.timing_violation = iteration < 30 ? 1 : 0;
        }
        if (mode >= 2) {
          PinTiming timing;
          timing.driver.assign(8, -1);
          timing.weight.assign(8, 0);
          timing.driver[3] = 2;
          timing.weight[3] = iteration % 10 ? 2 : 1;
          u.pin_timing = timing;
        }
        if (mode == 5) {
          u.pin_offset_x = d.pin_offset_x;
          u.pin_offset_x[1] = iteration % 10 ? .25 : .5;
        }
        if (mode == 4 && iteration == 20)
          u.stop_movable_updates = true;
        return u;
      };
      if (mode >= 3) {
        o.terminal_interval = 5;
        o.terminal_update = [](const Design &, int iteration) {
          return std::vector<TerminalPosition>{{4, 80, 10.003 + iteration}};
        };
        o.movable_interval = 5;
        o.movable_update = [](const Design &, int iteration) {
          return std::vector<MovablePosition>{
              {0, 12.001, 22.003, iteration == 20 ? 7 : 0}};
        };
      }
      auto full = fixture();
      GpState reference;
      check(place(full, o, &reference).iterations == 50);
      const auto expected_calls = calls;
      auto uninterrupted = fixture();
      calls.clear();
      place(uninterrupted, o);
      check(calls == expected_calls && uninterrupted.node_x == full.node_x &&
            uninterrupted.node_y == full.node_y);

      for (int split_at : {5, 17, 20, 21}) {
        calls.clear();
        auto split = fixture();
        GpState state;
        o.max_iterations = split_at;
        check(place(split, o, &state).iterations == split_at);
        check(state.iteration == split_at && finished_iteration == split_at);
        if (fillers)
          check(state.x.size() > split.node_x.size());
        o.max_iterations = 50 - split_at;
        check(place(split, o, &state).iterations == 50 - split_at);
        check(calls == expected_calls && finished_iteration == 50);
        check(full.node_x == split.node_x && full.node_y == split.node_y &&
              full.net_weights == split.net_weights &&
              full.pin_timing == split.pin_timing);
        same_state(reference, state);
      }
    }
  }
}
void checkpoint_milestones() {
  auto o = options();
  o.checkpoint_overflow_milestones = {1, .5};
  o.stop_overflow = 1;
  int calls = 0;
  o.checkpoint = [&](const Design &, int) {
    ++calls;
    return GpUpdate{};
  };
  auto full = fixture(), split = full;
  GpState reference, state;
  const auto result = place(full, o, &reference);
  check(calls == 1 && result.converged);
  calls = 0;
  o.max_iterations = 20;
  place(split, o, &state);
  o.max_iterations = 30;
  const auto resumed = place(split, o, &state);
  check(calls == 1 && resumed.converged);
  check(full.node_x == split.node_x && full.node_y == split.node_y);
  same_state(reference, state);
  check(state.convergence.checkpoint_milestone == 2);
}
void atomic_rejection() {
  auto o = options();
  o.max_iterations = 5;
  auto d = fixture();
  GpState state;
  place(d, o, &state);
  const auto saved_design = d;
  const auto saved_state = state;
  for (int mode = 0; mode < 4; ++mode) {
    o.max_iterations = 10;
    o.checkpoint = [&](const Design &snapshot, int iteration) {
      GpUpdate u;
      if (iteration == 10) {
        u.net_weights = snapshot.net_weights;
        u.net_weights[0] = 2;
      } else if (mode == 0) {
        throw std::runtime_error("STA failed after an accepted update");
      } else if (mode == 1) {
        u.node_size_x = snapshot.node_size_x;
        u.node_size_x[0] = 3;
      } else if (mode == 2) {
        u.replacement = snapshot;
      } else {
        u.density_reservations =
            std::vector<DensityReservation>{{"repair", 20, 20, 2, 2}};
      }
      return u;
    };
    bool rejected = false;
    try {
      place(d, o, &state);
    } catch (const std::exception &) {
      rejected = true;
    }
    check(rejected);
    check(d.node_x == saved_design.node_x && d.node_y == saved_design.node_y &&
          d.net_weights == saved_design.net_weights &&
          d.node_size_x == saved_design.node_size_x);
    same_state(saved_state, state);
  }
}
void refill_rebound_whitespace() {
  for (double density : {.002, .004}) {
    for (int cap : {8, 32}) {
      auto d = fixture();
      auto o = options();
      o.max_iterations = o.min_iterations = o.min_slice_iterations = 1;
      o.schedule_iterations = 50;
      o.target_density = density;
      o.filler_width = 2;
      o.filler_max_count = cap;
      o.fillers_in_whitespace = true;
      o.max_step =
          1e-30; // retain exact particle coordinates for the append audit
      GpState state;
      place(d, o, &state);
      const auto old = state;
      const size_t physical = d.node_x.size();
      const size_t retained = old.x.size() - physical;
      check(density == .002 ? retained == 0 : retained > 0);
      d.layout.xh *= 2;
      d.layout.yh *= 2;
      o.target_density = .004;
      auto repeated_design = d;
      auto repeated_state = old;
      const GpRebind binding;
      const auto result = place(d, o, &state, &binding);
      place(repeated_design, o, &repeated_state, &binding);
      same_state(state, repeated_state);
      check(state.iteration == 2);
      check(state.x.size() > old.x.size());
      check(state.x.size() - physical <= size_t(cap));
      for (size_t i = physical; i < old.x.size(); ++i)
        check(state.x[i] == old.x[i] && state.y[i] == old.y[i]);
      for (const auto *values :
           {&state.x, &state.y, &state.prev_x, &state.prev_y, &state.old_ex,
            &state.old_ey, &state.old_gx, &state.old_gy, &state.inflation})
        check(values->size() == state.x.size());
      const double charge = o.target_density * 160 * 80 - 6;
      check(std::abs(result.effective_filler_area - charge) < 1e-10);
      int calls = 0;
      o.checkpoint = [&](const Design &, int) {
        ++calls;
        return GpUpdate{};
      };
      o.max_iterations = o.min_slice_iterations = 20;
      o.max_step = .2;
      check(place(d, o, &state).iterations == 20);
      check(calls == 4 && state.iteration == 22);
    }
  }
}
void rebind_density_normalization() {
  auto d = fixture();
  auto o = options();
  o.max_iterations = 1;
  GpState state;
  place(d, o, &state);
  state.density_weight =
      1e9; // penalty accumulated under an infeasible old core
  const auto before = state;
  const GpRebind ordinary;
  o.max_iterations = 0;
  place(d, o, &state, &ordinary);
  check(state.density_weight == 1e9 && state.density_initialized);
  state = before;
  d.layout.xh *= 2;
  d.layout.yh *= 2;
  GpRebind recalibrate;
  recalibrate.renormalize_density = true;
  o.density_weight_seed =
      1e8; // explicit re-normalization uses the initial knob
  place(d, o, &state, &recalibrate);
  check(state.density_weight == o.initial_density_weight &&
        !state.density_initialized && state.iteration == before.iteration);
  check(state.x == before.x && state.y == before.y);
  const auto pending = state;
  const auto geometry = d;
  o.max_iterations = 5;
  o.min_slice_iterations = 5;
  o.checkpoint = [](const Design &, int) -> GpUpdate {
    throw std::runtime_error("failed timing after recalibration");
  };
  bool rejected = false;
  try {
    place(d, o, &state);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected && d.node_x == geometry.node_x && d.node_y == geometry.node_y);
  same_state(pending, state);
  check(!state.density_initialized);
  o.checkpoint = {};
  place(d, o, &state);
  check(state.density_initialized && state.density_weight > 0 &&
        state.density_weight != o.initial_density_weight &&
        state.density_weight < 1e8 && state.iteration == before.iteration + 5);
}
void reset_stale_routing_inflation() {
  auto d = fixture();
  auto o = options();
  o.max_iterations = 1;
  GpState state;
  place(d, o, &state);
  for (int i = 0; i < d.num_movable; ++i)
    state.inflation[i] = 2.4;
  const auto before = state;
  const auto physical = d;
  o.max_iterations = 0;
  GpRebind ordinary;
  place(d, o, &state, &ordinary);
  check(state.inflation == before.inflation);
  GpRebind recover;
  recover.reset_routing_inflation = true;
  recover.renormalize_density = true;
  place(d, o, &state, &recover);
  for (int i = 0; i < d.num_movable; ++i)
    check(state.inflation[i] == 1);
  for (size_t i = d.num_movable; i < state.inflation.size(); ++i)
    check(state.inflation[i] == before.inflation[i]);
  check(state.x == before.x && state.y == before.y &&
        state.iteration == before.iteration);
  check(d.node_size_x == physical.node_size_x &&
        d.node_size_y == physical.node_size_y &&
        d.pin_offset_x == physical.pin_offset_x &&
        d.pin_offset_y == physical.pin_offset_y);
  check(!state.density_initialized);
}
void restart_wirelength_schedule() {
  auto d = fixture();
  auto o = options();
  o.add_fillers = false;
  o.max_iterations = o.min_slice_iterations = o.schedule_iterations = 100;
  GpState old;
  check(place(d, o, &old).iterations == 100);
  const auto pose = d;
  o.max_iterations = 0;
  o.min_slice_iterations = 0;
  o.schedule_iterations = 50;
  // Compare fixed work budgets, independently of density convergence at a
  // slice boundary (whose early return intentionally ends that solve).
  o.min_iterations = 1000;
  GpRebind recover;
  recover.restart_wirelength = true;
  recover.renormalize_density = true;
  auto pending = old;
  place(d, o, &pending, &recover);
  check(pending.iteration == 100 && pending.wirelength_epoch == 100);
  check(pending.x == old.x && pending.y == old.y &&
        !pending.density_initialized);
  const auto retained = pending;
  o.max_iterations = o.min_slice_iterations = 50;
  std::vector<int> calls;
  o.checkpoint = [&](const Design &, int iteration) {
    calls.push_back(iteration);
    return GpUpdate{};
  };
  auto full = pose;
  auto reference = pending;
  check(place(full, o, &reference).iterations == 50);
  check(reference.iteration == 150 && reference.wirelength_epoch == 100);
  check(calls.front() == 105 && calls.back() == 150 && calls.size() == 10);
  const auto expected_calls = calls;
  auto split = pose;
  auto state = pending;
  calls.clear();
  o.max_iterations = o.min_slice_iterations = 20;
  check(place(split, o, &state).iterations == 20);
  o.max_iterations = o.min_slice_iterations = 30;
  check(place(split, o, &state).iterations == 30);
  check(calls == expected_calls && split.node_x == full.node_x &&
        split.node_y == full.node_y);
  same_state(reference, state);
  // A restarted wirelength schedule agrees with a fresh optimization at the
  // same physical pose; callback iteration numbers alone remain cumulative.
  auto cold = pose;
  GpState fresh;
  o.max_iterations = o.min_slice_iterations = 50;
  check(place(cold, o, &fresh).iterations == 50);
  check(cold.node_x == full.node_x && cold.node_y == full.node_y);
  // A failed timed slice never publishes its restarted state or moved cells.
  auto failed = pose;
  state = old;
  o.checkpoint = [](const Design &, int) -> GpUpdate {
    throw std::runtime_error("failed timing during wirelength restart");
  };
  bool rejected = false;
  try {
    place(failed, o, &state, &recover);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected && failed.node_x == pose.node_x &&
        failed.node_y == pose.node_y);
  same_state(old, state);
  same_state(retained, pending);
}
void trusted_macro_handoff() {
  auto d = fixture();
  d.node_x[3] = 5;
  d.node_y[3] = 25;
  d.node_size_x[3] = 8;
  d.node_size_y[3] = 5;
  d.node_macro[3] = 1;
  d.node_x[4] = 70;
  d.node_y[4] = 35;
  d.node_size_x[4] = 2;
  d.node_size_y[4] = 1; // a fixed standard cell must remain protected
  d.node_x.push_back(60);
  d.node_y.push_back(30);
  d.node_size_x.push_back(2);
  d.node_size_y.push_back(2);
  d.node_macro.push_back(1); // placement masks also have the macro bit
  d.node_placement_mask.assign(6, PlacementMask::none);
  d.node_placement_mask[5] = PlacementMask::hard;
  d.node_orientation.assign(6, Orientation::N);
  d.build_maps();
  auto o = options();
  o.max_iterations = 5;
  GpState state;
  place(d, o, &state);
  const auto original = d;
  const auto before = state;
  auto grown = d;
  grown.layout.xh = 100;
  grown.layout.yh = 50;
  grown.node_x[3] = 80;
  grown.node_y[3] = 40;
  set_orientation(grown, 3, Orientation::FN);
  const auto handoff = grown;
  o.max_iterations = 0;
  const GpRebind ordinary;
  bool rejected = false;
  try {
    place(grown, o, &state, &ordinary);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected && grown.node_x == handoff.node_x &&
        grown.node_orientation == handoff.node_orientation);
  same_state(before, state);
  GpRebind trusted;
  trusted.allow_fixed_macro_move = true;
  for (int mode = 0; mode < 5; ++mode) {
    auto bad = handoff;
    if (mode == 0)
      bad.node_x[4] += 1; // fixed standard cell
    if (mode == 1)
      set_orientation(bad, 4, Orientation::FS);
    if (mode == 2)
      bad.node_size_x[3] += 1; // macro size remains immutable
    if (mode == 3)
      bad.node_x[5] += 1; // fixed hard placement mask
    if (mode == 4)
      bad.node_macro[3] = 0; // retained macro identity
    const auto candidate = bad;
    rejected = false;
    try {
      place(bad, o, &state, &trusted);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected && bad.node_x == candidate.node_x &&
          bad.node_size_x == candidate.node_size_x &&
          bad.node_orientation == candidate.node_orientation);
    same_state(before, state);
  }
  // A failure after accepting the rebind locally cannot publish either pose
  // or optimizer state to the caller.
  o.max_iterations = o.min_slice_iterations = 5;
  o.checkpoint = [](const Design &, int) -> GpUpdate {
    throw std::runtime_error("STA failed after macro replanning");
  };
  rejected = false;
  try {
    place(grown, o, &state, &trusted);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected && grown.node_x == handoff.node_x);
  same_state(before, state);
  o.checkpoint = {};
  o.max_iterations = o.min_slice_iterations = 0;
  place(grown, o, &state, &trusted);
  check(state.iteration == before.iteration && state.restart &&
        state.momentum == 1);
  check(state.physical.node_x[3] == 80 && state.physical.node_y[3] == 40 &&
        state.physical.node_orientation[3] == Orientation::FN);
  for (int i = 0; i < d.num_movable; ++i)
    check(grown.node_x[i] == original.node_x[i] &&
          grown.node_y[i] == original.node_y[i]);
  check(state.x.size() == before.x.size());
  for (size_t i = grown.node_x.size(); i < before.x.size(); ++i)
    check(state.x[i] == before.x[i] && state.y[i] == before.y[i]);
  check(std::all_of(state.old_gx.begin(), state.old_gx.end(),
                    [](double v) { return v == 0; }));
  o.max_iterations = o.min_slice_iterations = 5;
  check(place(grown, o, &state).iterations == 5 &&
        state.iteration == before.iteration + 5);
  check(grown.node_x[3] == 80 && grown.node_y[3] == 40);
}
} // namespace
int main() {
  try {
    rebind_selects_native_route_model();
    external_route_controls_convergence();
    deferred_rebind_refills_atomically();
    split_equivalence();
    checkpoint_milestones();
    atomic_rejection();
    refill_rebound_whitespace();
    rebind_density_normalization();
    reset_stale_routing_inflation();
    restart_wirelength_schedule();
    external_route_retains_cadence_and_rolls_back_failures();
    trusted_macro_handoff();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
