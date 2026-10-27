// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/magnet.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("magnet flow assertion at line " +
                             std::to_string(at.line()));
}
Design pair() {
  Design d;
  d.layout = {0, 0, 200, 200, 1, 1};
  d.num_movable = 2;
  d.node_x = {20, 34};
  d.node_y = {20, 20};
  d.node_size_x = d.node_size_y = {10, 10};
  d.node_halo_x = d.node_halo_y = {2, 2};
  d.node_macro = {1, 1};
  d.node_align_class = d.node_snap_class = {0, 7};
  d.build_maps();
  return d;
}
Options controls() {
  Options o;
  o.magnet.enabled = true;
  o.magnet.settle_window = 4;
  o.max_iterations = 200;
  o.min_iterations = 20;
  o.bins_x = o.bins_y = 8;
  o.target_density = .9;
  o.density_upper_bound = true;
  o.add_fillers = o.detailed_place = false;
  o.macro_only_legalize = true;
  return o;
}
void real_gp() {
  auto d = pair();
  auto o = controls();
  int checkpoints = 0, grouped_pulls = 0, snapshots = 0;
  o.checkpoint_interval = 50;
  o.checkpoint = [&](const Design &snapshot, int iteration) {
    check(snapshot.num_movable == 2 && snapshot.node_macro.size() == 2);
    check(iteration % 50 == 0);
    ++checkpoints;
    return GpUpdate{};
  };
  o.macro_pull_interval = 70;
  o.macro_pull_update = [&](const Design &snapshot,
                            int) -> std::optional<std::vector<MacroPull>> {
    if (snapshot.num_movable == 1) {
      ++grouped_pulls;
      check(snapshot.node_size_x[0] == 24 && snapshot.node_size_y[0] == 10);
    }
    return std::nullopt;
  };
  o.observer_interval = 10;
  o.observer = [&](const Design &snapshot, PlacementStage, int) {
    check(snapshot.num_movable == 2 && snapshot.node_x.size() == 2);
    ++snapshots;
  };
  const auto r = place(d, o);
  check(r.magnet && r.magnet->mapping.groups.size() == 1);
  check(r.magnet->stopped_settled && r.converged);
  check(checkpoints && grouped_pulls && snapshots);
  check(macros_are_legal(d));
  check(std::abs(d.node_x[1] - d.node_x[0] - 14) < 1e-9);
  check(d.node_y[0] == d.node_y[1]);
  check(d.node_group_id[0] && d.node_group_id[0] == d.node_group_id[1]);
  int merged = 0;
  for (const auto &round : r.magnet->rounds) {
    if (round.reason != "merged")
      continue;
    ++merged;
    for (const auto &event : r.magnet->events)
      if (event.kind != MagnetEvent::round && event.kind != MagnetEvent::ready)
        check(event.iteration != round.iteration);
  }
  check(merged == 1);
  for (size_t i = 1; i < r.history.size(); ++i)
    check(r.history[i - 1].iteration < r.history[i].iteration);
  check(r.iterations <= o.max_iterations && r.magnet_seconds > 0);
  check(!r.magnet->frames.empty());
  bool saw_group = false, saw_snap = false;
  for (const auto &frame : r.magnet->frames) {
    check(frame.macros.size() == 2);
    check(frame.density.size() == size_t(frame.bins_x * frame.bins_y));
    saw_group |= frame.macros.front().group != 0;
    saw_snap |= frame.event == "snap_disabled";
  }
  check(saw_group && saw_snap);
}
void routing_model_latch_survives_group_restarts() {
  auto d = pair();
  auto o = controls();
  o.routability = true;
  o.route.model = RouteModel::automatic;
  o.route.hold = 1;
  const auto result = place(d, o);
  check(result.magnet && result.magnet->mapping.groups.size() == 1);
  check(result.route_model == RouteModel::rsmt);
  // The empty routing demand switches at the first sample. Contracting a
  // macro pair must not restart automatic routing or replace that switch time.
  check(result.route_switch_iteration == 1);
  check(!result.history.empty());
  for (const auto &sample : result.history) {
    check(sample.route_model == RouteModel::rsmt);
    check(sample.route_switch_iteration == 1);
  }
}
void accepted_routing_probe_is_not_repeated_after_a_merge() {
  auto d = pair();
  auto o = controls();
  o.routability = true;
  o.route.model = RouteModel::automatic;
  o.route.hold = 1000; // the explicit probe, not the hold latch, tests RSMT
  o.route.probe = true;
  const auto result = place(d, o);
  check(result.magnet && result.magnet->mapping.groups.size() == 1);
  check(result.route_model == RouteModel::rudy && result.route_probed);
  check(result.route_switch_iteration == -1);
  check(std::count_if(result.history.begin(), result.history.end(),
                      [](const auto &s) { return s.route_rsmt_seconds > 0; }) ==
        1);
}
void timing_observation_survives_merges_and_snap_relaxation() {
  for (bool recovered : {false, true}) {
    auto d = pair();
    auto o = controls();
    o.max_iterations = 600;
    o.post_density_factor = 1;
    o.checkpoint_after_density = true;
    o.checkpoint_interval = 60;
    int observations = 0;
    bool snapped = false, post_snap_gp = false;
    o.before_legalize = [&](Design &) { snapped = true; };
    o.observer_interval = 1;
    o.observer = [&](const Design &, PlacementStage stage, int) {
      post_snap_gp |= snapped && stage == PlacementStage::global_placement;
    };
    o.macro_merge_check = [&](const std::vector<MagnetMacroFrame> &) {
      return observations ? std::string{} : "await_first_timing";
    };
    o.checkpoint = [&](const Design &, int) {
      check(!snapped); // no STA merely to settle already frozen macros
      GpUpdate update;
      if (++observations == 1)
        update.timing_violation = 1;
      else if (observations == 2 && recovered)
        update.timing_violation = 0;
      // Movement-based deferral may return no new timing observation. Only
      // an actual new estimate may clear the last measured violation.
      return update;
    };
    const auto result = place(d, o);
    check(result.magnet && result.magnet->mapping.groups.size() == 1);
    check(observations > 1 && snapped && post_snap_gp);
    check(result.timing_observed && result.timing_satisfied == recovered);
    // Setup slack is an observation, never a placement feasibility condition.
    check(result.objective_converged);
    check(result.converged);
    // A new stage starts with its own timing observation lifetime.
    o = controls();
    o.magnet.enabled = false;
    const auto next_stage = place(d, o);
    check(!next_stage.timing_observed && next_stage.timing_satisfied);
  }
}
void density_continuation_deadline_is_stage_wide() {
  auto d = pair();
  auto o = controls();
  o.post_density_factor = 1;
  const auto result = place(d, o);
  check(result.magnet && result.magnet->mapping.groups.size() == 1);
  // The first density-feasible sample precedes the merge. Group contraction
  // and the final frozen-macro relaxation must not grant a fresh deadline.
  if (result.density_iteration != 20 || result.continuation_deadline != 40)
    throw std::runtime_error("density continuation reset to " +
                             std::to_string(result.density_iteration) + "/" +
                             std::to_string(result.continuation_deadline) +
                             " instead of 20/40");
}
void caps_and_invalid() {
  auto o = controls();
  auto d = pair();
  o.max_iterations = 2;
  const auto capped = place(d, o);
  check(capped.iterations == 2 && !capped.converged &&
        !capped.magnet->stopped_settled);
  // Virtual fillers extend optimizer arrays, but never the physical macro
  // mask. Adaptive macro momentum checks must not index that mask for fillers.
  auto filled = pair();
  o.add_fillers = true;
  o.filler_max_count = 32;
  const auto filled_result = place(filled, o);
  check(filled_result.filler_count > 0 && filled_result.filler_count <= 32);
  check(filled.node_x.size() == 2);
  o = controls();
  o.checkpoint_interval = 2; // each applied event resets the motion history
  o.checkpoint = [](const Design &, int) { return GpUpdate{}; };
  d = pair();
  const auto deferred = place(d, o);
  check(deferred.magnet->mapping.groups.empty() && !deferred.converged);
  d = pair();
  o = controls();
  o.legalize = false;
  try {
    (void)place(d, o);
    check(false);
  } catch (const std::invalid_argument &) {
  }
}
void stable_pair_is_rechecked_when_cells_leave_gap() {
  auto d = pair();
  d.node_x[1] = 38;
  d.num_movable = 3;
  d.node_x.push_back(32);
  d.node_y.push_back(18);
  d.node_size_x.push_back(4);
  d.node_size_y.push_back(14);
  d.node_halo_x.push_back(0);
  d.node_halo_y.push_back(0);
  d.node_macro.push_back(0);
  d.node_align_class.push_back(-1);
  d.node_snap_class.push_back(-1);
  d.build_maps();
  detail::MagnetRuntime runtime(controls());
  for (int i = 0; i <= 4; ++i) {
    runtime.observe(d, d.node_x, d.node_y, i);
    GpSample sample;
    sample.iteration = i;
    sample.overflow = .8;
    sample.hpwl = 0;
    runtime.sample(d, sample);
  }
  check(runtime.pending.empty() && runtime.report->rounds.size() == 1);
  check(std::any_of(runtime.report->rounds.back().search.candidates.begin(),
                    runtime.report->rounds.back().search.candidates.end(),
                    [](const auto &c) { return c.reason == "gap_density"; }));
  // The macros never move, and no feedback event occurs. Only the cell leaves
  // the corridor. Macro-only caching would incorrectly keep this pair apart.
  d.node_x[2] = 100;
  runtime.observe(d, d.node_x, d.node_y, 5);
  GpSample sample;
  sample.iteration = 5;
  sample.overflow = .8;
  sample.hpwl = 0;
  runtime.sample(d, sample);
  check(runtime.report->rounds.size() == 2 && runtime.pending.size() == 1);
}
void search_rearms_after_geometry_changes() {
  auto d = pair();
  d.node_x[1] = 80; // a settled pair initially outside the gap limit
  detail::MagnetRuntime runtime(controls());
  for (int i = 0; i <= 10; ++i) {
    runtime.observe(d, d.node_x, d.node_y, i);
    GpSample sample;
    sample.iteration = i;
    sample.overflow = .1;
    sample.hpwl = 0;
    runtime.sample(d, sample);
  }
  check(runtime.pending.empty() && runtime.report->rounds.size() == 7);
  check(runtime.report->rounds[0].reason == "no_eligible_pair");
  d.node_x[1] = 34;
  for (int i = 11; i <= 15; ++i) {
    runtime.observe(d, d.node_x, d.node_y, i);
    GpSample sample;
    sample.iteration = i;
    sample.overflow = .1;
    sample.hpwl = 0;
    runtime.sample(d, sample);
    if (i < 15)
      check(runtime.pending.empty());
  }
  check(runtime.report->rounds.size() == 12 && runtime.pending.size() == 1);
  check(runtime.report->rounds.back().iteration == 15);
}
void unrelated_macro_motion_does_not_block_a_stable_pair() {
  auto d = pair();
  d.num_movable = 3;
  d.node_x.push_back(100);
  d.node_y.push_back(100);
  d.node_size_x.push_back(10);
  d.node_size_y.push_back(10);
  d.node_halo_x.push_back(2);
  d.node_halo_y.push_back(2);
  d.node_macro.push_back(1);
  d.node_align_class.push_back(9);
  d.node_snap_class.push_back(9);
  d.build_maps();
  auto o = controls();
  // This check isolates local selection after a deliberately loose global gate.
  o.magnet.start_displacement = 1;
  detail::MagnetRuntime runtime(o);
  GpSample sample;
  for (int i = 0; i <= 4; ++i) {
    d.node_x[2] = i % 2 ? 110 : 100;
    runtime.observe(d, d.node_x, d.node_y, i);
    sample.iteration = i;
    sample.overflow = .8;
    sample.hpwl = 0;
    runtime.sample(d, sample);
    if (i < 4)
      check(runtime.pending.empty());
  }
  check(sample.macro_move_max > .05);
  check(runtime.pending.size() == 1);
  check(runtime.pending[0].members == (std::vector<int>{0, 1}));
}
void rejected_host_check_does_not_restart_gp() {
  auto d = pair();
  d.node_x[1] = 38;
  auto reference = d;
  auto o = controls();
  int calls = 0;
  o.macro_merge_check = [&](const std::vector<MagnetMacroFrame> &members) {
    ++calls;
    check(members.size() == 2 && members[0].source == 0 &&
          members[1].source == 1);
    check(members[1].x - members[0].x == 14 && members[1].y == members[0].y);
    return std::string{"routing_track_phase"};
  };
  const auto rejected = place(d, o);
  check(calls > 0 && rejected.magnet->mapping.groups.empty());
  for (const auto &round : rejected.magnet->rounds) {
    check(round.reason == "no_eligible_pair");
    check(std::any_of(round.search.candidates.begin(),
                      round.search.candidates.end(), [](const auto &c) {
                        return c.reason == "routing_track_phase" && !c.selected;
                      }));
  }
  o.macro_merge_check = {};
  o.magnet.max_gap = 0; // same solve, pair excluded by its positive gap
  const auto excluded = place(reference, o);
  check(excluded.magnet->mapping.groups.empty());
  check(rejected.iterations == excluded.iterations);
  check(d.node_x == reference.node_x && d.node_y == reference.node_y);
}
void host_checks_receive_flat_original_meta_macro_members() {
  auto d = pair();
  d.num_movable = 4;
  d.node_x = {20, 34, 20, 34};
  d.node_y = {20, 20, 34, 34};
  d.node_size_x = d.node_size_y = {10, 10, 10, 10};
  d.node_macro = {1, 1, 1, 1};
  d.node_halo_x = d.node_halo_y = {2, 2, 2, 2};
  d.node_align_class = d.node_snap_class = {0, 1, 2, 3};
  d.build_maps();
  auto o = controls();
  bool saw_meta = false;
  o.macro_merge_check = [&](const std::vector<MagnetMacroFrame> &members) {
    if (members.size() == 4)
      saw_meta = true;
    std::vector<int> ids;
    for (const auto &m : members) {
      check(m.width == 10 && m.height == 10 && m.orientation == Orientation::N);
      check(m.x == 0 || m.x == 14);
      check(m.y == 0 || m.y == 14);
      ids.push_back(m.source);
    }
    std::sort(ids.begin(), ids.end());
    check(std::adjacent_find(ids.begin(), ids.end()) == ids.end());
    return std::string{};
  };
  const auto result = place(d, o);
  check(saw_meta && result.magnet->mapping.groups.size() == 1);
  check(result.magnet->mapping.groups.front().members.size() == 4);
}
void invalid_batch_does_not_interrupt_the_solver() {
  auto d = pair();
  d.node_x[1] = 38;
  // The two-halo merge pose would move both members out of these bounds.
  // They remain movable macros, so ordinary geometry search finds the pair.
  d.node_boxes = {{20, 20, 20, 20}, {38, 20, 38, 20}};
  auto reference = d;
  auto o = controls();
  const auto rejected = place(d, o);
  check(rejected.magnet->mapping.groups.empty());
  check(std::any_of(rejected.magnet->rounds.begin(),
                    rejected.magnet->rounds.end(), [](const auto &r) {
                      return r.reason ==
                             "macro merge violates member position bounds";
                    }));
  // A rejection must have the same solver lifetime as excluding the pair
  // before search: no momentum/filler restart and no lost convergence window.
  o.magnet.max_gap = 0;
  const auto excluded = place(reference, o);
  check(rejected.iterations == excluded.iterations);
  check(rejected.converged == excluded.converged);
  check(d.node_x == reference.node_x && d.node_y == reference.node_y);
  for (const auto &e : rejected.magnet->events)
    check(e.kind != MagnetEvent::round);
}
void global_readiness_latches_without_resetting_local_motion() {
  auto d = pair();
  d.num_movable = 3;
  d.node_x.push_back(100);
  d.node_y.push_back(100);
  d.node_size_x.push_back(10);
  d.node_size_y.push_back(10);
  d.node_halo_x.push_back(2);
  d.node_halo_y.push_back(2);
  d.node_macro.push_back(1);
  d.node_align_class.push_back(9);
  d.node_snap_class.push_back(9);
  d.build_maps();
  auto o = controls();
  o.magnet.start_displacement = .1;
  detail::MagnetRuntime runtime(o);
  GpSample sample;
  for (int i = 0; i <= 4; ++i) {
    d.node_x[2] = i % 2 ? 102 : 100;
    runtime.observe(d, d.node_x, d.node_y, i);
    sample.iteration = i;
    sample.overflow = 0; // density cannot open the motion gate
    runtime.sample(d, sample);
    check(runtime.pending.empty() && runtime.report->rounds.empty());
  }
  sample.iteration = 5;
  sample.overflow = .9; // density cannot block a globally stable window
  d.node_x[2] = 101;    // exactly 10%, above the 5% local threshold
  runtime.observe(d, d.node_x, d.node_y, 5);
  runtime.sample(d, sample);
  check(runtime.report->start_iteration == 5);
  check(runtime.report->start_observed_displacement == .1);
  check(runtime.pending.size() ==
        1); // completed local history survived opening
  check(runtime.pending.front().members == (std::vector<int>{0, 1}));
  check(runtime.report->events.size() == 1);
  check(runtime.report->events.front().kind == MagnetEvent::ready);
  sample.iteration = 6;
  d.node_x[2] = 104;
  runtime.observe(d, d.node_x, d.node_y, 6);
  runtime.sample(d, sample);
  check(sample.macro_move_max > .1);
  check(runtime.pending.size() == 1); // later motion cannot close the latch
  check(runtime.report->events.size() == 1);
  check(runtime.report->overflow_history.size() == 7);
  check(runtime.report->motion_history.size() == 7);
  check(!runtime.report->motion_history.front().movement.window_complete);
  check(runtime.report->motion_history[4].movement.window_complete);
  for (double invalid : {-.5, 1.1}) {
    o.magnet.start_displacement = invalid;
    bool rejected = false;
    try {
      o.magnet.validate();
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected);
  }
}
void before_snap_preserves_ordinary_gp_and_does_not_restart_it() {
  auto baseline = pair();
  auto delayed = baseline;
  auto o = controls();
  o.magnet.enabled = false;
  o.legalize = false;
  o.macro_only_legalize = false;
  o.macro_pull_interval = 7;
  o.macro_pull_update =
      [](const Design &,
         int iteration) -> std::optional<std::vector<MacroPull>> {
    if (iteration == 7)
      return std::vector<MacroPull>{{0, 28, 27, .05}, {1, 42, 27, .05}};
    return std::nullopt;
  };
  const auto ordinary = place(baseline, o);
  o.magnet.enabled = true;
  o.magnet.before_snap = true;
  o.legalize = true;
  o.macro_only_legalize = true;
  const auto late = place(delayed, o);
  check(late.iterations == ordinary.iterations);
  check(late.history.size() == ordinary.history.size());
  for (size_t i = 0; i < late.history.size(); ++i) {
    check(late.history[i].iteration == ordinary.history[i].iteration);
    check(late.history[i].hpwl == ordinary.history[i].hpwl);
    check(late.history[i].overflow == ordinary.history[i].overflow);
  }
  check(late.magnet && late.magnet->mapping.groups.size() == 1);
  check(std::count_if(late.magnet->events.begin(), late.magnet->events.end(),
                      [](const auto &event) {
                        return event.kind == MagnetEvent::pull;
                      }) == 1);
  bool saw_before = false, merged = false;
  for (const auto &frame : late.magnet->frames) {
    if (frame.event == "before_magnet" && !saw_before) {
      saw_before = true;
      for (int i = 0; i < 2; ++i) {
        check(frame.macros[i].x == baseline.node_x[i]);
        check(frame.macros[i].y == baseline.node_y[i]);
      }
    }
    if (frame.event == "magnet")
      merged = true;
    if (merged)
      check(frame.event != "global placement");
  }
  check(saw_before && merged);
  check(late.magnet->motion_history.size() == size_t(late.iterations + 1));
  check(std::any_of(late.magnet->motion_history.begin(),
                    late.magnet->motion_history.end(), [](const auto &sample) {
                      return sample.center_step_valid &&
                             sample.center_step_sum > 0;
                    }));
  for (const auto &round : late.magnet->rounds)
    check(round.iteration == ordinary.iterations);
  check(macros_are_legal(delayed));
}
void absolute_motion_counts_physical_members_and_excludes_merge_jumps() {
  auto d = pair();
  auto o = controls();
  std::vector<double> observed;
  std::vector<bool> validity;
  o.macro_step_observer = [&](int, double sum, bool valid) {
    observed.push_back(sum);
    validity.push_back(valid);
  };
  detail::MagnetRuntime runtime(o);
  auto x = d.node_x, y = d.node_y;
  runtime.observe(d, x, y, 0);
  check(!runtime.report->motion_history.back().center_step_valid);
  x[0] += 3;
  y[0] += 4;
  y[1] += 2;
  runtime.observe(d, x, y, 1);
  check(runtime.report->motion_history.back().center_step_valid);
  check(runtime.report->motion_history.back().center_step_sum == 7);
  runtime.event(1, MagnetEvent::timing);
  y[0] += 1;
  x[1] += 3;
  y[1] += 4;
  runtime.observe(d, x, y, 2);
  check(runtime.report->motion_history.back().center_step_valid);
  check(runtime.report->motion_history.back().center_step_sum == 6);
  auto grouped = pair();
  RigidMacroMap mapping;
  const std::vector<MacroMerge> requests{{{0, 1}, {{20, 20}, {34, 20}}}};
  check(merge_macros(grouped, mapping, requests).applied);
  runtime.report->mapping = mapping;
  runtime.event(2, MagnetEvent::round);
  runtime.observe(grouped, grouped.node_x, grouped.node_y, 2);
  check(!runtime.report->motion_history.back().center_step_valid);
  grouped.node_x[0] += 3;
  grouped.node_y[0] += 4;
  runtime.observe(grouped, grouped.node_x, grouped.node_y, 3);
  const auto &sample = runtime.report->motion_history.back();
  check(sample.center_step_valid && sample.center_step_sum == 10);
  check(sample.physical_macros == 2 && sample.gp_macros == 1);
  check(!sample.movement.window_complete);
  check(observed == std::vector<double>({0, 7, 6, 0, 10}));
  check(validity == std::vector<bool>({false, true, true, false, true}));
}
void merge_retains_only_untouched_warpspace_targets() {
  auto grouped = pair();
  grouped.num_movable = 3;
  grouped.node_x.push_back(100);
  grouped.node_y.push_back(100);
  grouped.node_size_x.push_back(12);
  grouped.node_size_y.push_back(10);
  grouped.node_halo_x.push_back(2);
  grouped.node_halo_y.push_back(2);
  grouped.node_macro.push_back(1);
  grouped.node_align_class.push_back(9);
  grouped.node_snap_class.push_back(9);
  grouped.build_maps();
  auto o = controls();
  o.magnet.enabled = false;
  o.max_iterations = o.min_iterations = 40;
  o.legalize = o.macro_only_legalize = false;
  auto runtime = std::make_shared<detail::MagnetRuntime>(o);
  runtime->macro_pulls = {{0, 25, 25, .1}, {1, 39, 25, .1}, {2, 150, 140, .5}};
  runtime->macro_pull_changed = 5;
  RigidMacroMap mapping;
  const std::vector<MacroMerge> requests{{{0, 1}, {{20, 20}, {34, 20}}}};
  const auto merged = merge_macros(grouped, mapping, requests);
  check(merged.applied);
  runtime->remap_pulls(merged.previous_node);
  check(runtime->macro_pulls.size() == 1);
  const auto kept = runtime->macro_pulls.front();
  check(kept.node == 1 && kept.x == 150 && kept.y == 140 && kept.weight == .5);
  check(runtime->macro_pull_changed == 5);
  runtime->report->mapping = mapping;
  auto unanchored = grouped;
  o.iteration_offset = 20;
  o.magnet_runtime = runtime;
  o.macro_pull_interval = 1;
  o.macro_pull_update = [](const Design &,
                           int) -> std::optional<std::vector<MacroPull>> {
    return std::nullopt;
  };
  (void)place(grouped, o);
  check(runtime->macro_pulls.size() == 1 && runtime->macro_pull_changed == 5);
  // Retaining the target must retain its actual force, not just metadata.
  auto empty_runtime = std::make_shared<detail::MagnetRuntime>(o);
  empty_runtime->report->mapping = mapping;
  o.magnet_runtime = empty_runtime;
  (void)place(unanchored, o);
  check(std::hypot(grouped.node_x[1] - unanchored.node_x[1],
                   grouped.node_y[1] - unanchored.node_y[1]) > 1e-6);
  // An explicit empty replacement still clears targets at its own boundary.
  o.magnet_runtime = runtime;
  o.iteration_offset = 60;
  o.macro_pull_update = [](const Design &,
                           int) -> std::optional<std::vector<MacroPull>> {
    return std::vector<MacroPull>{};
  };
  (void)place(grouped, o);
  check(runtime->macro_pulls.empty() && runtime->macro_pull_changed >= 61);
}
} // namespace
int main() {
  try {
    real_gp();
    routing_model_latch_survives_group_restarts();
    accepted_routing_probe_is_not_repeated_after_a_merge();
    density_continuation_deadline_is_stage_wide();
    timing_observation_survives_merges_and_snap_relaxation();
    caps_and_invalid();
    search_rearms_after_geometry_changes();
    stable_pair_is_rechecked_when_cells_leave_gap();
    unrelated_macro_motion_does_not_block_a_stable_pair();
    rejected_host_check_does_not_restart_gp();
    host_checks_receive_flat_original_meta_macro_members();
    invalid_batch_does_not_interrupt_the_solver();
    global_readiness_latches_without_resetting_local_motion();
    before_snap_preserves_ordinary_gp_and_does_not_restart_it();
    absolute_motion_counts_physical_members_and_excludes_merge_jumps();
    merge_retains_only_untouched_warpspace_targets();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
