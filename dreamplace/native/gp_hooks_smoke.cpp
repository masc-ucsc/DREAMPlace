// BSD-3-Clause. See LICENSE.
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>

#include "dreamplace/native/placement.hpp"
using namespace DreamPlace::native;
void check(bool ok,
           std::source_location where = std::source_location::current()) {
  if (!ok) {
    throw std::runtime_error("GP hook assertion failed at line " +
                             std::to_string(where.line()));
  }
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 80, 40, 1, 1};
  d.num_movable = 3;
  d.node_x = {10, 20, 30, 0, 80};
  d.node_y = {10, 20, 10, 20, 20};
  d.node_size_x = {2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 0, 0};
  d.pin2node = {3, 0, 0, 1, 1, 2, 2, 4};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.build_maps();
  return d;
}
int main() {
  try {
    // Relative virtual pulls must track ordinary wire forces without needing
    // the unrelated local-density pull objective. Global net-weight units
    // should not change this one-macro two-step trajectory.
    const auto anchored = [](double net_weight) {
      Design d;
      d.layout = {0, 0, 100, 100, 1, 1};
      d.num_movable = 1;
      d.node_x = {40, 90};
      d.node_y = {40, 45};
      d.node_size_x = d.node_size_y = {10, 0};
      d.node_macro = {1, 0};
      d.pin2node = {0, 1};
      d.pin_offset_x = d.pin_offset_y = {5, 0};
      d.flat_net2pin = {0, 1};
      d.flat_net2pin_start = {0, 2};
      d.net_weights = {net_weight};
      d.build_maps();
      Options o;
      o.max_iterations = 2;
      o.bins_x = o.bins_y = 8;
      o.initial_density_weight = 1e-30;
      o.add_fillers = o.legalize = o.detailed_place = false;
      o.macro_pull_interval = 1;
      o.macro_pull_update =
          [](const Design &,
             int iteration) -> std::optional<std::vector<MacroPull>> {
        if (iteration == 1) {
          return std::vector<MacroPull>{{0, 30, 80, 1}};
        }
        return std::nullopt;
      };
      place(d, o);
      return std::pair{d.node_x[0], d.node_y[0]};
    };
    const auto light = anchored(100), heavy = anchored(10000);
    check(std::abs(light.first - heavy.first) < 1e-7);
    check(std::abs(light.second - heavy.second) < 1e-7);
    check(light.second > 45); // the anchor actually pulls upward
    auto full = fixture(), split = full;
    Options o;
    o.max_iterations = 40;
    o.schedule_iterations = 40;
    o.legalize = o.detailed_place = false;
    o.bins_x = o.bins_y = 8;
    o.stop_overflow = 0;
    place(full, o);
    GpState state;
    o.max_iterations = 20;
    place(split, o, &state);
    check(state.iteration == 20);
    place(split, o, &state);
    check(state.iteration == 40);
    check(full.node_x == split.node_x && full.node_y == split.node_y);
    const auto old = split;
    const auto history = state;
    split.node_size_x[0] = 3;
    bool bad = false;
    try {
      place(split, o, &state);
    } catch (const std::invalid_argument &) {
      bad = true;
    }
    check(bad && state.x == history.x && state.iteration == history.iteration);
    split = old;
    // A late refinement can hold final smoothing without reheating. Split
    // continuation is identical, and changing its smoothing invalidates state.
    auto settled = fixture(), settled_split = settled;
    auto smoothing_options = o;
    smoothing_options.wirelength_gamma_start = .05;
    smoothing_options.wirelength_gamma_end = .05;
    smoothing_options.max_iterations = 40;
    place(settled, smoothing_options);
    smoothing_options.max_iterations = 20;
    GpState settle_state;
    place(settled_split, smoothing_options, &settle_state);
    place(settled_split, smoothing_options, &settle_state);
    check(settled.node_x == settled_split.node_x &&
          settled.node_y == settled_split.node_y);
    smoothing_options.wirelength_gamma_start = .1;
    bool changed_start = false;
    try {
      place(settled_split, smoothing_options, &settle_state);
    } catch (const std::invalid_argument &) {
      changed_start = true;
    }
    check(changed_start);
    check(std::abs(density_overflow(settled, 8, 8, .8)) < 1e-9);
    auto soft = fixture();
    Options soft_options;
    soft_options.max_iterations = 20;
    soft_options.bins_x = soft_options.bins_y = 8;
    soft_options.add_fillers = soft_options.legalize =
        soft_options.detailed_place = false;
    GpState soft_state;
    place(soft, soft_options, &soft_state);
    soft_options.density_upper_bound = true;
    bool soft_changed = false;
    try {
      place(soft, soft_options, &soft_state);
    } catch (const std::invalid_argument &) {
      soft_changed = true;
    }
    check(soft_changed);
    soft = fixture();
    check(std::isfinite(place(soft, soft_options).hpwl));
    // A satisfied soft upper bound must not repel an isolated cell with no data
    // nets merely because a different part of the die contains a fixed macro.
    Design vacant;
    vacant.layout = {0, 0, 100, 100, 1, 1};
    vacant.num_movable = 1;
    vacant.node_x = {80, 20, 80.5};
    vacant.node_y = {50, 20, 50.5};
    vacant.node_size_x = {1, 40, 0};
    vacant.node_size_y = {1, 60, 0};
    vacant.node_macro = {0, 1, 0};
    vacant.flat_net2pin_start = {0};
    vacant.build_maps();
    auto vacancy_options = soft_options;
    vacancy_options.exact_macro_density = true;
    vacancy_options.target_density = .5;
    vacancy_options.bins_x = vacancy_options.bins_y = 10;
    const auto optimum = vacant;
    place(vacant, vacancy_options);
    check(vacant.node_x == optimum.node_x && vacant.node_y == optimum.node_y);
    auto near_macro = optimum;
    near_macro.node_x[0] = 42;
    near_macro.node_size_x[1] = 20;
    const auto near_optimum = near_macro;
    auto exact_cells = vacancy_options;
    exact_cells.exact_cell_density = true;
    place(near_macro, exact_cells);
    check(near_macro.node_x == near_optimum.node_x &&
          near_macro.node_y == near_optimum.node_y);
    auto exact_continuation = fixture();
    GpState exact_state;
    place(exact_continuation, soft_options, &exact_state);
    auto changed_footprint = soft_options;
    changed_footprint.exact_cell_density = true;
    bool footprint_changed = false;
    try {
      place(exact_continuation, changed_footprint, &exact_state);
    } catch (const std::invalid_argument &) {
      footprint_changed = true;
    }
    check(footprint_changed);
    auto changed_obstacles = soft_options;
    changed_obstacles.fixed_macro_overlap_weight = 5;
    bool obstacle_changed = false;
    try {
      place(exact_continuation, changed_obstacles, &exact_state);
    } catch (const std::invalid_argument &) {
      obstacle_changed = true;
    }
    check(obstacle_changed);
    auto changed_gamma = soft_options;
    changed_gamma.wirelength_gamma_end = .01;
    bool gamma_changed = false;
    try {
      place(exact_continuation, changed_gamma, &exact_state);
    } catch (const std::invalid_argument &) {
      gamma_changed = true;
    }
    check(gamma_changed);
    auto patient = fixture();
    auto patient_options = soft_options;
    patient_options.min_iterations = 80;
    patient_options.max_iterations = 80;
    patient_options.stop_overflow = 1;
    const auto patient_result = place(patient, patient_options);
    check(patient_result.iterations == 80);
    check(patient_result.selected_iteration >= 40);
    for (const auto &sample : patient_result.history) {
      if (sample.iteration >= 40 &&
          sample.overflow <= patient_options.stop_overflow) {
        check(patient_result.hpwl <= sample.hpwl + 1e-8);
      }
    }
    // A density-feasible but still moving objective must not stop at 40.
    auto settling = fixture();
    auto settle_options = soft_options;
    settle_options.stop_overflow = 1;
    settle_options.max_iterations = 240;
    settle_options.post_density_factor = 4;
    settle_options.convergence_tolerance = 0;
    settling.node_boxes.resize(5);
    for (int i = 0; i < 5; ++i) {
      settling.node_boxes[i] = {settling.node_x[i], settling.node_y[i],
                                settling.node_x[i], settling.node_y[i]};
    }
    const auto stable_result = place(settling, settle_options);
    check(stable_result.density_iteration == 40);
    check(stable_result.continuation_deadline == 200);
    check(stable_result.iterations == 80 && stable_result.converged);
    check(stable_result.objective_converged &&
          !stable_result.continuation_limit);
    // Fractional continuation is a hard iteration budget, including between
    // ordinary 20-step convergence samples.
    auto short_options = settle_options;
    short_options.post_density_factor = .25;
    short_options.routability = true;
    short_options.route_capacity_h = short_options.route_capacity_v = 1e-6;
    auto short_design = settling;
    const auto short_result = place(short_design, short_options);
    check(short_result.density_iteration == 40);
    check(short_result.continuation_deadline == 50 &&
          short_result.iterations == 50);
    check(short_result.continuation_limit && !short_result.converged);
    // Unmet timing is reported, but never vetoes geometric convergence.
    settle_options.checkpoint_interval = 20;
    settle_options.checkpoint = [](const Design &, int) {
      GpUpdate update;
      update.timing_violation = 1;
      return update;
    };
    auto failing_timing = settling;
    const auto bounded = place(failing_timing, settle_options);
    check(bounded.iterations == 80 && !bounded.continuation_limit);
    check(bounded.converged && bounded.timing_observed &&
          !bounded.timing_satisfied);
    // Objective updates reset stability, never the continuation deadline.
    settle_options.checkpoint = [](const Design &d, int iteration) {
      GpUpdate update;
      update.net_weights.assign(d.net_weights.size(), iteration % 40 ? 1. : 2.);
      return update;
    };
    auto changing = settling;
    const auto changing_result = place(changing, settle_options);
    check(changing_result.iterations == 200 &&
          changing_result.continuation_limit);
    // External continuation and topology checkpoints are separate APIs.
    // Use persistent routing overflow to verify the split-solve deadline.
    settle_options.checkpoint = {};
    settle_options.routability = true;
    settle_options.route_capacity_h = settle_options.route_capacity_v = 1e-6;
    auto continued = settling;
    GpState bounded_state;
    settle_options.max_iterations = 60;
    place(continued, settle_options, &bounded_state);
    check(bounded_state.convergence.deadline == 200);
    settle_options.max_iterations = 180;
    const auto resumed = place(continued, settle_options, &bounded_state);
    check(bounded_state.iteration == 200 && resumed.continuation_limit);
    auto incompatible = settle_options;
    incompatible.post_density_factor = 3;
    bool changed_stop = false;
    try {
      place(continued, incompatible, &bounded_state);
    } catch (const std::invalid_argument &) {
      changed_stop = true;
    }
    check(changed_stop);
    settle_options.checkpoint = {};
    settle_options.max_iterations = 60;
    auto capped_design = settling;
    const auto capped = place(capped_design, settle_options);
    check(capped.iterations == 60 && !capped.converged &&
          !capped.continuation_limit);
    auto route_design = settling;
    settle_options.max_iterations = 240;
    settle_options.routability = true;
    settle_options.route_capacity_h = settle_options.route_capacity_v = 1e-6;
    const auto congested = place(route_design, settle_options);
    // A static routing plateau may hand off before the continuation deadline;
    // it must still report the unsatisfied objective without convergence.
    check(congested.iterations < 200 && congested.route_stalled &&
          !congested.continuation_limit && !congested.converged);
    check(congested.continuation_deadline == 200);
    auto boxed = fixture();
    boxed.node_boxes.assign(5, {0, 0, 78, 39});
    boxed.node_boxes[0] = {4, 6, 4, 6};
    int snapshots = 0;
    o.observer_interval = 1;
    o.observer = [&](const Design &d, PlacementStage, int) {
      ++snapshots;
      check(d.node_x[0] == 4 && d.node_y[0] == 6);
    };
    place(boxed, o);
    check(snapshots > 2);
    auto order = fixture();
    order.node_x[0] = 50;
    order.node_x[1] = 10;
    order.order_edges.push_back({0, 1, Axis::x, 0, {false, 1}, "left"});
    o.observer = {};
    o.order_weight = 1000;
    o.max_iterations = 100;
    o.add_fillers = false;
    place(order, o);
    check(order.node_x[0] + 2 <= order.node_x[1] + .1);
    // Unused classes and alignment knobs must not perturb disabled results.
    auto disabled = fixture(), original = disabled;
    Options base;
    base.max_iterations = 40;
    base.legalize = base.detailed_place = false;
    base.bins_x = base.bins_y = 8;
    place(original, base);
    disabled.node_macro = {1, 1, 1, 0, 0};
    disabled.node_align_class = {0, 0, 0, -1, -1};
    base.macro_align_global = true; // zero weight must remain bit-identical
    base.macro_align_sigma = .25;
    base.macro_align_overlap_boost = 2;
    place(disabled, base);
    check(original.node_x == disabled.node_x &&
          original.node_y == disabled.node_y);
    base.macro_align_global = false;
    auto aligned = fixture();
    aligned.node_macro = {1, 1, 1, 0, 0};
    aligned.node_align_class = {0, 0, 0, -1, -1};
    base.macro_align_weight = 25;
    base.macro_align_start = 0;
    base.macro_align_full = 60;
    base.stop_overflow = 1;
    base.max_iterations = 80;
    auto alignment = place(aligned, base);
    check(alignment.iterations >= 60 && !alignment.history.empty());
    GpState aligned_state;
    base.max_iterations = 20;
    aligned = fixture();
    aligned.node_macro = {1, 1, 1, 0, 0};
    aligned.node_align_class = {0, 0, 0, -1, -1};
    place(aligned, base, &aligned_state);
    base.macro_align_sigma = .5;
    bool mismatch = false;
    try {
      place(aligned, base, &aligned_state);
    } catch (const std::invalid_argument &) {
      mismatch = true;
    }
    check(mismatch);
    // Active global alignment has deterministic continuation and a salted mode.
    auto global = fixture();
    global.node_macro = {1, 1, 1, 0, 0};
    global.node_align_class = {0, 0, 0, -1, -1};
    global.node_y[2] += .4;
    auto global_resumed = global;
    auto global_options = base;
    global_options.macro_align_global = true;
    global_options.macro_align_sigma = .25;
    global_options.schedule_iterations = 40;
    global_options.max_iterations = 40;
    global_options.stop_overflow = 0;
    const auto global_result = place(global, global_options);
    check(global_result.macro_alignment_seconds > 0);
    bool force_seen = false;
    for (const auto &sample : global_result.history) {
      force_seen |= sample.align_force > 0;
    }
    check(force_seen);
    global_options.max_iterations = 20;
    GpState global_state;
    place(global_resumed, global_options, &global_state);
    place(global_resumed, global_options, &global_state);
    check(global.node_x == global_resumed.node_x &&
          global.node_y == global_resumed.node_y);
    global_options.macro_align_global = false;
    bool global_changed = false;
    try {
      place(global_resumed, global_options, &global_state);
    } catch (const std::invalid_argument &) {
      global_changed = true;
    }
    check(global_changed);
    auto route = fixture();
    route.node_inflatable = {1, 1, 1, 0, 0};
    auto sizes = route.node_size_x;
    o.order_weight = 0;
    o.routability = true;
    o.route_capacity_h = o.route_capacity_v = .01;
    place(route, o);
    check(route.node_size_x == sizes);
    // Timing checkpoints do not split/reseed GP, and a no-op is bit-identical.
    auto untouched = fixture(), observed = untouched;
    untouched.node_macro = observed.node_macro = {0, 0, 0, 0, 0};
    Options cp;
    cp.max_iterations = 120;
    cp.bins_x = cp.bins_y = 8;
    cp.add_fillers = cp.legalize = cp.detailed_place = false;
    cp.target_density = .001;
    cp.stop_overflow = 0;
    place(untouched, cp);
    std::vector<int> calls;
    cp.checkpoint_interval = 50;
    cp.checkpoint = [&](const Design &, int iteration) {
      calls.push_back(iteration);
      return GpUpdate{};
    };
    place(observed, cp);
    check(calls == std::vector<int>({50, 100}));
    check(observed.node_x == untouched.node_x &&
          observed.node_y == untouched.node_y);
    // Warm-up is independent of density and leaves no-op GP unchanged.
    // Even a forced final checkpoint cannot bypass the in-loop warm-up.
    for (bool density_gate : {false, true}) {
      auto warming = fixture();
      warming.node_macro = {0, 0, 0, 0, 0};
      auto warm_options = cp;
      warm_options.checkpoint_start = 75;
      warm_options.checkpoint_final = true;
      warm_options.checkpoint_after_density = density_gate;
      calls.clear();
      place(warming, warm_options);
      check(calls == (density_gate ? std::vector<int>{120}
                                   : std::vector<int>{100, 120}));
      check(warming.node_x == untouched.node_x &&
            warming.node_y == untouched.node_y);
      warm_options.checkpoint_start = 200;
      calls.clear();
      place(warming, warm_options);
      check(calls.empty());
    }
    // Above-target density defers feedback, but the non-interval cap is timed.
    // The skipped no-op checkpoints must not disturb the optimization path.
    auto gated = fixture();
    gated.node_macro = {0, 0, 0, 0, 0};
    auto gate_options = cp;
    gate_options.checkpoint_after_density = true;
    calls.clear();
    const auto gated_result = place(gated, gate_options);
    check(calls == std::vector<int>({120}));
    check(gated.node_x == untouched.node_x && gated.node_y == untouched.node_y);
    check(gated_result.checkpoint_history.size() == 3);
    check(!gated_result.checkpoint_history[0].invoked &&
          !gated_result.checkpoint_history[1].invoked);
    check(gated_result.checkpoint_history.back().forced &&
          gated_result.checkpoint_history.back().overflow > 0);
    // An interval later than density convergence must not let GP finish
    // untimed.
    auto early = fixture();
    gate_options.target_density = .5;
    gate_options.stop_overflow = 1;
    gate_options.checkpoint_interval = 200;
    gate_options.max_iterations = 240;
    calls.clear();
    check(place(early, gate_options).iterations == 200);
    check(calls == std::vector<int>({200}));
    // The bounded post-density deadline also forces feedback off cadence.
    gate_options.post_density_factor = 4;
    gate_options.checkpoint_interval = 75;
    gate_options.checkpoint = [&](const Design &, int iteration) {
      calls.push_back(iteration);
      GpUpdate update;
      update.timing_violation = 1;
      return update;
    };
    calls.clear();
    auto deadline = settling;
    const auto gate_deadline = place(deadline, gate_options);
    check(calls == std::vector<int>({75}));
    check(!gate_deadline.continuation_limit && gate_deadline.converged &&
          !gate_deadline.checkpoint_history.back().forced);
    auto off_sample_cap = settling;
    gate_options.max_iterations = 113;
    gate_options.checkpoint_interval = 200;
    calls.clear();
    const auto sampled_cap = place(off_sample_cap, gate_options);
    check(calls == std::vector<int>({113}));
    check(sampled_cap.timing_observed && !sampled_cap.timing_satisfied &&
          !sampled_cap.converged);
    gate_options.checkpoint_interval = 0;
    calls.clear();
    place(deadline, gate_options);
    check(calls.empty());
    // The gate sees halo-padded density, while the host sees physical macros.
    auto gated_halo = fixture();
    gated_halo.node_macro = {1, 0, 0, 0, 0};
    gated_halo.node_halo_x = {20, 0, 0, 0, 0};
    gated_halo.node_halo_y = {10, 0, 0, 0, 0};
    gate_options = cp;
    gate_options.target_density = .01;
    gate_options.checkpoint_after_density = true;
    gate_options.checkpoint = [&](const Design &snapshot, int iteration) {
      check(snapshot.node_size_x[0] == 2 && snapshot.node_size_y[0] == 1);
      calls.push_back(iteration);
      return GpUpdate{};
    };
    calls.clear();
    const auto halo_gate = place(gated_halo, gate_options);
    check(calls == std::vector<int>({120}));
    check(halo_gate.checkpoint_history[0].overflow > 0);
    auto resized = fixture();
    resized.node_macro = {0, 0, 0, 0, 0};
    cp.checkpoint = [](const Design &d, int) {
      GpUpdate u;
      u.net_weights = d.net_weights;
      u.net_weights[0] = 5;
      u.node_size_x = d.node_size_x;
      u.node_size_x[1] = 3;
      u.pin_offset_x = d.pin_offset_x;
      u.pin_offset_x[3] = .25;
      return u;
    };
    place(resized, cp);
    check(resized.net_weights[0] == 5 && resized.node_size_x[1] == 3 &&
          resized.pin_offset_x[3] == .25);
    auto haloed = fixture();
    haloed.node_macro = {1, 0, 0, 0, 0};
    haloed.node_halo_x = {2, 0, 0, 0, 0};
    haloed.node_halo_y = {2, 0, 0, 0, 0};
    const auto resize_callback = cp.checkpoint;
    cp.checkpoint = [&](const Design &d, int iteration) {
      check(d.node_size_x[0] == 2 && d.node_size_y[0] == 1);
      return resize_callback(d, iteration);
    };
    place(haloed, cp);
    check(haloed.node_size_x[0] == 2 && haloed.node_size_x[1] == 3);
    check(haloed.node_halo_x[0] == 2);
    auto filled_resize = fixture();
    filled_resize.node_macro = {0, 0, 0, 0, 0};
    auto filler_resize_options = cp;
    filler_resize_options.checkpoint = resize_callback;
    filler_resize_options.add_fillers = true;
    filler_resize_options.target_density = .5;
    filler_resize_options.filler_max_count = 64;
    filler_resize_options.filler_area_fraction = .5;
    filler_resize_options.min_iterations = filler_resize_options.max_iterations;
    const auto resize_result = place(filled_resize, filler_resize_options);
    check(resize_result.filler_count > 0 && resize_result.filler_count <= 64);
    check(filled_resize.node_size_x[1] == 3 &&
          filled_resize.node_x.size() == 5);
    auto invalid = fixture();
    invalid.node_macro = {1, 0, 0, 0, 0};
    const auto saved = invalid;
    cp.checkpoint = [](const Design &d, int) {
      GpUpdate u;
      u.node_size_x = d.node_size_x;
      u.node_size_x[0] += 1;
      return u;
    };
    bad = false;
    try {
      place(invalid, cp);
    } catch (const std::invalid_argument &) {
      bad = true;
    }
    check(bad && invalid.node_x == saved.node_x &&
          invalid.node_size_x == saved.node_size_x);
    // A topology checkpoint may insert, reorder and remove only standard cells.
    // Keep a haloed macro, terminal and color candidate through both changes.
    auto topology = fixture();
    topology.node_macro = {1, 0, 0, 0, 0};
    topology.node_halo_x = {2, 0, 0, 0, 0};
    topology.node_halo_y = {2, 0, 0, 0, 0};
    topology.node_centroid_color = {0, 0, 0, -1, -1};
    topology.centroid_color_names = {"bank"};
    for (size_t i = 0; i < topology.node_x.size(); ++i) {
      topology.centroid_node_area.push_back(topology.node_size_x[i] *
                                            topology.node_size_y[i]);
    }
    CentroidGroup bank;
    bank.name = "bank";
    bank.colors = {0};
    topology.centroid_groups = {bank};
    cp.centroid_weight = .25;
    cp.centroid_repulsion_weight = 1;
    cp.centroid_adaptive_repulsion = true;
    cp.checkpoint_interval = 50;
    calls.clear();
    cp.checkpoint = [&](const Design &snapshot, int iteration) {
      calls.push_back(iteration);
      check(snapshot.node_size_x[0] == 2 && snapshot.node_halo_x[0] == 2);
      Design next = snapshot;
      GpUpdate u;
      if (iteration == 50) {
        // Insert a standard cell before the fixed suffix and split one net.
        const int at = next.num_movable;
        next.node_x.insert(next.node_x.begin() + at, 8);
        next.node_y.insert(next.node_y.begin() + at, 2);
        next.node_size_x.insert(next.node_size_x.begin() + at, 2);
        next.node_size_y.insert(next.node_size_y.begin() + at, 1);
        next.node_macro.insert(next.node_macro.begin() + at, 0);
        next.node_halo_x.insert(next.node_halo_x.begin() + at, 0);
        next.node_halo_y.insert(next.node_halo_y.begin() + at, 0);
        next.node_centroid_color.insert(next.node_centroid_color.begin() + at,
                                        0);
        next.centroid_node_area.insert(next.centroid_node_area.begin() + at, 2);
        for (auto &p : next.pin2node) {
          if (p >= at) {
            ++p;
          }
        }
        next.pin2node.push_back(at);
        next.pin_offset_x.push_back(1);
        next.pin_offset_y.push_back(.5);
        next.flat_net2pin.push_back(next.pin2node.size() - 1);
        next.flat_net2pin_start.back()++;
        ++next.num_movable;
        for (int i = 0; i < int(next.node_x.size()); ++i) {
          u.previous_node.push_back(i == at ? -1 : i - (i > at));
        }
      } else {
        // Remove the same inserted buffer, leaving the original pins/net
        // intact.
        const int at = next.num_movable - 1;
        next.node_x.erase(next.node_x.begin() + at);
        next.node_y.erase(next.node_y.begin() + at);
        next.node_size_x.erase(next.node_size_x.begin() + at);
        next.node_size_y.erase(next.node_size_y.begin() + at);
        next.node_macro.erase(next.node_macro.begin() + at);
        next.node_halo_x.erase(next.node_halo_x.begin() + at);
        next.node_halo_y.erase(next.node_halo_y.begin() + at);
        next.node_centroid_color.erase(next.node_centroid_color.begin() + at);
        next.centroid_node_area.erase(next.centroid_node_area.begin() + at);
        next.pin2node.pop_back();
        next.pin_offset_x.pop_back();
        next.pin_offset_y.pop_back();
        next.flat_net2pin.pop_back();
        next.flat_net2pin_start.back()--;
        for (auto &p : next.pin2node) {
          if (p > at) {
            --p;
          }
        }
        --next.num_movable;
        for (int i = 0; i < int(next.node_x.size()); ++i) {
          u.previous_node.push_back(i + (i >= at));
        }
      }
      next.build_maps();
      u.replacement = std::move(next);
      return u;
    };
    // An uncapped density budget scales with physical area, beyond 65k
    // particles.
    {
      auto sparse = fixture();
      sparse.layout.xh = sparse.layout.yh = 1000;
      Options resolution;
      resolution.add_fillers = true;
      resolution.filler_max_count = 0;
      resolution.bins_x = resolution.bins_y = 16;
      resolution.target_density = .5;
      resolution.max_iterations = resolution.min_iterations = 1;
      resolution.legalize = resolution.detailed_place = false;
      resolution.routability = false;
      const auto small = place(sparse, resolution);
      check(small.filler_count > 65536 && !small.filler_count_limited);
      check(small.filler_size_x <= 1000. / 16 * .25 &&
            small.filler_size_y <= 1000. / 16 * .25);
      sparse.layout.xh = 2000;
      const auto large = place(sparse, resolution);
      check(large.filler_count > small.filler_count * 1.9);
      resolution.filler_max_count = 64;
      const auto guarded = place(sparse, resolution);
      check(guarded.filler_count <= 64 && guarded.filler_count_limited);
      check(std::abs(guarded.filler_area - large.filler_area) <=
            small.filler_size_x * small.filler_size_y * 2);
    }
    // An empty/deferred checkpoint must not retire jumps. The first real
    // timing update can retire them and replace topology in the same callback,
    // including through the halo adapter and filler implementation.
    {
      auto jumping = topology;
      auto jumping_options = cp;
      jumping_options.add_fillers = true;
      jumping_options.filler_max_count = 128;
      jumping_options.checkpoint_interval = 25;
      jumping_options.movable_interval = 10;
      std::vector<int> jumps;
      jumping_options.movable_update = [&](const Design &, int iteration) {
        jumps.push_back(iteration);
        return std::vector<MovablePosition>{};
      };
      jumping_options.checkpoint = [&](const Design &snapshot, int iteration) {
        if (iteration % 50) {
          return GpUpdate{};
        }
        auto update = cp.checkpoint(snapshot, iteration);
        update.stop_movable_updates = true;
        return update;
      };
      place(jumping, jumping_options);
      check(jumps == std::vector<int>({10, 20, 30, 40, 50}));
      check(jumping.num_movable == topology.num_movable);
      calls.clear();
    }
    auto with_virtual = topology;
    auto virtual_options = cp;
    virtual_options.add_fillers = true;
    virtual_options.target_density = .5;
    virtual_options.min_iterations = virtual_options.max_iterations;
    virtual_options.filler_max_count = 128;
    virtual_options.local_density_weight = 2;
    virtual_options.macro_cell_density_weight = 2;
    virtual_options.local_density_bins = 32;
    virtual_options.observer = [&](const Design &snapshot, PlacementStage,
                                   int) {
      snapshot.validate();
      check(snapshot.node_x.size() == 5 || snapshot.node_x.size() == 6);
    };
    const auto virtual_result = place(with_virtual, virtual_options);
    check(virtual_result.filler_count > 0 &&
          virtual_result.filler_count <= 128);
    check(with_virtual.node_x.size() == 5 && with_virtual.num_movable == 3);
    check(calls == std::vector<int>({50, 100}));
    calls.clear();
    place(topology, cp);
    check(calls == std::vector<int>({50, 100}));
    check(topology.node_x.size() == 5 && topology.num_movable == 3);
    check(topology.node_x[3] == fixture().node_x[3] &&
          topology.node_x[4] == fixture().node_x[4]);
    check(topology.node_size_x[0] == 2 && topology.node_halo_x[0] == 2);
    check(topology.centroid_groups[0].name == "bank");
    // Final-GP memory schedules survive insertion/removal and fixed-node
    // index shifts, including halo conversion and virtual filler particles.
    for (bool grow : {false, true}) {
      auto scheduled = topology;
      scheduled.node_macro[0] = 0;
      scheduled.node_macro[3] = 1;
      scheduled.node_size_x[3] = scheduled.node_size_y[3] = 4;
      auto scheduled_options = virtual_options;
      scheduled_options.macro_cell_density_weight = 0;
      scheduled_options.centroid_weight = 0;
      scheduled_options.centroid_repulsion_weight = 0;
      scheduled_options.centroid_adaptive_repulsion = false;
      scheduled_options.max_iterations = scheduled_options.min_iterations = 225;
      scheduled_options.post_density_factor = 0;
      scheduled_options.vanished_memory = grow;
      if (!grow) {
        scheduled_options.macro_transparency.warmup = 20;
        scheduled_options.macro_transparency.period = 20;
        scheduled_options.macro_transparency.cycles = 2;
      }
      scheduled_options.checkpoint = [&](const Design &d, int iteration) {
        return iteration <= 100 ? cp.checkpoint(d, iteration) : GpUpdate{};
      };
      calls.clear();
      place(scheduled, scheduled_options);
      check(calls == std::vector<int>({50, 100}));
      scheduled.validate();
      check(scheduled.node_x.size() == 5 && scheduled.num_movable == 3);
      check(scheduled.node_x[3] == topology.node_x[3]);
      check(scheduled.node_y[3] == topology.node_y[3]);
      check(scheduled.node_size_x[3] == 4 && scheduled.node_size_y[3] == 4);
    }
    // Invalid mappings fail atomically, including through the halo wrapper.
    const auto prior = topology;
    cp.checkpoint = [](const Design &snapshot, int) {
      GpUpdate u;
      u.replacement = snapshot;
      u.previous_node = {1, 1, 2, 3,
                         4}; // loses the macro, repeats a standard cell
      return u;
    };
    bad = false;
    try {
      place(topology, cp);
    } catch (const std::invalid_argument &) {
      bad = true;
    }
    check(bad && topology.node_x == prior.node_x &&
          topology.node_size_x == prior.node_size_x);
    // A malformed replacement throws before padded_design indexes it.
    auto malformed = fixture();
    malformed.node_macro = {0, 0, 0, 0, 0};
    malformed.node_halo_x = {2, 0, 0, 0, 0};
    malformed.node_halo_y = {2, 0, 0, 0, 0};
    const auto malformed_before = malformed;
    cp.checkpoint = [](const Design &snapshot, int) {
      GpUpdate u;
      u.replacement = snapshot;
      u.replacement->node_halo_x.resize(2);
      u.previous_node = {0, 1, 2, 3, 4};
      return u;
    };
    bad = false;
    try {
      place(malformed, cp);
    } catch (const std::invalid_argument &) {
      bad = true;
    }
    check(bad && malformed.node_x == malformed_before.node_x &&
          malformed.node_halo_x == malformed_before.node_halo_x);
    // Halo padding of boxes and order edges must not leak into the caller's
    // Design or the checkpoint snapshot.
    auto shifted = fixture();
    shifted.node_macro = {0, 0, 0, 0, 0};
    shifted.node_halo_x = {2, 0, 0, 0, 0};
    shifted.node_halo_y = {2, 0, 0, 0, 0};
    shifted.node_boxes.assign(5, {0, 0, 78, 39});
    shifted.node_boxes[0] = {6, 6, 60, 30};
    shifted.order_edges.push_back({0, 1, Axis::x, 0, {false, 1}, "left"});
    const auto shifted_before = shifted;
    cp.checkpoint = [&](const Design &snapshot, int) {
      check(snapshot.node_boxes[0].xl == 6 && snapshot.node_boxes[0].yh == 30 &&
            snapshot.order_edges[0].separation == 0);
      return GpUpdate{};
    };
    place(shifted, cp);
    check(shifted.node_boxes[0].xl == 6 && shifted.node_boxes[0].xh == 60 &&
          shifted.node_boxes[0].yl == 6 && shifted.node_boxes[0].yh == 30 &&
          shifted.order_edges[0].separation == 0 &&
          shifted.node_size_x == shifted_before.node_size_x &&
          shifted.node_halo_x == shifted_before.node_halo_x);
    // A resize under halos keeps centroid areas physical, not halo-padded.
    auto resized_halo = fixture();
    resized_halo.node_macro = {0, 0, 0, 0, 0};
    resized_halo.node_halo_x = {2, 0, 0, 0, 0};
    resized_halo.node_halo_y = {2, 0, 0, 0, 0};
    for (size_t i = 0; i < resized_halo.node_x.size(); ++i) {
      resized_halo.centroid_node_area.push_back(resized_halo.node_size_x[i] *
                                                resized_halo.node_size_y[i]);
    }
    resized_halo.node_centroid_color = {0, 0, 0, -1, -1};
    resized_halo.centroid_color_names = {"bank"};
    CentroidGroup group;
    group.name = "bank";
    group.colors = {0};
    resized_halo.centroid_groups = {group};
    auto centroid_cp = cp;
    centroid_cp.centroid_weight = .25;
    centroid_cp.checkpoint = [](const Design &snapshot, int) {
      GpUpdate u;
      u.node_size_x = snapshot.node_size_x;
      u.node_size_x[0] = 3;
      return u;
    };
    int centroid_areas = 0;
    centroid_cp.observer_interval = 1;
    centroid_cp.observer = [&](const Design &snapshot, PlacementStage, int) {
      if (snapshot.node_size_x[0] == 3) {
        check(snapshot.centroid_node_area[0] == 3);
        ++centroid_areas;
      }
    };
    place(resized_halo, centroid_cp);
    check(centroid_areas > 0 && resized_halo.node_size_x[0] == 3);
    // Read-only center-step telemetry must preserve the ordinary GP trajectory.
    auto plain_motion = fixture(), measured_motion = plain_motion;
    plain_motion.node_macro = measured_motion.node_macro = {1, 0, 0, 0, 0};
    Options motion_options;
    motion_options.max_iterations = motion_options.min_iterations = 12;
    motion_options.bins_x = motion_options.bins_y = 8;
    motion_options.add_fillers = motion_options.legalize =
        motion_options.detailed_place = false;
    place(plain_motion, motion_options);
    std::vector<double> measured, expected;
    int previous_iteration = -1;
    double previous_x = 0, previous_y = 0;
    motion_options.macro_step_observer = [&](int iteration, double sum,
                                             bool valid) {
      check(valid == (iteration != 0));
      measured.push_back(sum);
    };
    motion_options.macro_motion_observer =
        [&](const Design &, const std::vector<double> &x,
            const std::vector<double> &y, int iteration) {
          expected.push_back(
              previous_iteration < 0
                  ? 0
                  : std::hypot(x[0] - previous_x, y[0] - previous_y));
          previous_x = x[0];
          previous_y = y[0];
          previous_iteration = iteration;
        };
    place(measured_motion, motion_options);
    check(measured.size() == 13 && measured == expected);
    check(plain_motion.node_x == measured_motion.node_x &&
          plain_motion.node_y == measured_motion.node_y);
    // Without a GP loop only the overflow test decides convergence.
    auto legal_only = fixture();
    Options legal_options;
    legal_options.max_iterations = 0;
    legal_options.post_density_factor = 4;
    legal_options.bins_x = legal_options.bins_y = 8;
    legal_options.add_fillers = legal_options.detailed_place = false;
    legal_options.stop_overflow = 1;
    check(place(legal_only, legal_options).converged);
    calls.clear();
    cp.checkpoint_interval = 0;
    cp.checkpoint = [&](const Design &, int i) {
      calls.push_back(i);
      return GpUpdate{};
    };
    auto off = fixture();
    place(off, cp);
    check(calls.empty());
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
