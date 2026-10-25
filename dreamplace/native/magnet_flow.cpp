// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/magnet.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace DreamPlace::native::detail {
namespace {
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
// A callback may resize ordinary cells, but topology resynthesis belongs to
// the earlier FP solve. Macro-stage STA retains the original pin identities.
template <class T>
void update_array(std::vector<T> &to, const std::vector<T> &from) {
  if (from.empty())
    return;
  require(from.size() == to.size(),
          "macro magnet checkpoint array size changed");
  to = from;
}
void add_costs(Result &to, const Result &from) {
  to.global_placement_seconds += from.global_placement_seconds;
  to.poisson_seconds += from.poisson_seconds;
  to.rudy_seconds += from.rudy_seconds;
  to.rudy_calls += from.rudy_calls;
  to.macro_routing_seconds += from.macro_routing_seconds;
  to.macro_alignment_seconds += from.macro_alignment_seconds;
  to.local_density_seconds += from.local_density_seconds;
  to.local_density_grid_seconds += from.local_density_grid_seconds;
  to.macro_cell_density_seconds += from.macro_cell_density_seconds;
  to.macro_density_pull_seconds += from.macro_density_pull_seconds;
  to.detailed_placement_seconds += from.detailed_placement_seconds;
  to.history.insert(to.history.end(), from.history.begin(), from.history.end());
  to.checkpoint_history.insert(to.checkpoint_history.end(),
                               from.checkpoint_history.begin(),
                               from.checkpoint_history.end());
}
} // namespace

MagnetRuntime::MagnetRuntime(const Options &options)
    : options_(options), stability_(options.magnet) {}

void MagnetRuntime::remap_pulls(const std::vector<int> &previous_node) {
  if (macro_pulls.empty())
    return;
  int old_count = 0;
  for (const auto &pull : macro_pulls)
    old_count = std::max(old_count, pull.node + 1);
  std::vector<int> next_node(old_count, -1);
  for (size_t i = 0; i < previous_node.size(); ++i)
    if (previous_node[i] >= 0 && previous_node[i] < old_count)
      next_node[previous_node[i]] = int(i);
  std::vector<MacroPull> kept;
  kept.reserve(macro_pulls.size());
  for (auto pull : macro_pulls)
    if (next_node[pull.node] >= 0) {
      pull.node = next_node[pull.node];
      kept.push_back(pull);
    }
  macro_pulls = std::move(kept);
}

void MagnetRuntime::event(int iteration, MagnetEvent kind) {
  stability_.event(iteration, kind);
  report->events.push_back({iteration, kind});
  previous_hpwl_ = -1;
  stable_hpwl_ = 0;
  if (kind == MagnetEvent::round || kind == MagnetEvent::snap) {
    // Rebase the next solver at the changed geometry. Packing/snap motion is
    // recorded by its frames, never attributed to a GP optimizer step.
    motion_iteration_ = -1;
    motion_nodes_.clear();
    motion_weights_.clear();
    motion_centers_.clear();
  }
}
void MagnetRuntime::observe(const Design &d, const std::vector<double> &x,
                            const std::vector<double> &y, int iteration) {
  const auto movement = stability_.observe(d, x, y, iteration);
  if (motion_nodes_.empty()) {
    for (int i = 0; i < d.num_movable; ++i)
      if (!d.node_macro.empty() && d.node_macro[i]) {
        motion_nodes_.push_back(i);
        motion_weights_.push_back(1);
      }
    for (const auto &group : report->mapping.groups) {
      const auto found = std::lower_bound(motion_nodes_.begin(),
                                          motion_nodes_.end(), group.node);
      require(found != motion_nodes_.end() && *found == group.node,
              "macro motion group has no movable node");
      motion_weights_[found - motion_nodes_.begin()] =
          int(group.members.size());
    }
  }
  MagnetReport::MotionSample sample;
  sample.iteration = iteration;
  sample.movement = movement;
  sample.gp_macros = int(motion_nodes_.size());
  sample.center_step_valid = motion_iteration_ + 1 == iteration &&
                             motion_centers_.size() == motion_nodes_.size();
  std::vector<std::pair<double, double>> centers;
  centers.reserve(motion_nodes_.size());
  for (size_t k = 0; k < motion_nodes_.size(); ++k) {
    const int i = motion_nodes_[k];
    const double cx = x[i] + d.node_size_x[i] / 2,
                 cy = y[i] + d.node_size_y[i] / 2;
    centers.emplace_back(cx, cy);
    sample.physical_macros += motion_weights_[k];
    if (sample.center_step_valid)
      sample.center_step_sum +=
          motion_weights_[k] * std::hypot(cx - motion_centers_[k].first,
                                          cy - motion_centers_[k].second);
  }
  motion_iteration_ = iteration;
  motion_centers_ = std::move(centers);
  report->motion_history.push_back(sample);
  if (options_.macro_step_observer)
    options_.macro_step_observer(iteration, sample.center_step_sum,
                                 sample.center_step_valid);
}
void MagnetRuntime::sample(const Design &d, GpSample &sample) {
  const int i = sample.iteration;
  report->overflow_history.emplace_back(i, sample.overflow);
  const auto movement = stability_.movement();
  if (!started_ && movement.window_complete &&
      movement.maximum <= options_.magnet.start_displacement) {
    started_ = true;
    report->start_iteration = i;
    report->start_observed_displacement = movement.maximum;
    // Readiness is observational: retain the existing per-object history.
    report->events.push_back({i, MagnetEvent::ready});
  }
  sample.magnet_round = round;
  sample.magnet_groups = groups;
  sample.macro_move_max = movement.maximum;
  sample.macro_move_rms = movement.rms;
  sample.last_instability_iteration = stability_.last_event();
  const bool hpwl_close =
      previous_hpwl_ >= 0 && std::abs(sample.hpwl - previous_hpwl_) <=
                                 options_.convergence_tolerance *
                                     std::max({1e-12, std::abs(sample.hpwl),
                                               std::abs(previous_hpwl_)});
  stable_hpwl_ = hpwl_close ? stable_hpwl_ + 1 : 0;
  previous_hpwl_ = sample.hpwl;
  if (!started_ || finalizing ||
      (options_.magnet.max_rounds && round >= options_.magnet.max_rounds) ||
      !movement.window_complete)
    return;
  const auto start = std::chrono::steady_clock::now();
  const int bx =
      options_.auxiliary_bins ? options_.auxiliary_bins : options_.bins_x;
  const int by =
      options_.auxiliary_bins ? options_.auxiliary_bins : options_.bins_y;
  auto search = magnet_candidates(d, options_.magnet, MagnetDensity(d, bx, by),
                                  &stability_, candidate_check);
  prepared.reset();
  std::string reason = "no_eligible_pair";
  if (!search.merges.empty()) {
    PreparedMagnetMerge transaction{
        prepare_snapshot ? prepare_snapshot(d) : d, report->mapping, {}};
    transaction.result =
        merge_macros(transaction.design, transaction.mapping, search.merges);
    if (transaction.result.applied) {
      prepared = std::move(transaction);
      reason = "pending";
    } else {
      reason = transaction.result.reason;
      for (auto &candidate : search.candidates)
        if (candidate.selected) {
          candidate.selected = false;
          candidate.reason = reason;
        }
      search.merges.clear();
    }
  }
  pending = search.merges;
  // Cells can empty a gap without moving its macros. Recheck at every stable
  // GP sample instead of caching a rejection solely by macro geometry.
  report->rounds.push_back(
      {i, sample.overflow, movement, std::move(search), {}, std::move(reason)});
  report->seconds +=
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
}
bool MagnetRuntime::ready(int iteration) const {
  if (!pending.empty() || !stability_.settled(iteration) || stable_hpwl_ < 2)
    return false;
  return true;
}

Result magnet_place(Design &design, const Options &options) {
  options.magnet.validate();
  require(design.node_pad_x.empty() && design.node_pad_y.empty(),
          "macro magnet requires physical unpadded input");
  require(!design.node_macro.empty() &&
              std::any_of(design.node_macro.begin(),
                          design.node_macro.begin() + design.num_movable,
                          [](auto m) { return m != 0; }),
          "macro magnet requires movable macros");
  auto runtime = std::make_shared<MagnetRuntime>(options);
  runtime->report->layout = design.layout;
  runtime->report->target_density = options.target_density;
  runtime->report->options = options.magnet;
  Design source = design, grouped = design;
  auto &mapping = runtime->report->mapping;
  mapping.source_to_node.resize(source.node_x.size());
  std::iota(mapping.source_to_node.begin(), mapping.source_to_node.end(), 0);
  runtime->prepare_snapshot = [&](const Design &snapshot) {
    return unpad_snapshot(snapshot, grouped);
  };
  // Only macro identities are indexed here. Scanning all standard cells for
  // every candidate would dominate reporting on designs with many memories.
  std::vector<std::vector<RigidMember>> node_members;
  const auto index_members = [&] {
    node_members.assign(grouped.num_movable, {});
    for (int i = 0; i < source.num_movable; ++i)
      if (source.node_macro[i])
        node_members[mapping.source_to_node[i]].push_back({i, 0, 0});
    for (const auto &group : mapping.groups)
      node_members[group.node] = group.members;
  };
  if (options.macro_merge_check) {
    runtime->candidate_check = [&](const MagnetCandidate &candidate) {
      std::vector<MagnetMacroFrame> members;
      const auto &request = candidate.merge;
      const double origin_x = request.poses.front().first;
      const double origin_y = request.poses.front().second;
      for (size_t k = 0; k < request.members.size(); ++k)
        for (const auto &member : node_members[request.members[k]]) {
          const int i = member.source;
          members.push_back({i, request.poses[k].first - origin_x + member.x,
                             request.poses[k].second - origin_y + member.y,
                             source.node_size_x[i], source.node_size_y[i],
                             source.node_orientation.empty()
                                 ? Orientation::N
                                 : source.node_orientation[i],
                             0, false});
        }
      return options.macro_merge_check(members);
    };
  }
  Options control = options;
  control.magnet.enabled = false;
  if (options.uroute && (options.route.model == RouteModel::uroute ||
                         options.route.uroute_ladder ||
                         (options.route.probe &&
                          options.route.probe_model == RouteModel::uroute))) {
    control.route_physical_snapshot = [&](const Design &snapshot) {
      const auto unpadded = unpad_snapshot(snapshot, grouped);
      auto expanded = expand_macros(unpadded, source, mapping);
      return options.route_physical_snapshot
                 ? options.route_physical_snapshot(expanded)
                 : expanded;
    };
  }
  // The shared rigid runtime already delivers physical-member telemetry.
  control.macro_step_observer = nullptr;
  if (!options.magnet.before_snap)
    control.magnet_runtime = runtime;
  else
    control.macro_motion_observer =
        [&](const Design &snapshot, const std::vector<double> &x,
            const std::vector<double> &y,
            int iteration) { runtime->observe(snapshot, x, y, iteration); };
  control.automatic_bins = control.dynamic_bins = false;
  control.schedule_iterations = options.schedule_iterations
                                    ? options.schedule_iterations
                                    : options.max_iterations;
  control.legalize = control.detailed_place = control.macro_only_legalize =
      false;
  control.before_legalize = nullptr;
  control.after_macro_legalize = nullptr;
  if (options.checkpoint) {
    control.checkpoint = [&](const Design &snapshot, int iteration) {
      auto expanded = expand_macros(snapshot, source, mapping);
      auto update = options.checkpoint(expanded, iteration);
      require(!update.replacement,
              "macro magnet does not support macro-stage topology replacement");
      update_array(expanded.net_weights, update.net_weights);
      update_array(expanded.node_size_x, update.node_size_x);
      update_array(expanded.node_size_y, update.node_size_y);
      update_array(expanded.pin_offset_x, update.pin_offset_x);
      update_array(expanded.pin_offset_y, update.pin_offset_y);
      const auto contracted =
          contract_macros(expanded, snapshot, source, mapping);
      if (!update.node_size_x.empty())
        update.node_size_x = contracted.node_size_x;
      if (!update.node_size_y.empty())
        update.node_size_y = contracted.node_size_y;
      if (!update.pin_offset_x.empty())
        update.pin_offset_x = contracted.pin_offset_x;
      if (!update.pin_offset_y.empty())
        update.pin_offset_y = contracted.pin_offset_y;
      source = std::move(expanded);
      if (options.magnet.before_snap)
        runtime->event(iteration, MagnetEvent::timing);
      return update;
    };
  }
  if (options.terminal_update) {
    control.terminal_update = [&](const Design &snapshot, int iteration) {
      const auto expanded = expand_macros(snapshot, source, mapping);
      auto updates = options.terminal_update(expanded, iteration);
      for (auto &terminal : updates) {
        require(terminal.node >= source.num_movable &&
                    terminal.node < int(source.node_x.size()) &&
                    source.node_size_x[terminal.node] == 0 &&
                    source.node_size_y[terminal.node] == 0,
                "invalid macro magnet terminal update");
        terminal.node = mapping.source_to_node[terminal.node];
      }
      if (options.magnet.before_snap && !updates.empty())
        runtime->event(iteration, MagnetEvent::terminal);
      return updates;
    };
  }
  if (options.magnet.before_snap && options.macro_pull_update) {
    control.macro_pull_update = [&](const Design &snapshot, int iteration) {
      auto update = options.macro_pull_update(snapshot, iteration);
      if (update)
        runtime->event(iteration, MagnetEvent::pull);
      return update;
    };
  }
  const auto capture = [&](const Design &expanded, int iteration,
                           const char *event) {
    const auto start = std::chrono::steady_clock::now();
    MagnetFrame frame;
    frame.iteration = iteration;
    frame.event = event;
    for (size_t i = 0; i < expanded.node_x.size(); ++i) {
      if (expanded.node_macro.empty() || !expanded.node_macro[i])
        continue;
      frame.macros.push_back(
          {int(i), expanded.node_x[i], expanded.node_y[i],
           expanded.node_size_x[i], expanded.node_size_y[i],
           expanded.node_orientation.empty() ? Orientation::N
                                             : expanded.node_orientation[i],
           expanded.node_group_id.empty() ? 0 : expanded.node_group_id[i],
           i >= size_t(expanded.num_movable)});
    }
    const MagnetDensity density(expanded, std::min(64, options.bins_x),
                                std::min(64, options.bins_y));
    frame.bins_x = density.bins_x();
    frame.bins_y = density.bins_y();
    frame.density = density.grid();
    runtime->report->frames.push_back(std::move(frame));
    runtime->report->seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
  };
  control.observer_interval = std::min(20, options.observer_interval);
  control.observer = [&](const Design &snapshot, PlacementStage stage,
                         int iteration) {
    if (stage == PlacementStage::finished ||
        (stage == PlacementStage::initial && iteration))
      return;
    const auto expanded = expand_macros(snapshot, source, mapping);
    capture(expanded, iteration, stage_name(stage));
    if (options.observer && (stage != PlacementStage::global_placement ||
                             iteration % options.observer_interval == 0))
      options.observer(expanded, stage, iteration);
  };
  Result costs, last;
  int iteration = 0;
  double initial_hpwl = hpwl(source);
  // Each successful batch contracts actual Design nodes. Geometry-dependent
  // scratch, density/filler maps and momentum restart at current poses. Stage
  // decisions and untouched pull targets survive in the shared runtime.
  if (options.magnet.before_snap) {
    // This mode has exactly the ordinary GP lifetime and force state. The
    // frozen final geometry is searched without a motion window; no GP solve
    // follows these merges. Geometry/density/host feasibility still apply.
    index_members();
    last = place(grouped, control);
    add_costs(costs, last);
    iteration = last.iterations;
    runtime->report->start_iteration = iteration;
    runtime->report->events.push_back({iteration, MagnetEvent::ready});
    while (!options.magnet.max_rounds ||
           runtime->round < options.magnet.max_rounds) {
      index_members();
      const auto begin = std::chrono::steady_clock::now();
      const int bx =
          options.auxiliary_bins ? options.auxiliary_bins : options.bins_x;
      const int by =
          options.auxiliary_bins ? options.auxiliary_bins : options.bins_y;
      auto search = magnet_candidates(grouped, options.magnet,
                                      MagnetDensity(grouped, bx, by), nullptr,
                                      runtime->candidate_check);
      for (auto &candidate : search.candidates) {
        for (const auto &m : node_members[candidate.a])
          candidate.source_a.push_back(m.source);
        for (const auto &m : node_members[candidate.b])
          candidate.source_b.push_back(m.source);
      }
      auto proposed = grouped;
      auto proposed_mapping = mapping;
      MacroMergeResult merged;
      if (!search.merges.empty())
        merged = merge_macros(proposed, proposed_mapping, search.merges);
      const auto reason = search.merges.empty() ? "no_eligible_pair"
                          : merged.applied      ? "merged"
                                                : merged.reason;
      if (!merged.applied)
        for (auto &candidate : search.candidates)
          if (candidate.selected) {
            candidate.selected = false;
            candidate.reason = reason;
          }
      runtime->report->rounds.push_back({iteration,
                                         last.overflow,
                                         {},
                                         std::move(search),
                                         merged.group_ids,
                                         reason});
      runtime->report->seconds += std::chrono::duration<double>(
                                      std::chrono::steady_clock::now() - begin)
                                      .count();
      if (!merged.applied)
        break;
      capture(expand_macros(grouped, source, mapping), iteration,
              "before_magnet");
      grouped = std::move(proposed);
      mapping = std::move(proposed_mapping);
      ++runtime->round;
      runtime->groups = int(mapping.groups.size());
      runtime->event(iteration, MagnetEvent::round);
      capture(expand_macros(grouped, source, mapping), iteration, "magnet");
    }
  } else
    while (iteration < options.max_iterations) {
      index_members();
      control.max_iterations = options.max_iterations - iteration;
      control.iteration_offset = iteration;
      last = place(grouped, control);
      add_costs(costs, last);
      iteration += last.iterations;
      for (auto &attempt : runtime->report->rounds) {
        for (auto &candidate : attempt.search.candidates) {
          if (!candidate.source_a.empty())
            continue;
          for (const auto &m : node_members[candidate.a])
            candidate.source_a.push_back(m.source);
          for (const auto &m : node_members[candidate.b])
            candidate.source_b.push_back(m.source);
        }
      }
      if (runtime->pending.empty())
        break;
      capture(expand_macros(grouped, source, mapping), iteration,
              "before_magnet");
      const auto start = std::chrono::steady_clock::now();
      require(runtime->prepared.has_value(),
              "macro merge paused without a prepared transaction");
      auto transaction = std::move(*runtime->prepared);
      runtime->prepared.reset();
      runtime->remap_pulls(transaction.result.previous_node);
      grouped = std::move(transaction.design);
      mapping = std::move(transaction.mapping);
      runtime->pending.clear();
      auto &round = runtime->report->rounds.back();
      round.reason = "merged";
      round.group_ids = std::move(transaction.result.group_ids);
      runtime->report->seconds += std::chrono::duration<double>(
                                      std::chrono::steady_clock::now() - start)
                                      .count();
      ++runtime->round;
      runtime->groups = int(mapping.groups.size());
      runtime->event(iteration, MagnetEvent::round);
      capture(expand_macros(grouped, source, mapping), iteration, "magnet");
      if (!last.history.empty())
        control.initial_density_weight = last.history.back().density_weight;
    }
  // Snap and LP operate on the contracted rectangles, never on members. The
  // host can refine their common translation, but cannot break rigid offsets.
  auto wrap_host = [&](const std::function<void(Design &)> &callback)
      -> std::function<void(Design &)> {
    if (!callback)
      return {};
    return [&, callback](Design &snapshot) {
      auto expanded = expand_macros(snapshot, source, mapping);
      callback(expanded);
      snapshot = contract_macros(expanded, snapshot, source, mapping);
      source = std::move(expanded);
    };
  };
  Options finish = control;
  finish.max_iterations = 0;
  finish.iteration_offset = iteration;
  finish.magnet_runtime.reset();
  finish.legalize = true;
  finish.macro_only_legalize = true;
  finish.detailed_place = false;
  finish.before_legalize = wrap_host(options.before_legalize);
  finish.after_macro_legalize = wrap_host(options.after_macro_legalize);
  auto legalized = place(grouped, finish);
  const auto snap_alignment = legalized.alignment;
  add_costs(costs, legalized);
  runtime->event(iteration, MagnetEvent::snap);
  capture(expand_macros(grouped, source, mapping), iteration,
          options.macro_snap ? "snap" : "snap_disabled");
  runtime->finalizing = true;
  // Freeze the snapped group poses during the last cell relaxation. A second
  // legalization cannot split a group or introduce an unobserved macro move.
  const auto original_boxes = grouped.node_boxes;
  if (grouped.node_boxes.empty()) {
    for (size_t i = 0; i < grouped.node_x.size(); ++i)
      grouped.node_boxes.push_back(
          {grouped.layout.xl, grouped.layout.yl,
           grouped.layout.xh - grouped.node_size_x[i],
           grouped.layout.yh - grouped.node_size_y[i]});
  }
  for (int i = 0; i < grouped.num_movable; ++i)
    if (grouped.node_macro[i])
      grouped.node_boxes[i] = {grouped.node_x[i], grouped.node_y[i],
                               grouped.node_x[i], grouped.node_y[i]};
  if (!options.magnet.before_snap && iteration < options.max_iterations) {
    // The selected pre-LP timing snapshot was already delivered. No additional
    // STA/anchor/IO events are introduced merely to settle the frozen macros.
    control.checkpoint = nullptr;
    control.macro_pull_update = nullptr;
    control.terminal_update = nullptr;
    control.max_iterations = options.max_iterations - iteration;
    control.iteration_offset = iteration;
    if (!last.history.empty())
      control.initial_density_weight = last.history.back().density_weight;
    last = place(grouped, control);
    add_costs(costs, last);
    iteration += last.iterations;
  }
  finish.before_legalize = nullptr;
  finish.iteration_offset = iteration;
  finish.macro_snap = false;
  finish.macro_only_legalize = options.macro_only_legalize;
  finish.detailed_place = options.detailed_place;
  const auto final_legalization = place(grouped, finish);
  add_costs(costs, final_legalization);
  auto final = last;
  final.detailed_history = final_legalization.detailed_history;
  grouped.node_boxes = original_boxes;
  auto expanded = expand_macros(grouped, source, mapping);
  runtime->report->stopped_settled = runtime->settled(iteration);
  final.magnet = runtime->report;
  final.magnet_seconds = runtime->report->seconds;
  final.alignment = snap_alignment;
  final.alignment.macro_legalization_seconds +=
      final_legalization.alignment.macro_legalization_seconds;
  final.alignment.cell_legalization_seconds +=
      final_legalization.alignment.cell_legalization_seconds;
  final.iterations = final.selected_iteration = iteration;
  final.converged = last.converged && runtime->report->stopped_settled;
  final.initial_hpwl = initial_hpwl;
  final.hpwl = hpwl(expanded);
  final.legal = is_legal(expanded);
  final.overflow = last.overflow;
  final.global_placement_seconds = costs.global_placement_seconds;
  final.poisson_seconds = costs.poisson_seconds;
  final.rudy_seconds = costs.rudy_seconds;
  final.rudy_calls = costs.rudy_calls;
  final.macro_routing_seconds = costs.macro_routing_seconds;
  final.macro_alignment_seconds = costs.macro_alignment_seconds;
  final.local_density_seconds = costs.local_density_seconds;
  final.local_density_grid_seconds = costs.local_density_grid_seconds;
  final.macro_cell_density_seconds = costs.macro_cell_density_seconds;
  final.macro_density_pull_seconds = costs.macro_density_pull_seconds;
  final.detailed_placement_seconds = costs.detailed_placement_seconds;
  final.history = std::move(costs.history);
  final.checkpoint_history = std::move(costs.checkpoint_history);
  capture(expanded, iteration, "finished");
  final.magnet_seconds = runtime->report->seconds;
  if (options.observer)
    options.observer(expanded, PlacementStage::finished, iteration);
  design = std::move(expanded);
  return final;
}
} // namespace DreamPlace::native::detail
