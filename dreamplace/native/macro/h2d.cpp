// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/h2d.hpp"
#include "dreamplace/native/macro/check.hpp"
#include "dreamplace/native/macro/fillers.hpp"
#include "dreamplace/native/rudy.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
namespace DreamPlace::native::macro {
RudyModel route_score_model(RouteModel active) {
  return active == RouteModel::rsmt || active == RouteModel::uroute
             ? RudyModel::rsmt
             : RudyModel::bbox;
}
Result place_h2d(Design &design, std::span<const NodeInfo> info,
                 const Parameters &parameters, const Topology &input,
                 const Options &cell_options) {
  design.validate();
  parameters.validate();
  if (info.size() != design.node_x.size())
    throw std::invalid_argument("H2D physical metadata mismatch");
  Parameters p = parameters;
  const auto &l = design.layout;
  if (!p.target_density)
    p.target_density = cell_options.target_density;
  if (!p.route_cap_h)
    p.route_cap_h = 1 / l.site_width;
  if (!p.route_cap_v)
    p.route_cap_v = 1 / l.row_height;
  p.validate();
  Design working = design;
  std::map<std::string, int> source_names;
  for (int i = 0; i < int(info.size()); ++i)
    source_names[info[i].name] = i;
  for (const auto &fixed : input.fixed) {
    auto found = source_names.find(fixed.name);
    if (found == source_names.end() || !info[found->second].macro)
      throw std::invalid_argument("unknown fixed H2D macro");
    int i = found->second;
    if (i >= design.num_movable &&
        (design.node_x[i] != fixed.x || design.node_y[i] != fixed.y ||
         (design.node_orientation.empty()
              ? Orientation::N
              : design.node_orientation[i]) != fixed.orientation))
      throw std::invalid_argument(
          "fixed constraint conflicts with fixed input: " + fixed.name);
    working.node_x[i] = fixed.x;
    working.node_y[i] = fixed.y;
    // Validating the whole design per constraint is O(constraints x design).
    set_orientation_unchecked(working, i, fixed.orientation);
  }
  if (!input.fixed.empty())
    working.validate();

  Topology t = seed_topology(working, info, input.dataflow, p, input);
  Result out;
  out.resolved = p;
  std::vector<int> macro_nodes, audit_nodes;
  std::vector<double> widths;
  for (int i = 0; i < int(info.size()); ++i) {
    if (info[i].macro)
      macro_nodes.push_back(i);
    else if (i < working.num_movable && working.node_size_y[i] > 0)
      widths.push_back(working.node_size_x[i]);
  }
  if (macro_nodes.empty()) {
    auto r = place(working, cell_options);
    out.placement = r;
    out.legal = r.legal;
    out.hpwl = r.hpwl;
    out.topology = t;
    design = std::move(working);
    return out;
  }
  audit_nodes = macro_nodes;
  for (int i = design.num_movable; i < int(info.size()); ++i)
    if (!info[i].macro && design.node_size_x[i] > 0 &&
        design.node_size_y[i] > 0)
      audit_nodes.push_back(i);
  double filler = 0;
  if (!widths.empty()) {
    std::sort(widths.begin(), widths.end());
    size_t low = widths.size() / 20, high = widths.size() - low;
    for (size_t i = low; i < high; ++i)
      filler += widths[i];
    filler /= high - low;
  }
  const double route_scale = (l.xh - l.xl) + (l.yh - l.yl);
  RouteModel active_route = cell_options.route.model == RouteModel::rsmt
                                ? RouteModel::rsmt
                                : RouteModel::rudy;
  const bool latches_route = route_latches(cell_options.route);
  bool route_probed = false;
  DreamPlace::native::Result route_progress;
  auto route_score = [&](const Design &d) {
    if (cell_options.route.model == RouteModel::off)
      return 0.;
    RudyOptions opts{p.bins, p.bins, p.route_cap_h, p.route_cap_v};
    opts.model = route_score_model(active_route);
    opts.degree = cell_options.route.degree;
    opts.accuracy = cell_options.route.accuracy;
    opts.initialize = cell_options.route_initialize;
    return rudy(d, opts).overflow;
  };
  auto merit = [&](const Design &d) {
    RudyResult route;
    route.overflow = route_score(d);
    return std::pair{hpwl(d) + p.alpha * route_scale * route.overflow,
                     route.overflow};
  };
  // Every legal cells-only pass (each round and the coarse-pitch recompaction)
  // uses the same resolved controls and keeps solver frames private.
  const auto final_cells_options = [&](const Abstracted &frozen,
                                       const Topology &topology,
                                       const Design &physical_source) {
    const auto &frozen_design = frozen.design;
    Options cells = cell_options;
    if (cell_options.route.probe &&
        cell_options.route.probe_model == RouteModel::uroute) {
      const auto expand_snapshot = [frozen, topology,
                                    physical_source](const Design &snapshot) {
        auto abstracted = frozen;
        abstracted.design = snapshot;
        auto expanded = physical_source;
        expand(abstracted, topology, expanded);
        return expanded;
      };
      if (cell_options.uroute)
        cells.route_physical_snapshot = expand_snapshot;
      else if (cell_options.route_probe_observer) {
        cells.route_probe_observer =
            [expand_snapshot, observer = cell_options.route_probe_observer](
                const Design &snapshot, int nx, int ny) {
              observer(expand_snapshot(snapshot), nx, ny);
            };
      }
    }
    cells.seed = p.seed;
    cells.num_threads = p.threads;
    cells.target_density = p.target_density;
    cells.filler_width = bounded_filler_width(frozen_design, cells, filler);
    cells.routability = p.routability;
    cells.route_capacity_h = p.route_cap_h;
    cells.route_capacity_v = p.route_cap_v;
    cells.order_weight = 0;
    cells.legalize = true;
    cells.detailed_place = true;
    cells.observer = {};
    cells.route.probe = cell_options.route.probe && !route_probed;
    if (latches_route && (active_route == RouteModel::rsmt ||
                          active_route == RouteModel::uroute))
      cells.route.model = active_route;
    return cells;
  };
  Parameters realization_parameters = p;
  realization_parameters.pitch_x = realization_parameters.pitch_y = 0;
  Design best;
  Topology best_topology;
  double best_merit = std::numeric_limits<double>::infinity();
  SequencePair previous = t.realized;
  GpState state;
  std::string state_shape;
  double last = std::numeric_limits<double>::infinity();
  std::optional<Design> last_design;
  // Auto and a failed probe have one routing-model latch for the entire H2D
  // stage. Every merit comparison must use that model, including incumbents
  // from earlier rounds. A successful probe also runs only once per stage.
  // A GP with routability disabled reports "off"; it must not downgrade an
  // explicit RSMT score or an already latched auto stage.
  const auto observe_route = [&](const DreamPlace::native::Result &result) {
    route_probed |= result.route_probed;
    if (result.route_probed && result.route_probe_overflow >= 0)
      route_progress = result;
    const auto rank = [](RouteModel model) {
      return model == RouteModel::uroute ? 2
             : model == RouteModel::rsmt ? 1
                                         : 0;
    };
    if (!latches_route || rank(result.route_model) <= rank(active_route))
      return;
    active_route = result.route_model;
    if (std::isfinite(best_merit)) {
      const auto [value, route] = merit(best);
      best_merit = value;
      for (auto &old : out.iterations)
        if (old.selected) {
          old.merit = value;
          old.rudy_overflow = route;
        }
    }
    if (last_design)
      last = merit(*last_design).first;
    last_design.reset();
    // Subsequent abstract passes use explicit RSMT controls. Their signature
    // differs from the prior continuation, so start a matching solver state.
    state = {};
    state_shape.clear();
  };
  const int rounds = p.loop ? p.K : 1;
  for (int round = 0; round < rounds; ++round) {
    auto start = std::chrono::steady_clock::now();
    Iteration stat;
    stat.round = round;
    auto checkpoint = start;
    auto phase = [&](const std::string &name) {
      auto now = std::chrono::steady_clock::now();
      stat.phases[name] =
          std::chrono::duration<double>(now - checkpoint).count();
      checkpoint = now;
    };
    auto a = abstract_design(working, info, t, p);
    const size_t physical_nets = a.design.net_weights.size();
    add_pseudonets(a, t, p);
    Options control = cell_options;
    control.max_iterations = p.iters_per_outer;
    control.bins_x = control.bins_y = p.bins;
    control.num_threads = p.threads;
    control.seed = p.seed;
    control.target_density = p.target_density;
    control.legalize = control.detailed_place = false;
    control.order_weight = p.order ? p.lambda0 * std::pow(p.gamma, round) : 0;
    control.order_growth = 1;
    control.schedule_iterations = p.iters_per_outer;
    control.filler_width = bounded_filler_width(a.design, control, filler);
    control.routability = p.routability;
    control.route_capacity_h = p.route_cap_h;
    control.route_capacity_v = p.route_cap_v;
    control.observer = {};
    // The abstract macro geometry contains meta-macros and pseudo pins. The
    // learned feature contract applies only to expanded physical candidates;
    // final_cells_options performs that one-shot probe and bounded feedback.
    control.route.probe = cell_options.route.probe && !route_probed &&
                          cell_options.route.probe_model != RouteModel::uroute;
    if (latches_route && (active_route == RouteModel::rsmt ||
                          active_route == RouteModel::uroute))
      control.route.model = RouteModel::rsmt;
    // A refreshed pseudo-net endpoint can change the abstract connectivity.
    // Reuse full GP history only when that solver topology still matches.
    std::string shape = std::to_string(a.design.node_x.size()) + ":" +
                        std::to_string(a.design.pin2node.size()) + ":" +
                        std::to_string(a.design.num_movable);
    for (int v : a.design.pin2node)
      shape += ":" + std::to_string(v);
    if (shape != state_shape) {
      state = {};
      state_shape = shape;
    }
    auto gp = place(a.design, control, &state);
    observe_route(gp);
    stat.gp_iterations = gp.iterations;
    stat.gp_overflow = gp.overflow;
    stat.gp_converged = gp.converged;
    phase("abstract_pseudonets_gp");
    update_groups(a, t);
    expand(a, t, working);
    if (cell_options.observer)
      cell_options.observer(working, PlacementStage::global_placement,
                            gp.iterations);
    // Pseudo-nets never enter physical compaction or the final cells pass.
    (void)physical_nets;
    auto physical = abstract_design(working, info, t, p);
    Realization realized = group_legalize(physical, t, realization_parameters);
    stat.lp_solves += realized.lp_solves;
    phase("group_legalize");
    if (p.toposearch && p.search != "none" && round < p.freeze_after) {
      auto search = toposearch(working, info, realized, realization_parameters,
                               previous, route_score);
      stat.lp_solves += search.lp_solves;
      stat.evaluations = search.evaluations;
      stat.screened = search.screened;
      stat.accepted_moves = search.accepted;
      if (search.cost <=
          search.seed_cost -
              p.hysteresis * std::max(1., std::abs(search.seed_cost)))
        realized = std::move(search.best);
      previous = std::move(search.sequence);
    }
    phase("toposearch");
    auto compacted_input = abstract_design(working, info, realized.topology, p);
    auto tightened = compact(compacted_input, realized.topology,
                             realization_parameters, realized.graph);
    stat.lp_solves += tightened.lp_solves;
    if (tightened.legal)
      realized = std::move(tightened);
    phase("compact");
    stat.slack = realized.slack;
    stat.relaxed = realized.relaxed;
    Topology candidate = realized.topology;
    candidate.realized_graph = realized.graph;
    auto frozen = abstract_design(working, info, candidate, p, true);
    const Options cells = final_cells_options(frozen, candidate, working);
    DreamPlace::native::Result cell_placement;
    try {
      cell_placement = place(frozen.design, cells);
    } catch (const std::runtime_error &e) {
      if (!std::isfinite(best_merit))
        throw;
      best_topology.diagnostics.push_back(
          std::string("rejected later cells pass: ") + e.what());
      stat.seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start)
                         .count();
      out.iterations.push_back(stat);
      break;
    }
    observe_route(cell_placement);
    phase("cells");
    Design placed = working;
    expand(frozen, candidate, placed);
    stat.legal = is_legal(placed) && check(placed, audit_nodes).legal();
    if (!stat.legal) {
      if (!std::isfinite(best_merit))
        throw std::runtime_error("H2D expanded placement is illegal");
      best_topology.diagnostics.push_back("rejected illegal later placement");
      stat.seconds = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start)
                         .count();
      out.iterations.push_back(stat);
      break;
    }
    auto [value, route] = merit(placed);
    stat.hpwl = hpwl(placed);
    stat.rudy_overflow = route;
    stat.merit = value;
    stat.anchor_deviation = anchor_deviation(candidate, l);
    stat.order_violation = order_violation(candidate);
    if (value < best_merit) {
      for (auto &old : out.iterations)
        old.selected = false;
      stat.selected = true;
      out.placement = cell_placement;
      best_merit = value;
      best = placed;
      best_topology = candidate;
    }
    stat.seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    out.iterations.push_back(stat);
    working = std::move(placed);
    t = std::move(candidate);

    if (round + 1 >= p.freeze_after && std::isfinite(last) &&
        (last - value) <= p.eps * std::max(1., std::abs(last)))
      break;
    last = value;
    if (latches_route && active_route != RouteModel::uroute)
      last_design = working;
  }
  // L22: coarse-pitch compaction on the accepted topology, followed by a fresh
  // cells-only legal pass. Accept only if legal and the shared merit improves.
  if (p.pitch_x > l.site_width || p.pitch_y > l.row_height) {
    try {
      auto a = abstract_design(best, info, best_topology, p);
      auto snapped =
          compact(a, best_topology, p, position_graph(a, best_topology));
      if (snapped.legal) {
        auto frozen = abstract_design(best, info, snapped.topology, p, true);
        const Options cells =
            final_cells_options(frozen, snapped.topology, best);
        const auto snapped_placement = place(frozen.design, cells);
        observe_route(snapped_placement);
        auto candidate = best;
        expand(frozen, snapped.topology, candidate);
        auto [value, route] = merit(candidate);
        (void)route;
        if (is_legal(candidate) && check(candidate, audit_nodes).legal() &&
            value <= best_merit) {
          out.placement = snapped_placement;
          best = std::move(candidate);
          best_topology = std::move(snapped.topology);
          best_merit = value;
        }
      }
    } catch (const std::runtime_error &e) {
      best_topology.diagnostics.push_back(
          std::string("coarse-pitch candidate rejected: ") + e.what());
    }
  }
  out.topology = std::move(best_topology);

  out.topology.macro_fingerprint = fingerprint(design, info);
  out.topology.netlist_hash = fingerprint(design, info, true);
  for (const auto &anchor : out.topology.anchors) {
    auto single = out.topology;
    single.anchors = {anchor};
    double amount = anchor_deviation(single, l);
    if (amount > 1e-6)
      out.violations.push_back({anchor.id, "anchor", amount});
  }
  for (const auto &edge : out.topology.order) {
    auto single = out.topology;
    single.order = {edge};
    double amount = order_violation(single);
    if (amount > 1e-6)
      out.violations.push_back({edge.id, "order", amount});
  }
  out.legal = true;
  out.hpwl = hpwl(best);
  out.rudy_overflow = route_score(best);
  out.merit = best_merit;
  if (latches_route)
    out.placement.route_model = active_route;
  out.placement.route_probed |= route_probed;
  if (route_probed) {
    out.placement.route_probe_overflow = route_progress.route_probe_overflow;
    out.placement.route_uroute_switch_iteration =
        route_progress.route_uroute_switch_iteration;
    out.placement.route_switch_uroute_rsmt =
        route_progress.route_switch_uroute_rsmt;
    out.placement.route_switch_uroute = route_progress.route_switch_uroute;
  }
  design = std::move(best);
  return out;
}
} // namespace DreamPlace::native::macro
