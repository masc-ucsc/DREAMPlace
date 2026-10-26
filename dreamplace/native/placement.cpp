// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <future>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

#include "dreamplace/native/local_density.hpp"
#include "dreamplace/native/magnet.hpp"
#include "dreamplace/native/rudy.hpp"
#include "dreamplace/native/spectral.hpp"
#include "dreamplace/native/uroute_weights.hpp"

namespace DreamPlace::native {
bool macro_transparent(const MacroTransparency &o, int completed_steps) {
  if (o.always) {
    return true;
  }
  const auto end =
      std::int64_t(o.warmup) + 2 * std::int64_t(o.period) * o.cycles;
  return o.cycles > 0 && o.period > 0 && completed_steps >= o.warmup &&
         completed_steps < end &&
         ((completed_steps - o.warmup) / o.period) % 2 == 0;
}
double vanished_memory_scale(int iteration, double initial_scale,
                             int growth_iterations, int start_iteration) {
  return initial_scale + (1 - initial_scale) *
                             std::clamp((double(iteration) - start_iteration) /
                                            growth_iterations,
                                        0., 1.);
}
double vanished_memory_local_weight(double initial, int iteration, double floor,
                                    int start_iteration) {
  return std::max(
      floor, initial - .05 * std::max(0., double(iteration) - start_iteration));
}
int select_gp_trial(const std::vector<GpTrialScore> &scores,
                    const GpSearch &policy, std::string *reason) {
  if (scores.empty()) {
    throw std::invalid_argument("empty GP candidate set");
  }
  for (double t : {policy.route_tolerance, policy.timing_tolerance,
                   policy.hpwl_tolerance}) {
    if (!std::isfinite(t) || t < 0 || t > .1) {
      throw std::invalid_argument("invalid GP trial comparison tolerance");
    }
  }
  std::vector<int> active(scores.size());
  std::iota(active.begin(), active.end(), 0);
  for (const auto &s : scores) {
    for (double v : {s.route_overflow, s.wns, s.hpwl, s.tns}) {
      if (!std::isfinite(v)) {
        throw std::invalid_argument("nonfinite GP trial score");
      }
    }
  }
  const auto filter = [&](auto cost, double tolerance, const char *name) {
    if (active.size() == 1) {
      return;
    }
    double best = std::numeric_limits<double>::infinity();
    for (int i : active) {
      best = std::min(best, cost(scores[i]));
    }
    const double band = tolerance * std::max(std::abs(best), 1e-12);
    std::erase_if(active, [&](int i) { return cost(scores[i]) > best + band; });
    if (reason && active.size() == 1) {
      *reason = name;
    }
  };
  if (reason) {
    *reason = "tie_keep_first";
  }
  filter([](const auto &s) { return s.route_overflow; }, policy.route_tolerance,
         "rsmt_rudy");
  filter([](const auto &s) { return -s.wns; }, policy.timing_tolerance, "wns");
  filter([](const auto &s) { return s.hpwl; }, policy.hpwl_tolerance, "hpwl");
  filter([](const auto &s) { return -s.tns; }, policy.timing_tolerance, "tns");
  return active.front();
}
namespace {
void require(bool condition, const char *message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}
bool finite(double v) { return std::isfinite(v); }
bool aligned(double v) { return std::abs(v - std::round(v)) <= 1e-6; }
int count(size_t n) {
  require(n <= static_cast<size_t>(std::numeric_limits<int>::max()),
          "array exceeds 32-bit kernel indices");
  return static_cast<int>(n);
}
double placement_preconditioner(const Design &d, int i, double density_weight) {
  const double degree =
      i < int(d.flat_node2pin_start.size()) - 1
          ? d.flat_node2pin_start[i + 1] - d.flat_node2pin_start[i]
          : 0;
  return std::max(1., degree +
                          density_weight * d.node_size_x[i] * d.node_size_y[i]);
}
void csr_check(const std::vector<int> &flat, const std::vector<int> &start,
               int objects, int pins, const std::vector<int> &inverse) {
  require(start.size() == static_cast<size_t>(objects) + 1 &&
              start.front() == 0 && start.back() == pins,
          "invalid CSR start/end");
  require(flat.size() == static_cast<size_t>(pins),
          "CSR must contain every pin exactly once");
  std::vector<bool> seen(pins, false);
  for (int n = 0; n < objects; ++n) {
    require(start[n] <= start[n + 1] && start[n] >= 0 && start[n + 1] <= pins,
            "nonmonotonic CSR");
    for (int j = start[n]; j < start[n + 1]; ++j) {
      int p = flat[j];
      require(p >= 0 && p < pins && !seen[p], "invalid or duplicate CSR pin");
      seen[p] = true;
      require(inverse[p] == n, "CSR disagrees with inverse map");
    }
  }
}

// Exact physical overflow, excluding fillers; fixed area reduces bin capacity.
double overflow(const Design &d, int bx, int by, double target_density,
                bool scale_macro = false, bool ignore_fixed_macros = false,
                double fixed_macro_scale = 1,
                bool include_reservations = true) {
  const auto &l = d.layout;
  double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  std::vector<double> movable(bx * by), fixed(bx * by);
  double total = 0;
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    if (ignore_fixed_macros && i >= d.num_movable && !d.node_macro.empty() &&
        d.node_macro[i]) {
      continue;
    }
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
           h = d.node_size_y[i];
    if (fixed_macro_scale != 1 && i >= d.num_movable && !d.node_macro.empty() &&
        d.node_macro[i]) {
      x += w * (1 - fixed_macro_scale) / 2;
      y += h * (1 - fixed_macro_scale) / 2;
      w *= fixed_macro_scale;
      h *= fixed_macro_scale;
    }
    const double charge = scale_macro && i < d.num_movable &&
                                  !d.node_macro.empty() && d.node_macro[i]
                              ? target_density
                              : 1;
    if (i < d.num_movable) {
      total += charge * w * h;
    }
    // Clamp before the int conversion: fixed nodes may lie far outside the
    // core.
    const auto cell = [](double v, int count) {
      return static_cast<int>(std::clamp(v, 0., double(count)));
    };
    int x0 = cell(std::floor((x - l.xl) / dx), bx);
    int x1 = cell(std::ceil((x + w - l.xl) / dx), bx);
    int y0 = cell(std::floor((y - l.yl) / dy), by);
    int y1 = cell(std::ceil((y + h - l.yl) / dy), by);
    auto &map = i < d.num_movable ? movable : fixed;
    for (int k = x0; k < x1; ++k) {
      for (int j = y0; j < y1; ++j) {
        double a = std::max(0., std::min(x + w, l.xl + (k + 1) * dx) -
                                    std::max(x, l.xl + k * dx));
        double b = std::max(0., std::min(y + h, l.yl + (j + 1) * dy) -
                                    std::max(y, l.yl + j * dy));
        map[k * by + j] += charge * a * b;
      }
    }
  }
  if (include_reservations) {
    for (const auto &r : d.density_reservations) {
      if (r.materialized_node >= 0)
        continue;
      total += r.width * r.height;
      const int x0 = std::clamp(int(std::floor((r.x - l.xl) / dx)), 0, bx);
      const int x1 =
          std::clamp(int(std::ceil((r.x + r.width - l.xl) / dx)), 0, bx);
      const int y0 = std::clamp(int(std::floor((r.y - l.yl) / dy)), 0, by);
      const int y1 =
          std::clamp(int(std::ceil((r.y + r.height - l.yl) / dy)), 0, by);
      for (int x = x0; x < x1; ++x)
        for (int y = y0; y < y1; ++y) {
          const double w =
              std::max(0., std::min(r.x + r.width, l.xl + (x + 1) * dx) -
                               std::max(r.x, l.xl + x * dx));
          const double h =
              std::max(0., std::min(r.y + r.height, l.yl + (y + 1) * dy) -
                               std::max(r.y, l.yl + y * dy));
          movable[x * by + y] += w * h;
        }
    }
  }
  double excess = 0;
  for (int i = 0; i < bx * by; ++i) {
    excess += std::max(0., movable[i] - target_density *
                                            std::max(0., dx * dy - fixed[i]));
  }
  return total ? excess / total : 0;
}

double memory_settled_iteration(const Options &o) {
  return o.vanished_memory_start + 20. +
         std::max(double(o.vanished_memory_growth_iterations),
                  std::ceil(
                      (o.local_density_weight - o.vanished_memory_local_floor) /
                      .05));
}

void options_check(const Options &o) {
  const auto &search = o.gp_search;
  require(search.start >= 220 && search.interval > search.steps &&
              search.steps >= 1 && search.steps <= 3 && search.rounds >= 0 &&
              search.rounds <= 100 && finite(search.filler_step) &&
              search.filler_step > 0 && search.filler_step <= 1 &&
              finite(search.local_step) && search.local_step > 0 &&
              search.local_step <= 1,
          "invalid GP A/B search controls");
  for (double t : {search.route_tolerance, search.timing_tolerance,
                   search.hpwl_tolerance}) {
    require(finite(t) && t >= 0 && t <= .1,
            "invalid GP A/B comparison tolerance");
  }
  if (o.trial_tuning) {
    require(!search.enabled && finite(o.trial_tuning->filler_charge) &&
                o.trial_tuning->filler_charge >= 0 &&
                o.trial_tuning->filler_charge <= 1 &&
                finite(o.trial_tuning->local_weight) &&
                o.trial_tuning->local_weight >= 0 &&
                o.trial_tuning->local_weight <= 1,
            "invalid internal GP trial tuning");
  }
  const auto &cycle = o.macro_transparency;
  require(cycle.cycles >= 0 && cycle.warmup >= 0 && cycle.period > 0,
          "invalid macro transparency schedule");
  require(!cycle.always || cycle.cycles == 0,
          "always-transparent mode cannot be combined with cycling");
  if (cycle.cycles > 0) {
    require(std::int64_t(cycle.warmup) +
                    (2 * std::int64_t(cycle.cycles) + 1) * cycle.period <=
                o.max_iterations,
            "macro transparency requires a final fully opaque settling period");
  }
  require(finite(o.vanished_memory_local_floor) &&
              o.vanished_memory_local_floor >= 0 &&
              o.vanished_memory_local_floor <= o.local_density_weight &&
              (o.vanished_memory || o.vanished_memory_local_floor == 0),
          "vanished_memory local floor requires the schedule and a value from "
          "zero to initial local weight");
  require(
      finite(o.vanished_memory_initial_scale) &&
          o.vanished_memory_initial_scale >= 0 &&
          o.vanished_memory_initial_scale <= 1 &&
          o.vanished_memory_growth_iterations > 0 &&
          o.vanished_memory_growth_iterations <= 1000000 &&
          o.vanished_memory_start >= 0 && o.vanished_memory_start <= 1000000 &&
          (o.vanished_memory || (o.vanished_memory_initial_scale == .6 &&
                                 o.vanished_memory_growth_iterations == 25 &&
                                 o.vanished_memory_start == 70)),
      "invalid vanished_memory initial scale, growth duration or start");
  if (o.vanished_memory) {
    require(!search.enabled || search.start >= memory_settled_iteration(o),
            "GP A/B search must start after memory growth and local-pressure "
            "settling");
    require(!cycle.enabled(),
            "vanished_memory cannot be combined with transparency cycling");
    require((o.trial_tuning ? o.schedule_iterations : o.max_iterations) >=
                memory_settled_iteration(o),
            "vanished_memory requires settling after the local-pressure ramp");
  }
  if (cycle.enabled() || o.vanished_memory) {
    require(
        o.macro_cell_density_weight == 0 && o.macro_density_pull_weight == 0,
        "macro transparency requires movable-macro density forces disabled");
  }
  require(o.max_iterations >= 0 && o.bins_x >= 2 && o.bins_x <= 1024 &&
              o.bins_y >= 2 && o.bins_y <= 1024,
          "invalid iteration limit or bin count (2..1024)");
  require(finite(o.local_density_weight) && o.local_density_weight >= 0 &&
              finite(o.local_density_range) && o.local_density_range > 0 &&
              finite(o.local_density_jitter) && o.local_density_jitter >= 0 &&
              o.local_density_jitter <= .1 && o.local_density_bins >= 8 &&
              o.local_density_bins <= 1024 &&
              std::has_single_bit(unsigned(o.local_density_bins)) &&
              (o.auxiliary_bins == 0 ||
               (o.auxiliary_bins >= 2 && o.auxiliary_bins <= 1024)),
          "invalid local density controls");
  require(finite(o.local_density_target) &&
              (o.local_density_target == -1 ||
               (o.local_density_target >= 0 && o.local_density_target <= 1)),
          "local density target must be -1 (inherit) or in [0,1]");
  require(finite(o.timing_spread_weight) && o.timing_spread_weight >= 0 &&
              o.timing_spread_weight <= 1,
          "timing spreading weight must be in [0,1]");
  require(finite(o.macro_cell_density_weight) &&
              o.macro_cell_density_weight >= 0,
          "invalid macro cell density weight");
  require(finite(o.macro_cell_density_reach) && o.macro_cell_density_reach > 0,
          "invalid macro cell density reach");
  require(
      finite(o.macro_density_pull_weight) && o.macro_density_pull_weight >= 0 &&
          o.macro_density_pull_weight <= 1 &&
          finite(o.macro_density_pull_range) &&
          o.macro_density_pull_range > 0 && o.macro_density_pull_range <= 32,
      "invalid macro density pull controls");
  require(o.filler_max_count >= 0 && o.filler_max_count <= 100000000 &&
              finite(o.filler_macro_fraction) && o.filler_macro_fraction >= 0 &&
              o.filler_macro_fraction <= 1 &&
              !(o.filler_macro_fraction > 0 && o.filler_width > 0) &&
              finite(o.filler_cell_scale) && o.filler_cell_scale > 0 &&
              finite(o.filler_bin_fraction) && o.filler_bin_fraction > 0 &&
              o.filler_bin_fraction <= 1 && finite(o.filler_area_fraction) &&
              o.filler_area_fraction > 0 && o.filler_area_fraction <= 2,
          "invalid virtual filler budget");
  require(o.num_threads > 0 && o.num_threads <= 128,
          "thread count must be 1..128");
  require(finite(o.target_density) && o.target_density > 0 &&
              o.target_density <= 1,
          "density must be in (0,1]");
  require(o.macro_pull_settle_iterations >= 0,
          "negative macro pull settling budget");
  require(finite(o.stop_overflow) && o.stop_overflow >= 0 &&
              o.stop_overflow <= 1,
          "invalid overflow threshold");
  require(finite(o.initial_density_weight) && o.initial_density_weight > 0,
          "invalid density weight");
  require(finite(o.max_step) && o.max_step > 0 && o.max_step <= 1,
          "max_step must be in (0,1]");
  require(o.min_iterations >= 0, "negative minimum GP iterations");
  require(o.post_density_factor >= 0 && o.post_density_factor <= 100 &&
              finite(o.convergence_tolerance) && o.convergence_tolerance >= 0 &&
              o.convergence_tolerance <= 1,
          "invalid post-density convergence controls");
  require(finite(o.wirelength_gamma_end) && o.wirelength_gamma_end > 0 &&
              o.wirelength_gamma_end <= 1,
          "wirelength gamma ratio must be in (0,1]");
  require(finite(o.wirelength_gamma_start) && o.wirelength_gamma_start > 0 &&
              o.wirelength_gamma_start <= 1,
          "wirelength initial gamma ratio must be in (0,1]");
  require(finite(o.fixed_macro_overlap_weight) &&
              o.fixed_macro_overlap_weight >= 0,
          "invalid fixed-macro overlap weight");
  require(o.detailed_iterations >= 1 && o.detailed_iterations <= 100,
          "detailed placement iterations must be 1..100");
  require(!o.detailed_place || o.legalize,
          "detailed placement requires legalization");
  require(!o.macro_only_legalize || (o.legalize && !o.detailed_place),
          "macro-only legalization requires legalization without cell detail");
  require(finite(o.filler_width) && o.filler_width >= 0 &&
              finite(o.order_weight) && o.order_weight >= 0 &&
              finite(o.order_growth) && o.order_growth > 0 &&
              o.schedule_iterations >= 0,
          "invalid H2D GP hook options");
  validate_route_policy(o.route);
  require(!o.routability ||
              (finite(o.route_capacity_h) && o.route_capacity_h > 0 &&
               finite(o.route_capacity_v) && o.route_capacity_v > 0),
          "invalid routing capacity");
  require(finite(o.routability_macro_weight) &&
              o.routability_macro_weight >= 0 &&
              finite(o.routability_macro_range) &&
              o.routability_macro_range > 0 && o.routability_macro_range <= 16,
          "invalid macro routing pressure controls");
  require(finite(o.macro_overlap_weight) && o.macro_overlap_weight >= 0 &&
              finite(o.macro_jitter) && o.macro_jitter >= 0 &&
              o.macro_jitter <= .1,
          "invalid macro repulsion controls");
  require(
      finite(o.macro_align_weight) && o.macro_align_weight >= 0 &&
          finite(o.macro_align_sigma) && o.macro_align_sigma >= 0 &&
          finite(o.macro_align_overlap_boost) &&
          o.macro_align_overlap_boost >= 0 &&
          (!o.macro_align_weight || (o.macro_align_start >= 0 &&
                                     o.macro_align_start < o.macro_align_full)),
      "invalid macro alignment controls");
  require(finite(o.macro_snap_tolerance) && o.macro_snap_tolerance >= 0 &&
              o.macro_snap_weight > 0 && o.macro_snap_weight <= 1000000,
          "invalid macro snap controls");
  for (double value : {o.macro_snap_options.alignment_tolerance,
                       o.macro_snap_options.spacing_tolerance,
                       o.macro_snap_options.spacing_penalty,
                       o.macro_snap_options.symmetry_tolerance,
                       o.macro_snap_options.max_displacement}) {
    require(finite(value) && value >= 0, "invalid macro regularity control");
  }
  require(o.macro_pull_interval >= 0, "invalid virtual macro pull interval");
  require(finite(o.cohesion_weight) && o.cohesion_weight >= 0,
          "invalid cohesion weight");
  require(o.terminal_interval >= 0, "negative terminal update interval");
  require(o.movable_interval >= 0, "negative movable update interval");
  require(o.checkpoint_interval >= 0, "negative GP checkpoint interval");
  require(o.checkpoint_start >= 0, "negative GP checkpoint start");
  require(!o.observer || o.observer_interval > 0,
          "observer interval must be positive");
}
bool has_halo(const Design &d) {
  const auto positive = [](double value) { return value > 0; };
  return std::any_of(d.node_halo_x.begin(), d.node_halo_x.end(), positive) ||
         std::any_of(d.node_halo_y.begin(), d.node_halo_y.end(), positive);
}
double halo(const std::vector<double> &values, int i, double pitch) {
  return values.empty() ? 0 : std::ceil(values[i] / pitch - 1e-9) * pitch;
}
Design padded_design(const Design &source) {
  Design d = source;
  d.node_pad_x.assign(d.node_x.size(), 0);
  d.node_pad_y.assign(d.node_x.size(), 0);
  for (auto &group : d.cohesion_groups) {
    if (group.areas.empty()) {
      for (int i : group.nodes) {
        group.areas.push_back(source.node_size_x[i] * source.node_size_y[i]);
      }
    }
  }
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    const double x = halo(source.node_halo_x, i, d.layout.site_width);
    const double y = halo(source.node_halo_y, i, d.layout.row_height);
    if (!d.node_boxes.empty()) {
      d.node_boxes[i].xl -= x;
      d.node_boxes[i].xh -= x;
      d.node_boxes[i].yl -= y;
      d.node_boxes[i].yh -= y;
    }
    d.node_x[i] -= x;
    d.node_y[i] -= y;
    if (!d.node_legalize_anchor_x.empty()) {
      d.node_legalize_anchor_x[i] -= x;
      d.node_legalize_anchor_y[i] -= y;
    }
    d.node_pad_x[i] = x;
    d.node_pad_y[i] = y;
    d.node_size_x[i] += 2 * x;
    d.node_size_y[i] += 2 * y;
    for (int k = d.flat_node2pin_start[i]; k < d.flat_node2pin_start[i + 1];
         ++k) {
      const int p = d.flat_node2pin[k];
      d.pin_offset_x[p] += x;
      d.pin_offset_y[p] += y;
    }
  }
  for (auto &edge : d.order_edges) {
    const double a =
        edge.axis == Axis::x
            ? halo(source.node_halo_x, edge.from, d.layout.site_width)
            : halo(source.node_halo_y, edge.from, d.layout.row_height);
    const double b =
        edge.axis == Axis::x
            ? halo(source.node_halo_x, edge.to, d.layout.site_width)
            : halo(source.node_halo_y, edge.to, d.layout.row_height);
    // Preserve a relation on physical footprints under origin padding.
    edge.separation -= a + b;
  }
  d.node_halo_x.clear();
  d.node_halo_y.clear();
  return d;
}
} // namespace

namespace detail {
std::vector<MacroPullForce>
macro_pull_forces(const Design &d, const std::vector<MacroPull> &pulls,
                  const Gradient &wire, const Gradient &density,
                  double density_weight) {
  std::vector<MacroPullForce> result;
  if (pulls.empty()) {
    return result;
  }
  result.reserve(pulls.size());
  double reference = 0, minimum_reference = 0, length = 0;
  for (const auto &pull : pulls) {
    const int i = pull.node;
    const double pre = placement_preconditioner(d, i, density_weight);
    reference += (std::abs(wire.x[i] + density_weight * density.x[i]) +
                  std::abs(wire.y[i] + density_weight * density.y[i])) /
                 pre;
    minimum_reference += 1 / pre;
    length += std::max(d.node_size_x[i], d.node_size_y[i]);
  }
  reference = std::max(reference, minimum_reference) / pulls.size();
  length = std::max(1e-12, length / pulls.size());
  for (const auto &pull : pulls) {
    const int i = pull.node;
    const double dx = d.node_x[i] + d.node_size_x[i] / 2 - pull.x;
    const double dy = d.node_y[i] + d.node_size_y[i] / 2 - pull.y;
    // Shared length and identical far-target cap: no per-macro size in the
    // displacement factor. Multiplying by pre cancels the optimizer's division.
    const double scale = pull.weight * reference *
                         placement_preconditioner(d, i, density_weight) /
                         std::max(length, std::hypot(dx, dy));
    result.push_back({i, scale * dx, scale * dy});
  }
  return result;
}
Design unpad_snapshot(const Design &snapshot, const Design &reference) {
  if (snapshot.node_pad_x.empty()) {
    return snapshot;
  }
  require(snapshot.node_x.size() == reference.node_x.size() &&
              snapshot.node_pad_x.size() == snapshot.node_x.size() &&
              snapshot.node_pad_y.size() == snapshot.node_x.size(),
          "physical snapshot reference size mismatch");
  Design physical = snapshot;
  for (int i = 0; i < count(physical.node_x.size()); ++i) {
    const double hx = snapshot.node_pad_x[i], hy = snapshot.node_pad_y[i];
    physical.node_x[i] += hx;
    physical.node_y[i] += hy;
    if (!physical.node_legalize_anchor_x.empty()) {
      physical.node_legalize_anchor_x[i] += hx;
      physical.node_legalize_anchor_y[i] += hy;
    }
    physical.node_size_x[i] -= 2 * hx;
    physical.node_size_y[i] -= 2 * hy;
    for (int k = physical.flat_node2pin_start[i];
         k < physical.flat_node2pin_start[i + 1]; ++k) {
      const int pin = physical.flat_node2pin[k];
      physical.pin_offset_x[pin] -= hx;
      physical.pin_offset_y[pin] -= hy;
    }
  }
  physical.node_halo_x = reference.node_halo_x;
  physical.node_halo_y = reference.node_halo_y;
  physical.node_pad_x.clear();
  physical.node_pad_y.clear();
  physical.node_boxes = reference.node_boxes;
  physical.order_edges = reference.order_edges;
  physical.cohesion_groups = reference.cohesion_groups;
  return physical;
}
} // namespace detail

void Design::build_maps() {
  int nodes = count(node_x.size()), pins = count(pin2node.size());
  require(!flat_net2pin_start.empty(), "missing net CSR");
  int nets = count(flat_net2pin_start.size() - 1);
  // Build into temporary arrays so a malformed design does not leave partial
  // maps.
  std::vector<int> p2n(pins, -1);
  require(flat_net2pin.size() == static_cast<size_t>(pins) &&
              flat_net2pin_start.front() == 0 &&
              flat_net2pin_start.back() == pins,
          "invalid net CSR");
  for (int n = 0; n < nets; ++n) {
    int a = flat_net2pin_start[n], b = flat_net2pin_start[n + 1];
    require(a >= 0 && a <= b && b <= pins, "nonmonotonic net CSR");
    for (int j = a; j < b; ++j) {
      int p = flat_net2pin[j];
      require(p >= 0 && p < pins && p2n[p] == -1,
              "invalid or duplicate net pin");
      p2n[p] = n;
    }
  }
  std::vector<int> start(nodes + 1), flat(pins);
  for (int i : pin2node) {
    require(i >= 0 && i < nodes, "invalid pin2node");
    ++start[i + 1];
  }
  std::partial_sum(start.begin(), start.end(), start.begin());
  auto next = start;
  for (int p = 0; p < pins; ++p) {
    flat[next[pin2node[p]]++] = p;
  }
  if (net_weights.empty()) {
    net_weights.assign(nets, 1);
  }
  if (net_mask.empty()) {
    net_mask.resize(nets);
    for (int n = 0; n < nets; ++n) {
      net_mask[n] = flat_net2pin_start[n + 1] - flat_net2pin_start[n] > 1;
    }
  }
  pin2net = std::move(p2n);
  flat_node2pin_start = std::move(start);
  flat_node2pin = std::move(flat);
}

void Design::validate() const {
  int n = count(node_x.size()), p = count(pin2node.size()),
      nets = count(net_weights.size());
  require(num_movable >= 0 && num_movable <= n, "invalid movable count");
  const auto &l = layout;
  require(finite(l.xl) && finite(l.yl) && finite(l.xh) && finite(l.yh) &&
              l.xh > l.xl && l.yh > l.yl && finite(l.xh - l.xl) &&
              finite(l.yh - l.yl),
          "invalid core bounds");
  require(finite(l.site_width) && l.site_width > 0 && finite(l.row_height) &&
              l.row_height > 0,
          "invalid site/row pitch");
  require((l.xh - l.xl) / l.site_width < std::numeric_limits<int>::max() &&
              (l.yh - l.yl) / l.row_height < std::numeric_limits<int>::max(),
          "layout exceeds kernel site indices");
  require(aligned((l.xh - l.xl) / l.site_width) &&
              aligned((l.yh - l.yl) / l.row_height),
          "core must contain whole sites/rows");
  require(node_y.size() == node_x.size() &&
              node_size_x.size() == node_x.size() &&
              node_size_y.size() == node_x.size(),
          "node array size mismatch");
  require(node_legalize_weights.empty() ||
              node_legalize_weights.size() == node_x.size(),
          "legalization weight array size mismatch");
  require(node_legalize_last.empty() ||
              node_legalize_last.size() == node_x.size(),
          "legalization deferred array size mismatch");
  require((node_legalize_anchor_x.empty() && node_legalize_anchor_y.empty()) ||
              (node_legalize_anchor_x.size() == node_x.size() &&
               node_legalize_anchor_y.size() == node_x.size()),
          "legalization anchor array size mismatch");
  for (double weight : node_legalize_weights) {
    require(finite(weight) && weight >= 1,
            "invalid legalization displacement weight");
  }
  require((node_halo_x.empty() || node_halo_x.size() == node_x.size()) &&
              (node_halo_y.empty() || node_halo_y.size() == node_x.size()),
          "halo array size mismatch");
  require((node_pad_x.empty() || node_pad_x.size() == node_x.size()) &&
              (node_pad_y.empty() || node_pad_y.size() == node_x.size()),
          "pad array size mismatch");
  require(node_orientation.empty() || node_orientation.size() == node_x.size(),
          "orientation array size mismatch");
  require(node_orientation_fixed.empty() ||
              (node_orientation_fixed.size() == node_x.size() &&
               std::all_of(node_orientation_fixed.begin(),
                           node_orientation_fixed.end(),
                           [](auto v) { return v <= 1; })),
          "invalid frozen orientation mask");
  require(node_group_id.empty() || node_group_id.size() == node_x.size(),
          "group identity array size mismatch");
  require(node_boxes.empty() || node_boxes.size() == node_x.size(),
          "position box size mismatch");
  for (const auto *capacity :
       {&node_blocked_capacity_h, &node_blocked_capacity_v}) {
    require(capacity->empty() ||
                (capacity->size() == node_x.size() &&
                 std::all_of(capacity->begin(), capacity->end(),
                             [](double v) { return finite(v) && v >= 0; })),
            "invalid macro blocked capacity");
  }
  require(node_align_class.empty() ||
              (node_align_class.size() == node_x.size() &&
               std::all_of(node_align_class.begin(), node_align_class.end(),
                           [](int v) { return v >= -1; })),
          "invalid macro alignment classes");
  require(node_snap_class.empty() ||
              (node_snap_class.size() == node_x.size() &&
               std::all_of(node_snap_class.begin(), node_snap_class.end(),
                           [](int c) { return c >= -1; })),
          "invalid macro snap classes");
  require(node_macro.empty() || node_macro.size() == node_x.size(),
          "macro mask size mismatch");
  require(std::all_of(node_macro.begin(), node_macro.end(),
                      [](auto v) { return v <= 1; }),
          "invalid macro mask");
  require(node_placement_mask.empty() ||
              node_placement_mask.size() == node_x.size(),
          "placement mask size mismatch");
  for (int i = 0; i < n; ++i) {
    const auto mask = placement_mask(i);
    require(mask <= PlacementMask::soft, "invalid placement mask kind");
    if (mask != PlacementMask::none) {
      require(i >= num_movable && !node_macro.empty() && node_macro[i] &&
                  flat_node2pin_start.size() == node_x.size() + 1 &&
                  flat_node2pin_start[i] == flat_node2pin_start[i + 1],
              "placement masks must be fixed pinless macro geometry");
      require((node_halo_x.empty() || node_halo_x[i] == 0) &&
                  (node_halo_y.empty() || node_halo_y[i] == 0) &&
                  (node_orientation.empty() ||
                   node_orientation[i] == Orientation::N),
              "placement masks must have zero halo and N orientation");
    }
  }
  require(node_inflatable.empty() || node_inflatable.size() == node_x.size(),
          "inflation mask size mismatch");
  for (const auto &e : order_edges) {
    require(e.from >= 0 && e.from < n && e.to >= 0 && e.to < n &&
                e.from != e.to && (e.axis == Axis::x || e.axis == Axis::y) &&
                finite(e.separation) && finite(e.strength.weight) &&
                e.strength.weight >= 0,
            "invalid relative order");
  }
  for (int i = 0; i < n; ++i) {
    double w = node_size_x[i], h = node_size_y[i];
    if (!node_orientation.empty()) {
      require(valid_orientation(node_orientation[i]),
              "invalid node orientation");
      require((w > 0 && h > 0) || node_orientation[i] == Orientation::N,
              "terminal cannot be flipped");
    }
    const double hx = node_halo_x.empty() ? 0 : node_halo_x[i];
    const double hy = node_halo_y.empty() ? 0 : node_halo_y[i];
    require(finite(hx) && finite(hy) && hx >= 0 && hy >= 0,
            "invalid node halo");
    require((hx == 0 && hy == 0) || (w > 0 && h > 0),
            "halo requires a positive-area node");
    require(finite(node_x[i]) && finite(node_y[i]) && finite(w) && finite(h) &&
                w >= 0 && h >= 0 && finite(w * h),
            "nonfinite position or invalid node size");
    require(finite(node_x[i] + w) && finite(node_y[i] + h),
            "node bounds overflow");
    if (!node_boxes.empty()) {
      const auto &b = node_boxes[i];
      require(finite(b.xl) && finite(b.yl) && finite(b.xh) && finite(b.yh) &&
                  b.xl <= b.xh && b.yl <= b.yh,
              "invalid position box");
      if (i < num_movable) {
        require(std::max(l.xl, b.xl) <= std::min(l.xh - w, b.xh) &&
                    std::max(l.yl, b.yl) <= std::min(l.yh - h, b.yh),
                "position box outside core");
      }
    }
    if (!node_inflatable.empty()) {
      require(node_inflatable[i] <= 1, "invalid inflation mask");
    }
    if (i < num_movable) {
      require(w > 0 && h > 0 && w <= l.xh - l.xl && h <= l.yh - l.yl,
              "movable node does not fit core");
      // Macro legalization uses private site/row-rounded occupancy envelopes.
      // Keep the caller's physical dimensions and pin offsets unchanged.
    }
  }
  std::vector<bool> grouped(node_x.size(), false);
  require(centroid_targets.empty() ||
              centroid_targets.size() == centroid_groups.size(),
          "centroid target size mismatch");
  for (auto [x, y] : centroid_targets) {
    require(finite(x) && finite(y), "nonfinite frozen centroid target");
  }
  if (!centroid_groups.empty()) {
    require(centroid_node_area.size() == node_x.size() &&
                node_centroid_color.size() == node_x.size(),
            "centroid physical area/color metadata size mismatch");
    std::vector<bool> populated(centroid_color_names.size(), false);
    for (size_t n = 0; n < node_x.size(); ++n) {
      require(finite(centroid_node_area[n]) && centroid_node_area[n] >= 0,
              "invalid centroid physical area");
      const int c = node_centroid_color[n];
      require(c >= -1 && c < static_cast<int>(centroid_color_names.size()) &&
                  (c < 0 || centroid_node_area[n] > 0),
              "invalid centroid color identity");
      if (c >= 0) {
        populated[c] = true;
      }
    }
    for (size_t g = 0; g < centroid_groups.size(); ++g) {
      const auto &group = centroid_groups[g];
      require(group.parent >= -1 && group.parent < static_cast<int>(g) &&
                  group.family >= 0 && group.repulsion_scope >= -1,
              "invalid centroid family/parent/repulsion scope");
      require(finite(group.repulsion_area) &&
                  (group.repulsion_area == -1 || group.repulsion_area >= 0),
              "invalid centroid repulsion area");
      require(!group.nodes.empty() || !group.colors.empty(),
              "empty centroid group");
      require(
          (group.areas.empty() || group.areas.size() == group.nodes.size()) &&
              (group.offsets.empty() ||
               group.offsets.size() == group.nodes.size()),
          "centroid logical member metadata size mismatch");
      std::set<int> nodes, colors;
      for (size_t k = 0; k < group.nodes.size(); ++k) {
        const int n = group.nodes[k];
        const bool unique = nodes.insert(n).second;
        require(
            n >= 0 && n < static_cast<int>(node_x.size()) &&
                (unique || (!group.areas.empty() && !group.offsets.empty())) &&
                std::find(group.colors.begin(), group.colors.end(),
                          node_centroid_color[n]) == group.colors.end() &&
                centroid_node_area[n] > 0,
            "invalid or duplicate explicit centroid member");
        if (!group.areas.empty()) {
          require(finite(group.areas[k]) && group.areas[k] > 0,
                  "invalid centroid logical member area");
        }
        if (!group.offsets.empty()) {
          require(finite(group.offsets[k].first) &&
                      finite(group.offsets[k].second),
                  "nonfinite centroid logical member offset");
        }
      }
      for (int c : group.colors) {
        require(c >= 0 && c < static_cast<int>(centroid_color_names.size()) &&
                    colors.insert(c).second && populated[c],
                "invalid or duplicate centroid anchor color");
      }
    }
  }
  for (const auto &group : cohesion_groups) {
    require(finite(group.weight) && group.weight >= 0,
            "invalid cohesion group weight");
    require(group.areas.empty() || group.areas.size() == group.nodes.size(),
            "cohesion area size mismatch");
    require(group.offsets.empty() || group.offsets.size() == group.nodes.size(),
            "cohesion offset size mismatch");
    for (size_t k = 0; k < group.nodes.size(); ++k) {
      int i = group.nodes[k];
      require(i >= 0 && i < n && (!grouped[i] || !group.offsets.empty()),
              "invalid or overlapping cohesion membership");
      if (!group.offsets.empty()) {
        require(finite(group.offsets[k].first) &&
                    finite(group.offsets[k].second),
                "invalid cohesion member offset");
      }
      grouped[i] = true;
      double a = group.areas.empty() ? node_size_x[i] * node_size_y[i]
                                     : group.areas[k];
      require(finite(a) && a > 0 && node_size_x[i] > 0 && node_size_y[i] > 0,
              "cohesion requires physical positive-area nodes");
    }
  }
  require(pin2net.size() == pin2node.size() &&
              pin_offset_x.size() == pin2node.size() &&
              pin_offset_y.size() == pin2node.size() &&
              net_mask.size() == net_weights.size(),
          "pin/net array size mismatch");
  for (int i = 0; i < p; ++i) {
    require(pin2node[i] >= 0 && pin2node[i] < n && pin2net[i] >= 0 &&
                pin2net[i] < nets,
            "invalid pin map index");
    require(finite(pin_offset_x[i]) && finite(pin_offset_y[i]),
            "nonfinite pin offset");
  }
  for (int i = 0; i < nets; ++i) {
    require(finite(net_weights[i]) && net_weights[i] >= 0 && net_mask[i] <= 1,
            "invalid net weight/mask");
  }
  csr_check(flat_net2pin, flat_net2pin_start, nets, p, pin2net);
  csr_check(flat_node2pin, flat_node2pin_start, n, p, pin2node);
  for (int i = 0; i < nets; ++i) {
    require(!net_mask[i] ||
                flat_net2pin_start[i + 1] - flat_net2pin_start[i] >= 2,
            "active net needs at least two pins");
  }
  detail::validate_pin_timing(*this);
  reservation_charge(*this);
}

double density_overflow(const Design &d, int bx, int by, double density,
                        bool scale_macro) {
  d.validate();
  require(bx > 0 && by > 0 && bx <= 4096 && by <= 4096,
          "invalid density metric bins");
  require(finite(density) && density > 0 && density <= 1,
          "invalid density metric target");
  return overflow(d, bx, by, density, scale_macro, false, 1, false);
}

double reservation_density_overflow(const Design &d, int bx, int by,
                                    double density, bool scale_macro) {
  d.validate();
  require(bx > 0 && by > 0 && bx <= 4096 && by <= 4096,
          "invalid density metric bins");
  require(finite(density) && density > 0 && density <= 1,
          "invalid density metric target");
  return overflow(d, bx, by, density, scale_macro, false, 1, true);
}
double hpwl(const Design &d) {
  d.validate();
  double wl = 0;
  for (int n = 0; n < count(d.net_weights.size()); ++n) {
    if (!d.net_mask[n]) {
      continue;
    }
    double xmin = std::numeric_limits<double>::infinity(), ymin = xmin,
           xmax = -xmin, ymax = -xmin;
    for (int j = d.flat_net2pin_start[n]; j < d.flat_net2pin_start[n + 1];
         ++j) {
      int p = d.flat_net2pin[j], i = d.pin2node[p];
      double x = d.node_x[i] + d.pin_offset_x[p],
             y = d.node_y[i] + d.pin_offset_y[p];
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    wl += d.net_weights[n] * (xmax - xmin + ymax - ymin);
  }
  return wl;
}

bool macros_are_legal(const Design &d, std::string *reason) {
  d.validate();
  Design macros;
  macros.layout = d.layout;
  // Retain all fixed obstacles. Movable cells are provisional and must neither
  // obstruct the macros nor require row/site legality during macro planning.
  const auto append = [&](int i) {
    macros.node_x.push_back(d.node_x[i]);
    macros.node_y.push_back(d.node_y[i]);
    macros.node_size_x.push_back(d.node_size_x[i]);
    macros.node_size_y.push_back(d.node_size_y[i]);
    macros.node_macro.push_back(d.node_macro.empty() ? 0 : d.node_macro[i]);
    macros.node_placement_mask.push_back(d.placement_mask(i));
    if (!d.node_halo_x.empty()) {
      macros.node_halo_x.push_back(d.node_halo_x[i]);
    }
    if (!d.node_halo_y.empty()) {
      macros.node_halo_y.push_back(d.node_halo_y[i]);
    }
    if (!d.node_orientation.empty()) {
      macros.node_orientation.push_back(d.node_orientation[i]);
    }
  };
  for (int i = 0; i < d.num_movable; ++i) {
    if (!d.node_macro.empty() && d.node_macro[i]) {
      append(i);
    }
  }
  macros.num_movable = macros.node_x.size();
  for (int i = d.num_movable; i < count(d.node_x.size()); ++i) {
    append(i);
  }
  macros.build_maps();
  return is_legal(macros, reason);
}

bool is_legal(const Design &d, std::string *reason) {
  if (reason) {
    reason->clear();
  }
  d.validate();
  if (has_halo(d)) {
    return is_legal(padded_design(d), reason);
  }
  const auto &l = d.layout;
  constexpr double eps = 1e-6;
  struct Interval {
    double x, end;
    int id;
  };
  std::vector<std::vector<Interval>> rows(
      static_cast<size_t>(std::round((l.yh - l.yl) / l.row_height)));
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
           h = d.node_size_y[i];
    if (!w || !h || d.placement_mask(i) == PlacementMask::soft) {
      continue;
    }
    if (i < d.num_movable && (x < l.xl - eps || y < l.yl - eps ||
                              x + w > l.xh + eps || y + h > l.yh + eps ||
                              ((d.node_macro.empty() || !d.node_macro[i]) &&
                               (!aligned((x - l.xl) / l.site_width) ||
                                !aligned((y - l.yl) / l.row_height))))) {
      if (reason) {
        *reason = "bounds or grid: node " + std::to_string(i) +
                  " x=" + std::to_string(x) + " y=" + std::to_string(y) +
                  " site=" + std::to_string(l.site_width) +
                  " row=" + std::to_string(l.row_height);
      }
      return false;
    }
    // Fixed nodes may lie anywhere (validate() only needs finite coordinates).
    // Clamp in floating point: converting an out-of-range double to int is
    // undefined and indexed rows[INT_MIN] for a node far above the core.
    const double row_count = static_cast<double>(rows.size());
    const int first = static_cast<int>(
        std::clamp(std::floor((y - l.yl + eps) / l.row_height), 0., row_count));
    const int last = static_cast<int>(std::clamp(
        std::ceil((y + h - l.yl - eps) / l.row_height), 0., row_count));
    for (int r = first; r < last; ++r) {
      rows[r].push_back({x, x + w, i});
    }
  }
  const bool movable_macros =
      !d.node_macro.empty() &&
      std::any_of(d.node_macro.begin(), d.node_macro.begin() + d.num_movable,
                  [](auto macro) { return macro != 0; });
  for (auto &row : rows) {
    std::sort(row.begin(), row.end(), [](auto a, auto b) { return a.x < b.x; });
    if (movable_macros) {
      // Off-row macro rectangles can share a bucket without sharing physical Y.
      // Keep the usual fixed/fixed exemption, but check exact Y intersections.
      std::multimap<double, int> active, movable;
      for (const auto v : row) {
        while (!active.empty() && active.begin()->first <= v.x + eps) {
          active.erase(active.begin());
        }
        while (!movable.empty() && movable.begin()->first <= v.x + eps) {
          movable.erase(movable.begin());
        }
        const auto &candidates = v.id < d.num_movable ? active : movable;
        for (const auto &[end, other] : candidates) {
          (void)end;
          if (d.node_y[other] < d.node_y[v.id] + d.node_size_y[v.id] - eps &&
              d.node_y[v.id] < d.node_y[other] + d.node_size_y[other] - eps) {
            if (reason) {
              *reason = "overlap: node " + std::to_string(other) + "; node " +
                        std::to_string(v.id);
            }
            return false;
          }
        }
        active.emplace(v.end, v.id);
        if (v.id < d.num_movable) {
          movable.emplace(v.end, v.id);
        }
      }
      continue;
    }
    double end_all = -std::numeric_limits<double>::infinity(),
           end_movable = end_all;
    int end_all_node = -1, end_movable_node = -1;
    for (auto v : row) {
      if (v.id < d.num_movable ? v.x < end_all - eps
                               : v.x < end_movable - eps) {
        if (reason) {
          const int other =
              v.id < d.num_movable ? end_all_node : end_movable_node;
          *reason = "overlap: node " + std::to_string(other) +
                    " x=" + std::to_string(d.node_x[other]) +
                    " y=" + std::to_string(d.node_y[other]) +
                    " w=" + std::to_string(d.node_size_x[other]) +
                    " h=" + std::to_string(d.node_size_y[other]) + "; node " +
                    std::to_string(v.id) +
                    " x=" + std::to_string(d.node_x[v.id]) +
                    " y=" + std::to_string(d.node_y[v.id]) +
                    " w=" + std::to_string(d.node_size_x[v.id]) +
                    " h=" + std::to_string(d.node_size_y[v.id]);
        }
        return false;
      }
      if (v.end > end_all) {
        end_all = v.end;
        end_all_node = v.id;
      }
      if (v.id < d.num_movable && v.end > end_movable) {
        end_movable = v.end;
        end_movable_node = v.id;
      }
    }
  }
  return true;
}

const char *stage_name(PlacementStage stage) {
  switch (stage) {
  case PlacementStage::centroid_start:
    return "centroid attraction start";
  case PlacementStage::initial:
    return "initial";
  case PlacementStage::global_placement:
    return "global placement";
  case PlacementStage::legalization:
    return "legalization";
  case PlacementStage::detailed_placement:
    return "detailed placement";
  case PlacementStage::finished:
    return "finished";
  }
  return "unknown";
}

const char *route_model_name(RouteModel model) {
  switch (model) {
  case RouteModel::off:
    return "off";
  case RouteModel::rudy:
    return "rudy";
  case RouteModel::rsmt:
    return "rsmt";
  case RouteModel::automatic:
    return "auto";
  case RouteModel::uroute:
    return "uroute";
  }
  throw std::invalid_argument("invalid route model");
}
RouteModel route_model_from_name(std::string_view name) {
  for (auto model : {RouteModel::off, RouteModel::rudy, RouteModel::rsmt,
                     RouteModel::automatic, RouteModel::uroute})
    if (name == route_model_name(model))
      return model;
  throw std::invalid_argument("invalid route model: " + std::string(name));
}
void validate_route_policy(const RoutePolicy &p) {
  require(finite(p.switch_overflow) && finite(p.target_overflow) &&
              p.target_overflow >= 0 &&
              p.switch_overflow >= p.target_overflow && p.hold >= 1 &&
              p.hold <= 1000 && p.degree >= 2 && p.degree <= 1000 &&
              p.accuracy >= 1 && p.accuracy <= 10 && p.probe_iterations >= 1 &&
              p.probe_iterations <= 10000,
          "invalid route switch, target, hold, degree, accuracy or probe "
          "budget");
  route_model_name(p.model);
  require(p.probe_model == RouteModel::rsmt ||
              p.probe_model == RouteModel::uroute,
          "route probe model must be rsmt or uroute");
  require(finite(p.uroute_switch_overflow) && p.uroute_switch_overflow >= 0 &&
              p.uroute_hold >= 1 && p.uroute_hold <= 1000 &&
              (!(p.model == RouteModel::automatic && p.uroute_ladder) ||
               p.uroute_switch_overflow <= p.switch_overflow),
          "invalid URoute switch or hold");
}
bool uroute_requested(const RoutePolicy &p) {
  return p.model == RouteModel::uroute ||
         (p.model == RouteModel::automatic && p.uroute_ladder) ||
         (p.probe && p.probe_model == RouteModel::uroute);
}
bool route_latches(const RoutePolicy &p) {
  return p.model == RouteModel::automatic ||
         (p.probe &&
          (p.model == RouteModel::rudy || p.probe_model == RouteModel::uroute));
}
Result place(Design &design, const Options &options, GpState *state) {
  design.validate();
  Options o = options;
  if (o.route.model == RouteModel::off) {
    o.routability = false;
  }
  if (o.automatic_bins) {
    const auto bins = detail::automatic_bins(design);
    o.bins_x = bins.first;
    o.bins_y = bins.second;
    o.automatic_bins = false;
  }
  if (o.dynamic_bins && design.node_pad_x.empty()) {
    require(o.bin_max >= 32 && o.bin_max <= 1024 &&
                std::has_single_bit(unsigned(o.bin_max)),
            "invalid dynamic bin cap");
    const auto bins = detail::dynamic_bins(design, o.target_density, o.bin_max);
    o.bins_x = bins.first;
    o.bins_y = bins.second;
  }
  if (o.gp_search.enabled && !o.schedule_iterations) {
    o.schedule_iterations = o.max_iterations;
  }
  options_check(o);
  o.magnet.validate();
  if (o.magnet.enabled) {
    require(!state && !o.magnet_runtime,
            "macro magnet requires one solve without external continuation");
    require(o.legalize, "macro magnet requires placement with legalization");
    return detail::magnet_place(design, o);
  }
  if (o.gp_search.enabled) {
    require(!state && !o.trial_tuning && o.trial_timing && o.add_fillers &&
                o.routability && o.improved_routability &&
                o.route.model == RouteModel::rsmt &&
                !o.macro_transparency.enabled() && !o.terminal_update &&
                !o.movable_update && !o.macro_pull_update &&
                o.centroid_weight == 0 && o.macro_align_weight == 0 &&
                o.macro_cell_density_weight == 0 &&
                o.macro_density_pull_weight == 0 &&
                o.local_density_jitter == 0 && o.local_density_weight > 0 &&
                o.local_density_weight <= 1 &&
                std::int64_t(o.gp_search.start) + o.gp_search.steps + 20 <=
                    o.max_iterations,
            "GP A/B search requires final fixed-macro RSMT placement, fillers, "
            "timing and settling budget");
    require(design.node_macro.empty() ||
                std::none_of(design.node_macro.begin(),
                             design.node_macro.begin() + design.num_movable,
                             [](auto v) { return v != 0; }),
            "GP A/B search cannot move macros");
  }
  if (o.checkpoint_rollback) {
    require(o.checkpoint && o.checkpoint_interval > 0 && !state &&
                !o.gp_search.enabled && !o.trial_tuning && !o.magnet_runtime &&
                !o.macro_transparency.enabled() && !o.vanished_memory &&
                !o.terminal_update && !o.movable_update &&
                !o.macro_pull_update && !o.routability &&
                !o.post_density_factor && o.local_density_weight == 0 &&
                o.centroid_weight == 0 && o.centroid_parent_weight == 0 &&
                o.centroid_align_weight == 0 && o.macro_align_weight == 0 &&
                o.macro_cell_density_weight == 0 &&
                o.macro_density_pull_weight == 0 && o.order_weight == 0 &&
                o.cohesion_weight == 0 && !o.macro_motion_observer &&
                !o.macro_step_observer,
            "GP rollback requires a fixed-footprint checkpoint solve without "
            "mutable auxiliary policies");
    require(design.node_macro.empty() ||
                std::none_of(design.node_macro.begin(),
                             design.node_macro.begin() + design.num_movable,
                             [](auto v) { return v != 0; }),
            "GP rollback cannot move macros");
  }
  if (o.trial_tuning) {
    require(state && state->signature &&
                state->iteration >=
                    std::max(220., o.vanished_memory
                                       ? memory_settled_iteration(o)
                                       : 220.),
            "GP trial requires an opaque full optimizer checkpoint");
  }
  const int aux_x = o.auxiliary_bins ? o.auxiliary_bins : o.bins_x;
  const int aux_y = o.auxiliary_bins ? o.auxiliary_bins : o.bins_y;
  require(
      finite(o.centroid_weight) && o.centroid_weight >= 0 &&
          finite(o.centroid_attraction_outside) &&
          o.centroid_attraction_outside >= 0 &&
          finite(o.centroid_repulsion_weight) &&
          o.centroid_repulsion_weight >= 0 && o.centroid_ramp > 0 &&
          (!o.centroid_repulsion_weight || o.centroid_weight > 0) &&
          finite(o.centroid_stop_overlap) && o.centroid_stop_overlap >= 0 &&
          o.centroid_stop_overlap < 1 &&
          (!o.centroid_adaptive_repulsion || o.centroid_repulsion_weight > 0),
      "invalid centroid attraction/repulsion controls");
  const bool centroids =
      o.centroid_weight > 0 && !design.centroid_groups.empty();
  const bool separate_circles = centroids && o.centroid_adaptive_repulsion;
  require(!centroids || !o.centroid_anchor_macros_only ||
              !design.centroid_targets.empty(),
          "macro-only centroid anchors require frozen targets");
  require(finite(o.centroid_align_weight) && o.centroid_align_weight >= 0 &&
              finite(o.centroid_align_sigma) && o.centroid_align_sigma >= 0 &&
              (!o.centroid_align_weight ||
               (o.centroid_weight > 0 && design.centroid_targets.empty() &&
                o.centroid_align_start >= 0 &&
                o.centroid_align_full > o.centroid_align_start)),
          "invalid centroid alignment controls");
  const bool align_centroids = centroids && o.centroid_align_weight > 0;
  require(finite(o.centroid_parent_weight) && o.centroid_parent_weight >= 0 &&
              finite(o.centroid_parent_outside) &&
              o.centroid_parent_outside >= 0 &&
              (!o.centroid_parent_weight ||
               (o.centroid_weight > 0 && design.centroid_targets.empty())),
          "invalid parent-centroid attraction controls");
  const bool parent_centroids = centroids && o.centroid_parent_weight > 0;
  require(
      !centroids || !state,
      "centroid placement requires one GP solve without external continuation");
  require(!o.checkpoint || o.checkpoint_interval == 0 || !state,
          "GP checkpoints require no external continuation");
  require(!o.terminal_update || o.terminal_interval == 0 || !state,
          "terminal updates require no external continuation");
  require(!o.movable_update || o.movable_interval == 0 || !state,
          "movable updates require no external continuation");
  require(
      !(o.macro_transparency.enabled() || o.vanished_memory) || !state ||
          (o.trial_tuning && !o.macro_transparency.enabled()),
      "macro transparency requires one GP solve without external continuation");
  if (o.macro_transparency.enabled() || o.vanished_memory) {
    require(design.node_macro.empty() ||
                std::none_of(design.node_macro.begin(),
                             design.node_macro.begin() + design.num_movable,
                             [](auto value) { return value != 0; }),
            "macro transparency is final-placement only: macros must be fixed");
  }
  if (has_halo(design)) {
    Design padded = padded_design(design);
    Options control = o;
    if (o.observer) {
      control.observer = [&](const Design &snapshot, PlacementStage stage,
                             int iteration) {
        Design physical = design;
        for (int i = 0; i < physical.num_movable; ++i) {
          physical.node_x[i] =
              snapshot.node_x[i] +
              halo(design.node_halo_x, i, design.layout.site_width);
          physical.node_y[i] =
              snapshot.node_y[i] +
              halo(design.node_halo_y, i, design.layout.row_height);
        }
        if (o.checkpoint && o.checkpoint_interval > 0) {
          physical.net_weights = snapshot.net_weights;
          physical.pin_timing = snapshot.pin_timing;
          physical.pin_offset_x = snapshot.pin_offset_x;
          physical.pin_offset_y = snapshot.pin_offset_y;
          for (int i = 0; i < count(physical.node_x.size()); ++i) {
            const double hx =
                halo(design.node_halo_x, i, design.layout.site_width);
            const double hy =
                halo(design.node_halo_y, i, design.layout.row_height);
            physical.node_size_x[i] = snapshot.node_size_x[i] - 2 * hx;
            physical.node_size_y[i] = snapshot.node_size_y[i] - 2 * hy;
            for (int k = physical.flat_node2pin_start[i];
                 k < physical.flat_node2pin_start[i + 1]; ++k) {
              const int pin = physical.flat_node2pin[k];
              physical.pin_offset_x[pin] -= hx;
              physical.pin_offset_y[pin] -= hy;
            }
          }
        }
        o.observer(physical, stage, iteration);
      };
    }
    Design halo_model = design;
    const auto unpad = [&](const Design &snapshot) {
      return detail::unpad_snapshot(snapshot, halo_model);
    };
    if (o.trial_timing) {
      control.trial_timing = [&](const Design &snapshot) {
        return o.trial_timing(unpad(snapshot));
      };
    }
    if (o.macro_pull_update && o.macro_pull_interval > 0) {
      control.macro_pull_update = [&](const Design &snapshot, int iteration) {
        return o.macro_pull_update(unpad(snapshot), iteration);
      };
    }
    if (o.checkpoint && o.checkpoint_interval > 0) {
      if (o.observer) {
        control.observer = [&](const Design &snapshot, PlacementStage stage,
                               int iteration) {
          o.observer(unpad(snapshot), stage, iteration);
        };
      }
      control.checkpoint = [&](const Design &snapshot, int iteration) {
        auto physical = unpad(snapshot);
        auto update = o.checkpoint(physical, iteration);
        if (update.replacement) {
          // Fail on a malformed host netlist before padded_design indexes it.
          Design candidate = std::move(*update.replacement);
          candidate.validate();
          halo_model = std::move(candidate);
          update.replacement = padded_design(halo_model);
        } else {
          for (int i = 0; i < count(physical.node_x.size()); ++i) {
            const double hx =
                halo(halo_model.node_halo_x, i, design.layout.site_width);
            const double hy =
                halo(halo_model.node_halo_y, i, design.layout.row_height);
            if (!update.node_size_x.empty()) {
              update.node_size_x.at(i) += 2 * hx;
            }
            if (!update.node_size_y.empty()) {
              update.node_size_y.at(i) += 2 * hy;
            }
            for (int k = physical.flat_node2pin_start[i];
                 k < physical.flat_node2pin_start[i + 1]; ++k) {
              const int pin = physical.flat_node2pin[k];
              if (!update.pin_offset_x.empty()) {
                update.pin_offset_x.at(pin) += hx;
              }
              if (!update.pin_offset_y.empty()) {
                update.pin_offset_y.at(pin) += hy;
              }
            }
          }
        }
        return update;
      };
    }
    if (o.before_legalize) {
      control.before_legalize = [&](Design &snapshot) {
        auto physical = unpad(snapshot);
        o.before_legalize(physical);
        physical.validate();
        halo_model = std::move(physical);
        snapshot = padded_design(halo_model);
      };
    }
    if (o.after_macro_legalize) {
      control.after_macro_legalize = [&](Design &snapshot) {
        auto physical = unpad(snapshot);
        o.after_macro_legalize(physical);
        physical.validate();
        halo_model = std::move(physical);
        snapshot = padded_design(halo_model);
      };
    }
    if (o.terminal_update && o.terminal_interval > 0) {
      control.terminal_update = [&](const Design &snapshot, int iteration) {
        return o.terminal_update(unpad(snapshot), iteration);
      };
      if (o.observer) {
        control.observer = [&](const Design &snapshot, PlacementStage stage,
                               int iteration) {
          o.observer(unpad(snapshot), stage, iteration);
        };
      }
    }
    if (o.movable_update && o.movable_interval > 0) {
      control.movable_update = [&](const Design &snapshot, int iteration) {
        auto updates = o.movable_update(unpad(snapshot), iteration);
        for (auto &position : updates) {
          require(position.node >= 0 && position.node < snapshot.num_movable,
                  "invalid movable update index");
          position.x -= halo(halo_model.node_halo_x, position.node,
                             snapshot.layout.site_width);
          position.y -= halo(halo_model.node_halo_y, position.node,
                             snapshot.layout.row_height);
        }
        return updates;
      };
    }
    const auto result = place(padded, control, state);
    if ((o.checkpoint && o.checkpoint_interval > 0) ||
        (o.terminal_update && o.terminal_interval > 0) ||
        (o.movable_update && o.movable_interval > 0) || o.before_legalize ||
        o.after_macro_legalize) {
      design = unpad(padded);
      return result;
    }
    for (int i = 0; i < design.num_movable; ++i) {
      design.node_x[i] = padded.node_x[i] +
                         halo(design.node_halo_x, i, design.layout.site_width);
      design.node_y[i] = padded.node_y[i] +
                         halo(design.node_halo_y, i, design.layout.row_height);
    }
    return result;
  }
  GpState next;
  const bool warm = state && state->signature;
  // FNV over geometry/connectivity. Coordinates, pin offsets and weights may
  // evolve between H2D checkpoints; changing the solver topology cannot reuse
  // it.
  std::uint64_t signature = 14695981039346656037ULL;
  auto salt = [&](std::uint64_t v) {
    signature = (signature ^ v) * 1099511628211ULL;
  };
  salt(design.num_movable);
  salt(o.bins_x);
  salt(o.bins_y);
  salt(o.add_fillers);
  if (o.add_fillers) {
    if (o.filler_macro_fraction > 0) {
      salt(std::bit_cast<std::uint64_t>(o.filler_macro_fraction));
    }
    salt(std::bit_cast<std::uint64_t>(o.filler_cell_scale));
    salt(std::bit_cast<std::uint64_t>(o.filler_bin_fraction));
    salt(o.filler_max_count);
    salt(std::bit_cast<std::uint64_t>(o.filler_area_fraction));
    if (o.filler_cell_space) {
      salt(0x63656c6c73706163ULL);
    }
  }
  if (o.dynamic_bins || o.local_density_weight || o.macro_cell_density_weight ||
      o.macro_density_pull_weight || o.auxiliary_bins) {
    salt(o.dynamic_bins);
    salt(o.bin_max);
    salt(o.auxiliary_bins);
    salt(o.local_density_bins);
    for (double v : {o.local_density_weight, o.local_density_range,
                     o.local_density_target, o.local_density_jitter}) {
      salt(std::bit_cast<std::uint64_t>(v));
    }
  }
  if (o.macro_cell_density_weight) {
    salt(0x6d6163726f726563ULL); // long-range, stage-calibrated cell pressure
    salt(std::bit_cast<std::uint64_t>(o.macro_cell_density_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_cell_density_reach));
  }
  if (o.timing_spread_weight) {
    salt(0x74696d6573707264ULL);
    salt(std::bit_cast<std::uint64_t>(o.timing_spread_weight));
  }
  if (o.macro_density_pull_weight) {
    salt(std::bit_cast<std::uint64_t>(o.macro_density_pull_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_density_pull_range));
  }
  if (o.density_upper_bound) {
    salt(o.density_upper_bound);
  }
  if (o.min_iterations != 40) {
    salt(o.min_iterations);
  }
  if (o.post_density_factor) {
    salt(0x706f737464656e73ULL);
    salt(o.post_density_factor);
    salt(std::bit_cast<std::uint64_t>(o.convergence_tolerance));
  }
  if (o.wirelength_gamma_end != 0.05) {
    salt(std::bit_cast<std::uint64_t>(o.wirelength_gamma_end));
  }
  if (o.wirelength_gamma_start != 1) {
    salt(0x67616d6d61737461ULL);
    salt(std::bit_cast<std::uint64_t>(o.wirelength_gamma_start));
  }
  if (o.vanished_memory) {
    salt(std::bit_cast<std::uint64_t>(o.vanished_memory_initial_scale));
    salt(o.vanished_memory_growth_iterations);
    salt(o.vanished_memory_start);
  }
  salt(o.seed);
  for (double v :
       {design.layout.xl, design.layout.yl, design.layout.xh, design.layout.yh,
        design.layout.site_width, design.layout.row_height, o.target_density,
        o.filler_width}) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  for (double v : design.node_size_x) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  for (double v : design.node_size_y) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  if (!design.node_orientation_fixed.empty()) {
    salt(0x7269676964666978ULL);
    for (auto fixed : design.node_orientation_fixed) {
      salt(fixed);
    }
  }
  for (int v : design.pin2node) {
    salt(v);
  }
  for (int v : design.flat_net2pin) {
    salt(v);
  }
  for (int v : design.flat_net2pin_start) {
    salt(v);
  }
  if (o.cohesion_weight > 0) {
    salt(std::bit_cast<std::uint64_t>(o.cohesion_weight));
    for (const auto &group : design.cohesion_groups) {
      salt(group.nodes.size());
      salt(std::bit_cast<std::uint64_t>(group.weight));
      for (int i : group.nodes) {
        salt(i);
      }
      for (double a : group.areas) {
        salt(std::bit_cast<std::uint64_t>(a));
      }
      if (!group.offsets.empty()) {
        salt(0x72696769646f6666ULL);
        for (auto [x, y] : group.offsets) {
          salt(std::bit_cast<std::uint64_t>(x));
          salt(std::bit_cast<std::uint64_t>(y));
        }
      }
    }
  }
  for (auto mask : design.node_placement_mask) {
    salt(static_cast<unsigned>(mask));
  }
  if (o.fixed_macro_overlap_weight > 0) {
    salt(std::uint64_t{0x66697865646d6163});
    salt(std::bit_cast<std::uint64_t>(o.fixed_macro_overlap_weight));
    for (auto v : design.node_macro) {
      salt(v);
    }
  }
  if (o.exact_cell_density) {
    salt(std::uint64_t{0x657861637463656c});
  }
  if (o.scale_macro_density || o.exact_macro_density || o.adaptive_density ||
      o.macro_overlap_weight || o.macro_jitter || o.improved_routability) {
    salt(o.scale_macro_density);
    salt(o.exact_macro_density);
    salt(o.adaptive_density);
    salt(o.improved_routability);
    if (o.improved_routability) {
      salt(std::bit_cast<std::uint64_t>(o.route_capacity_h));
      salt(std::bit_cast<std::uint64_t>(o.route_capacity_v));
      for (const auto *values :
           {&design.node_blocked_capacity_h, &design.node_blocked_capacity_v}) {
        for (double v : *values) {
          salt(std::bit_cast<std::uint64_t>(v));
        }
      }
      for (auto v : design.node_inflatable) {
        salt(v);
      }
    }
    salt(std::bit_cast<std::uint64_t>(o.macro_overlap_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_jitter));
    for (auto m : design.node_macro) {
      salt(m);
    }
  }
  if (o.improved_routability && !o.routability_charge_compensation) {
    salt(0x524f555445434847ULL);
  }
  if (o.routability && o.routability_macro_weight > 0) {
    salt(0x525544594d414352ULL);
    salt(std::bit_cast<std::uint64_t>(o.routability_macro_weight));
    salt(std::bit_cast<std::uint64_t>(o.routability_macro_range));
  }
  if (o.macro_align_weight > 0) {
    if (o.macro_align_global) {
      salt(0x5354524950414c49ULL);
    }
    salt(std::bit_cast<std::uint64_t>(o.macro_align_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_align_sigma));
    salt(std::bit_cast<std::uint64_t>(o.macro_align_overlap_boost));
    salt(o.macro_align_start);
    salt(o.macro_align_full);
    salt(o.macro_align_nearest);
    for (auto c : design.node_align_class) {
      salt(c);
    }
    for (auto m : design.node_macro) {
      salt(m);
    }
  }
  salt(static_cast<std::uint64_t>(o.route.model));
  salt(std::bit_cast<std::uint64_t>(o.route.switch_overflow));
  salt(std::bit_cast<std::uint64_t>(o.route.target_overflow));
  salt(o.route.hold);
  salt(o.route.degree);
  salt(o.route.accuracy);
  salt(o.route.probe);
  salt(o.route.probe_iterations);
  salt(static_cast<std::uint64_t>(o.route.probe_model));
  salt(std::bit_cast<std::uint64_t>(o.route.uroute_switch_overflow));
  salt(o.route.uroute_hold);
  salt(o.route.uroute_ladder);
  salt(bool(o.route_physical_snapshot));
  // A continuation must keep its immutable stage context, including support.
  salt(reinterpret_cast<std::uintptr_t>(o.uroute.get()));
  for (unsigned char c : o.uroute_identity)
    salt(c);
  if (warm && state->signature != signature) {
    throw std::invalid_argument("GP continuation options changed");
  }
  require(!warm || state->pin_timing == design.pin_timing,
          "GP continuation pin timing changed outside a checkpoint");
  require(!warm || state->density_reservations == design.density_reservations,
          "GP continuation density reservations changed outside a checkpoint");
  Result r;
  r.uroute_skipped = o.uroute_skipped;
  if (o.uroute)
    r.uroute_sample_ids = o.uroute->sample_ids();
  else if (r.uroute_skipped.empty() && uroute_requested(o.route))
    r.uroute_skipped = "no usable URoute context";
  // The automatic third rung: RSMT may latch to URoute after its own hold.
  const bool learned_ladder = o.routability &&
                              o.route.model == RouteModel::automatic &&
                              o.route.uroute_ladder && o.uroute;
  r.bins_x = o.bins_x;
  r.bins_y = o.bins_y;
  r.initial_hpwl = hpwl(design);
  Design d = design;
  Design model = design;
  int original = count(d.node_x.size());
  const bool fixed_macro_obstacles =
      !d.node_macro.empty() &&
      std::any_of(d.node_macro.begin() + d.num_movable, d.node_macro.end(),
                  [](auto m) { return m != 0; });
  const double fixed_macro_weight =
      fixed_macro_obstacles && !o.macro_transparency.always
          ? o.fixed_macro_overlap_weight
          : 0;
  if (!d.node_boxes.empty()) {
    for (int i = 0; i < d.num_movable; ++i) {
      const auto &b = d.node_boxes[i];
      d.node_x[i] = std::clamp(d.node_x[i], std::max(d.layout.xl, b.xl),
                               std::min(d.layout.xh - d.node_size_x[i], b.xh));
      d.node_y[i] = std::clamp(d.node_y[i], std::max(d.layout.yl, b.yl),
                               std::min(d.layout.yh - d.node_size_y[i], b.yh));
    }
  }
  if (o.observer) {
    o.observer(d, PlacementStage::initial,
               warm ? state->iteration : o.iteration_offset);
  }
  const auto &l = d.layout;
  const double width = l.xh - l.xl, height = l.yh - l.yl;
  double dx = width / o.bins_x, dy = height / o.bins_y;
  // A local engine owns all temporary/filler state. No global RNG or settings.
  std::mt19937 rng(o.seed);
  std::uniform_real_distribution<double> unit(0, 1);
  double area = 0, fixed_area = 0, macro_area = 0, fixed_macro_area = 0;
  for (int i = 0; i < original; ++i) {
    double a = d.node_size_x[i] * d.node_size_y[i];
    if (i < d.num_movable) {
      area += a;
      if (!d.node_macro.empty() && d.node_macro[i]) {
        macro_area += a;
      }
    } else {
      double w = std::max(0., std::min(l.xh, d.node_x[i] + d.node_size_x[i]) -
                                  std::max(l.xl, d.node_x[i]));
      double h = std::max(0., std::min(l.yh, d.node_y[i] + d.node_size_y[i]) -
                                  std::max(l.yl, d.node_y[i]));
      fixed_area += w * h;
      if (!d.node_macro.empty() && d.node_macro[i]) {
        fixed_macro_area += w * h;
      }
    }
  }
  // Reserve physical macro bodies, not their halos or their current overlap
  // union: these bodies must eventually fit even when initial poses overlap.
  double cell_space = width * height;
  if (o.filler_cell_space) {
    for (int i = 0; i < original; ++i) {
      if (!d.node_macro.empty() && d.node_macro[i]) {
        const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
        const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
        cell_space -= std::max(0., d.node_size_x[i] - 2 * px) *
                      std::max(0., d.node_size_y[i] - 2 * py);
      } else if (i >= d.num_movable) {
        cell_space -=
            std::max(0., std::min(l.xh, d.node_x[i] + d.node_size_x[i]) -
                             std::max(l.xl, d.node_x[i])) *
            std::max(0., std::min(l.yh, d.node_y[i] + d.node_size_y[i]) -
                             std::max(l.yl, d.node_y[i]));
      }
    }
    cell_space = std::max(0., cell_space);
  }
  bool reservation_accounting =
      warm ? state->reservation_accounting : !d.density_reservations.empty();
  double filler_rounding_loss = 0;
  const auto gp_start = std::chrono::steady_clock::now();
  if (d.num_movable && o.max_iterations) {
    auto active = [&](int i) { return i < d.num_movable || i >= original; };
    std::vector<int> movable_holds(d.num_movable, 0);
    const auto can_move = [&](int i) {
      return active(i) && (i >= d.num_movable || movable_holds[i] == 0);
    };
    // Area without the halo padding a place() wrapper added to node_size.
    const auto physical_area = [&](int i) {
      const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
      const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
      return (d.node_size_x[i] - 2 * px) * (d.node_size_y[i] - 2 * py);
    };
    auto project = [&](int i, double &x, double &y) {
      const PositionBox box =
          i < original && !d.node_boxes.empty()
              ? d.node_boxes[i]
              : PositionBox{l.xl, l.yl, l.xh - d.node_size_x[i],
                            l.yh - d.node_size_y[i]};
      x = std::clamp(x, std::max(l.xl, box.xl),
                     std::min(l.xh - d.node_size_x[i], box.xh));
      y = std::clamp(y, std::max(l.yl, box.yl),
                     std::min(l.yh - d.node_size_y[i], box.yh));
    };
    for (int i = 0; i < d.num_movable; ++i) {
      if (!warm && !o.iteration_offset && o.macro_jitter > 0 &&
          !d.node_macro.empty() && d.node_macro[i]) {
        // Index-local RNG avoids changing fillers or other nodes when a mask
        // changes.
        std::mt19937 local(o.seed ^ (0x9e3779b9U * (unsigned(i) + 1)));
        std::uniform_real_distribution<double> jitter(-o.macro_jitter,
                                                      o.macro_jitter);
        d.node_x[i] += jitter(local) * d.node_size_x[i];
        d.node_y[i] += jitter(local) * d.node_size_y[i];
      }
      if (!warm && o.local_density_weight > 0 && o.local_density_jitter > 0 &&
          (d.node_macro.empty() || !d.node_macro[i])) {
        std::mt19937 local(o.seed ^ (0x85ebca6bU * (unsigned(i) + 1)));
        std::uniform_real_distribution<double> jitter(-o.local_density_jitter,
                                                      o.local_density_jitter);
        d.node_x[i] += jitter(local) * d.node_size_x[i];
        d.node_y[i] += jitter(local) * d.node_size_y[i];
      }
      project(i, d.node_x[i], d.node_y[i]);
    }
    if (centroids) {
      r.centroid.started = true;
      if (o.observer) {
        o.observer(d, PlacementStage::centroid_start, 0);
      }
    }
    if (o.add_fillers) {
      double target_area =
          std::max(0., o.target_density * (width * height - fixed_area) - area +
                           reservation_charge(d).materialized_area +
                           (o.scale_macro_density
                                ? (1 - o.target_density) * macro_area
                                : 0)) *
          o.filler_area_fraction;
      if (o.filler_cell_space) {
        target_area =
            (1 - o.target_density) * cell_space * o.filler_area_fraction;
      }
      double fw = o.filler_width > 0 ? o.filler_width
                                     : area / d.num_movable / l.row_height;
      double fh = l.row_height;
      if (o.filler_macro_fraction > 0) {
        double smallest_area = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < d.node_x.size(); ++i) {
          if (d.node_macro.empty() || !d.node_macro[i]) {
            continue;
          }
          const double w = d.node_size_x[i] -
                           (d.node_pad_x.empty() ? 0 : 2 * d.node_pad_x[i]);
          const double h = d.node_size_y[i] -
                           (d.node_pad_y.empty() ? 0 : 2 * d.node_pad_y[i]);
          if (w > 0 && h > 0 && w * h < smallest_area) {
            smallest_area = w * h;
            fw = o.filler_macro_fraction * w;
            fh = o.filler_macro_fraction * h;
          }
        }
        require(finite(smallest_area),
                "macro-sized fillers require a physical macro");
      } else if (o.filler_width == 0) {
        // Use the central 90% of real standard-cell dimensions, independently
        // on each axis. SRAM area must not turn virtual particles into long
        // row-height strips, including with no count guard.
        std::vector<double> widths, heights;
        for (int i = 0; i < d.num_movable; ++i) {
          if ((d.node_macro.empty() || !d.node_macro[i]) &&
              d.node_size_x[i] > 0 && d.node_size_y[i] > 0) {
            widths.push_back(d.node_size_x[i]);
            heights.push_back(d.node_size_y[i]);
          }
        }
        const auto trimmed_mean = [](std::vector<double> &v) {
          std::sort(v.begin(), v.end());
          const size_t trim = v.size() / 20;
          return std::accumulate(v.begin() + trim, v.end() - trim, 0.) /
                 (v.size() - 2 * trim);
        };
        if (!widths.empty()) {
          fw = o.filler_cell_scale * trimmed_mean(widths);
          fh = o.filler_cell_scale * trimmed_mean(heights);
        }
      }
      // Macro-stage particles intentionally bypass the cell-resolution bin
      // cap. Their entire footprint remains private density charge.
      const double max_fw = o.filler_macro_fraction > 0
                                ? width
                                : std::min(width, dx * o.filler_bin_fraction);
      const double max_fh = o.filler_macro_fraction > 0
                                ? height
                                : std::min(height, dy * o.filler_bin_fraction);
      fw = std::clamp(fw, l.site_width, std::max(l.site_width, max_fw));
      fh = std::clamp(fh, l.row_height, std::max(l.row_height, max_fh));
      if (o.filler_max_count > 0 &&
          target_area > fw * fh * o.filler_max_count) {
        r.filler_count_limited = true;
        const double scale =
            std::sqrt(target_area / (fw * fh * o.filler_max_count));
        fw = std::min(width, fw * scale);
        fh = std::min(height, fh * scale);
        // Compensate if one dimension hit the outline before the other.
        fw = std::min(width,
                      std::max(fw, target_area / (fh * o.filler_max_count)));
        fh = std::min(height,
                      std::max(fh, target_area / (fw * o.filler_max_count)));
      }
      if (o.magnet_runtime) {
        // Contraction preserves padded macro area, but changes node count.
        // Keep the first solve's particle geometry rather than changing its
        // charge/discretization just because two nodes became one.
        if (o.magnet_runtime->filler_width > 0) {
          fw = o.magnet_runtime->filler_width;
          fh = o.magnet_runtime->filler_height;
        } else {
          o.magnet_runtime->filler_width = fw;
          o.magnet_runtime->filler_height = fh;
        }
      }
      double desired = std::floor(target_area / (fw * fh));
      if (o.filler_max_count > 0) {
        desired = std::min(desired, double(o.filler_max_count));
      }
      require(desired <= double(std::numeric_limits<int>::max() - original),
              "filler index capacity exceeded; configure filler size or an "
              "explicit count guard");
      int fillers = desired;
      r.filler_count = fillers;
      r.filler_area = fillers * fw * fh;
      filler_rounding_loss = std::max(0., target_area - r.filler_area);
      r.filler_size_x = fw;
      r.filler_size_y = fh;
      for (int i = 0; i < fillers; ++i) {
        d.node_x.push_back(l.xl + unit(rng) * (width - fw));
        d.node_y.push_back(l.yl + unit(rng) * (height - fh));
        d.node_size_x.push_back(fw);
        d.node_size_y.push_back(fh);
      }
      d.build_maps();
    }
    int n = count(d.node_x.size());
    if (warm) {
      for (const auto *a :
           {&state->x, &state->y, &state->prev_x, &state->prev_y,
            &state->old_gx, &state->old_gy, &state->old_ex, &state->old_ey,
            &state->inflation}) {
        require(a->size() == static_cast<size_t>(n) &&
                    std::all_of(a->begin(), a->end(),
                                [](double v) { return finite(v); }),
                "invalid GP continuation arrays");
      }
      require(state->iteration >= 0 && finite(state->density_weight) &&
                  state->density_weight > 0 && finite(state->momentum) &&
                  state->momentum >= 1 && finite(state->step) &&
                  state->step >= 0,
              "invalid GP continuation scalars");
    }
    std::vector<double> inflation(n, 1);
    if (warm) {
      inflation = state->inflation;
    }
    // Refresh only when weights/connectivity change, not on every GP step.
    // The maximum protects a critical branch even if other incident nets have
    // slack. Masked nets, macros, fixed nodes and virtual fillers do not enter.
    const auto make_spread_mobility = [&] {
      std::vector<double> mobility;
      r.timing_spread_cells = 0;
      r.timing_spread_min_mobility = 1;
      if (o.timing_spread_weight == 0) {
        return mobility;
      }
      mobility.assign(n, 1);
      for (int i = 0; i < d.num_movable; ++i) {
        if (!d.node_macro.empty() && d.node_macro[i]) {
          continue;
        }
        double weight = 1;
        for (int k = d.flat_node2pin_start[i]; k < d.flat_node2pin_start[i + 1];
             ++k) {
          const int net = d.pin2net[d.flat_node2pin[k]];
          if (d.net_mask[net]) {
            weight = std::max(weight, d.net_weights[net]);
          }
        }
        mobility[i] =
            1 - o.timing_spread_weight + o.timing_spread_weight / weight;
        if (mobility[i] < 1) {
          ++r.timing_spread_cells;
          r.timing_spread_min_mobility =
              std::min(r.timing_spread_min_mobility, mobility[i]);
        }
      }
      return mobility;
    };
    auto spread_mobility = make_spread_mobility();
    std::vector<double> sx(n), sy(n), ox(n), oy(n), ratio(n);
    const auto sharp = [&](int i) {
      return (i < d.num_movable && o.exact_cell_density &&
              (d.node_macro.empty() || !d.node_macro[i])) ||
             (i < original && o.exact_macro_density && !d.node_macro.empty() &&
              d.node_macro[i]);
    };
    for (int i = 0; i < n; ++i) {
      sx[i] = active(i) && !sharp(i)
                  ? std::max(d.node_size_x[i], std::sqrt(2.) * dx)
                  : d.node_size_x[i];
      sy[i] = active(i) && !sharp(i)
                  ? std::max(d.node_size_y[i], std::sqrt(2.) * dy)
                  : d.node_size_y[i];
      ox[i] = (d.node_size_x[i] - sx[i]) / 2;
      oy[i] = (d.node_size_y[i] - sy[i]) / 2;
      ratio[i] = sx[i] * sy[i] != 0
                     ? d.node_size_x[i] * d.node_size_y[i] / (sx[i] * sy[i])
                     : 0;
      if (!active(i) || (o.scale_macro_density && i < original &&
                         !d.node_macro.empty() && d.node_macro[i])) {
        ratio[i] *= o.target_density;
      }
    }
    PoissonSolver solver(o.bins_x, o.bins_y, width, height, o.dynamic_bins, 0,
                         o.num_threads);
    std::unique_ptr<detail::LocalDensity> local_density;
    std::unique_ptr<detail::LocalDensity> macro_cell_density;
    double macro_cell_scale = -1;
    const auto rebuild_macro_cell_density = [&] {
      if (o.macro_cell_density_weight <= 0) {
        return;
      }
      const double reach = detail::macro_cell_screening_length(
          model, o.macro_cell_density_reach);
      r.macro_cell_density_screening_length = reach;
      const double pitch = std::sqrt(width * height) / o.local_density_bins;
      macro_cell_density = std::make_unique<detail::LocalDensity>(
          model, o.local_density_bins, std::max(reach / pitch, 1e-12),
          o.target_density, o.num_threads);
      macro_cell_scale = -1;
    };
    rebuild_macro_cell_density();

    const double local_target =
        o.local_density_target < 0 ? o.target_density : o.local_density_target;
    if (o.local_density_weight > 0 || o.macro_cell_density_weight > 0 ||
        o.macro_density_pull_weight > 0) {
      require(!state || o.trial_tuning,
              "local spreading continuation requires a GP trial");
      local_density = std::make_unique<detail::LocalDensity>(
          model, o.local_density_bins, o.local_density_range, local_target,
          o.num_threads);
    }
    std::vector<double> map, phi, fx, fy, reserved_map;
    bool reservation_map_dirty = true;
    auto x = d.node_x, y = d.node_y, prev_x = x, prev_y = y;
    std::vector<double> old_gx(n), old_gy(n), old_ex(n), old_ey(n);
    double density_weight = o.initial_density_weight, t = 1, alpha = 0;
    bool density_initialized =
        warm ? state->density_initialized
             : (o.magnet_runtime && o.magnet_runtime->density_initialized);
    const int offset = warm ? state->iteration : o.iteration_offset;
    double previous_overflow = warm ? state->previous_overflow : -1;
    int macro_pull_changed =
        o.magnet_runtime ? o.magnet_runtime->macro_pull_changed : -1;
    double previous_hpwl = warm ? state->previous_hpwl : r.initial_hpwl;
    double route_overflow = 0, inflated_ratio = 1;
    const auto *rigid_route = o.magnet_runtime && o.magnet_runtime->routing
                                  ? &*o.magnet_runtime->routing
                                  : nullptr;
    const bool cycling = o.macro_transparency.cycles > 0;
    const auto transparency_end =
        std::int64_t(o.macro_transparency.warmup) +
        2 * std::int64_t(o.macro_transparency.period) *
            o.macro_transparency.cycles;
    bool transparent = false;
    double exclusion_scale = 1;
    RouteModel active_route =
        warm ? state->route_model
        : rigid_route
            ? rigid_route->model
            : (o.route.model == RouteModel::uroute
                   ? (o.uroute ? RouteModel::uroute : RouteModel::rsmt)
                   : (o.route.model == RouteModel::rsmt ? RouteModel::rsmt
                                                        : RouteModel::rudy));
    if (!o.routability) {
      active_route = RouteModel::off;
    }
    int route_hold = warm ? state->route_hold : 0;
    int uroute_hold = warm ? state->route_uroute_hold : 0;
    int probe_deadline = warm          ? state->route_probe_deadline
                         : rigid_route ? rigid_route->probe_deadline
                                       : 0;
    if (warm) {
      r.route_switch_iteration = state->route_switch_iteration;
      r.route_switch_rudy = state->route_switch_rudy;
      r.route_switch_rsmt = state->route_switch_rsmt;
      r.route_probed = state->route_probed;
      r.route_uroute_switch_iteration = state->route_uroute_switch_iteration;
      r.route_switch_uroute_rsmt = state->route_switch_uroute_rsmt;
      r.route_switch_uroute = state->route_switch_uroute;
      r.route_probe_overflow = state->route_probe_overflow;
    } else if (rigid_route) {
      r.route_switch_iteration = rigid_route->switch_iteration;
      r.route_switch_rudy = rigid_route->switch_rudy;
      r.route_switch_rsmt = rigid_route->switch_rsmt;
      r.route_probed = rigid_route->probed;
      r.route_uroute_switch_iteration = rigid_route->uroute_switch_iteration;
      r.route_switch_uroute_rsmt = rigid_route->switch_uroute_rsmt;
      r.route_switch_uroute = rigid_route->switch_uroute;
      r.route_probe_overflow = rigid_route->probe_overflow;
    }
    double route_sample_seconds = 0, route_peak_utilization = 0;
    bool route_sampled = false;
    int route_blocked_pins = 0, route_unroutable_nets = 0, route_calls = 0;
    double route_rudy_seconds = 0, route_rsmt_seconds = 0,
           route_uroute_seconds = 0;
    r.route_model = active_route;
    auto measure_route = [&](const Design &physical, RouteModel model,
                             bool record_history = true) {
      const auto started = std::chrono::steady_clock::now();
      RudyOptions route_options{aux_x,
                                aux_y,
                                o.route_capacity_h,
                                o.route_capacity_v,
                                o.improved_routability,
                                o.improved_routability};
      route_options.model = model == RouteModel::uroute ? RudyModel::uroute
                            : model == RouteModel::rsmt ? RudyModel::rsmt
                                                        : RudyModel::bbox;
      route_options.uroute = o.uroute;
      route_options.degree = o.route.degree;
      route_options.accuracy = o.route.accuracy;
      route_options.initialize = o.route_initialize;
      route_options.fixed_macro_scale = exclusion_scale;
      if (transparent) {
        route_options.macro_blockages = false;
        route_options.pin_access = false;
      }
      std::optional<Design> expanded;
      if (o.route_physical_snapshot)
        expanded = o.route_physical_snapshot(physical);
      auto result = rudy(expanded ? *expanded : physical, route_options);
      const double seconds = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - started)
                                 .count();
      r.rudy_seconds += seconds;
      ++r.rudy_calls;
      // Trial scores have their own trace. Do not let a losing candidate
      // overwrite the selected trajectory's most recent routing sample.
      if (record_history) {
        route_sample_seconds += seconds;
        ++route_calls;
        route_sampled = true;
        if (model == RouteModel::uroute) {
          route_uroute_seconds += seconds;
        } else if (model == RouteModel::rsmt) {
          route_rsmt_seconds += seconds;
        } else {
          route_rudy_seconds += seconds;
        }
        route_blocked_pins = result.blocked_pins;
        route_unroutable_nets = result.unroutable_nets;
      }
      return result;
    };
    auto latch_route = [&](const Design &physical, RudyResult old,
                           int iteration) {
      auto result = measure_route(physical, RouteModel::rsmt);
      r.route_switch_iteration = iteration;
      r.route_switch_rudy = old.overflow;
      r.route_switch_rsmt = result.overflow;
      active_route = r.route_model = RouteModel::rsmt;
      return result;
    };
    auto estimate_route = [&](const Design &physical, int iteration) {
      auto result = measure_route(physical, active_route);
      if (active_route == RouteModel::rudy &&
          o.route.model == RouteModel::automatic) {
        route_hold =
            result.overflow <= o.route.switch_overflow ? route_hold + 1 : 0;
        if (route_hold >= o.route.hold) {
          result = latch_route(physical, std::move(result), iteration);
        }
      }
      // Do not count the first-rung switch sample toward the second hold.
      // Consecutive samples are measured on RSMT before a one-way URoute latch.
      else if (active_route == RouteModel::rsmt && learned_ladder) {
        uroute_hold = result.overflow <= o.route.uroute_switch_overflow
                          ? uroute_hold + 1
                          : 0;
        if (uroute_hold >= o.route.uroute_hold) {
          r.route_switch_uroute_rsmt = result.overflow;
          result = measure_route(physical, RouteModel::uroute);
          r.route_switch_uroute = result.overflow;
          r.route_uroute_switch_iteration = iteration;
          active_route = r.route_model = RouteModel::uroute;
        }
      }
      route_peak_utilization = result.max_utilization;
      return result;
    };
    auto route_history = [&](GpSample &sample) {
      sample.route_blocked_pins = route_blocked_pins;
      sample.route_unroutable_nets = route_unroutable_nets;
      sample.route_peak_utilization = route_peak_utilization;
      sample.route_calls = route_calls;
      sample.route_rudy_seconds = route_rudy_seconds;
      sample.route_rsmt_seconds = route_rsmt_seconds;
      sample.route_uroute_seconds = route_uroute_seconds;
      sample.route_uroute_switch_iteration = r.route_uroute_switch_iteration;
      sample.route_switch_uroute_rsmt = r.route_switch_uroute_rsmt;
      sample.route_switch_uroute = r.route_switch_uroute;
      sample.route_probe_overflow = r.route_probe_overflow;
      route_calls = 0;
      route_rudy_seconds = route_rsmt_seconds = route_uroute_seconds = 0;
      sample.route_model = active_route;
      sample.route_sampled = route_sampled;
      sample.route_seconds = route_sample_seconds;
      sample.route_switch_iteration = r.route_switch_iteration;
      sample.route_switch_rudy = r.route_switch_rudy;
      sample.route_switch_rsmt = r.route_switch_rsmt;
      sample.route_probe = r.route_probed;
      sample.route_bins = aux_x * aux_y;
      route_sample_seconds = 0;
      route_sampled = false;
    };

    const bool macro_routing =
        o.routability && o.improved_routability &&
        o.routability_macro_weight > 0 && !d.node_macro.empty() &&
        std::any_of(d.node_macro.begin(), d.node_macro.begin() + d.num_movable,
                    [](auto macro) { return macro != 0; });
    auto routing_utilization =
        warm ? state->routing_utilization : std::vector<double>{};
    require(
        routing_utilization.empty() ||
            (routing_utilization.size() == size_t(aux_x * aux_y) &&
             std::all_of(routing_utilization.begin(), routing_utilization.end(),
                         [](double u) { return finite(u) && u >= 0; })),
        "invalid GP routing continuation");
    double macro_routing_force = 0, macro_routing_force_peak = 0;
    std::uint64_t macro_routing_channels = 0, macro_routing_channels_peak = 0;
    // With a soft upper bound, density pressure may overshoot a better feasible
    // placement while waiting for the minimum iteration window. Preserve that
    // placement for a static solve, without changing continuation/checkpoint
    // state or any of the classic density path.
    const bool keep_feasible =
        o.density_upper_bound && !state && !o.magnet_runtime && !o.checkpoint &&
        !centroids && !(o.order_weight > 0 && !design.order_edges.empty()) &&
        !o.macro_align_weight && !macro_routing &&
        !(o.macro_pull_update && o.macro_pull_interval > 0);
    double best_feasible_hpwl = std::numeric_limits<double>::infinity();
    bool pin_timing_seen = d.pin_timing && d.pin_timing->enabled();
    std::vector<MacroPull> macro_pulls = o.magnet_runtime
                                             ? o.magnet_runtime->macro_pulls
                                             : std::vector<MacroPull>{};
    std::vector<double> best_feasible_x, best_feasible_y;
    int best_feasible_iteration = 0;
    bool restart = o.trial_tuning.has_value();
    int last_checkpoint = -1;
    bool movable_updates_active = o.movable_update && o.movable_interval > 0;
    bool checkpoint_invoked =
        o.magnet_runtime && o.magnet_runtime->checkpoint_invoked;
    const bool density_gated_checkpoints =
        o.checkpoint_after_density && o.checkpoint && o.checkpoint_interval > 0;
    auto convergence = warm ? state->convergence
                            : (o.magnet_runtime && o.magnet_runtime->convergence
                                   ? *o.magnet_runtime->convergence
                                   : GpConvergence{});
    if (!warm && o.magnet_runtime) {
      // Current poses seed a cold solver after contraction/snap. Objective
      // samples from the previous geometry cannot establish convergence.
      convergence.previous_hpwl = convergence.previous_inflation = -1;
      convergence.stable_samples = 0;
    }
    std::vector<AlignmentPair> align_pairs =
        warm ? state->align_pairs : std::vector<AlignmentPair>{};
    double align_residual = 0, align_ramp = 0, align_force_peak = 0;
    const auto alignment_pairs = [&](const Design &view) {
      return o.macro_align_global
                 ? detail::macro_grid_align_pairs(view, o.macro_align_nearest)
                 : detail::macro_align_pairs(view, o.macro_align_nearest);
    };
    double centroid_repulsion_weight = o.centroid_repulsion_weight;
    std::uint64_t centroid_parent_edges = 0;
    double centroid_parent_distance = 0, centroid_parent_force = 0;
    std::uint64_t centroid_align_pairs = 0;
    double centroid_align_residual = 0, centroid_align_ramp = 0,
           centroid_align_force = 0;
    auto centroid_alignment_schedule = o;
    centroid_alignment_schedule.macro_align_start = o.centroid_align_start;
    centroid_alignment_schedule.macro_align_full = o.centroid_align_full;
    double previous_circle_overlap = 1;
    if (warm && o.improved_routability) {
      double inflated = 0;
      for (int i = 0; i < d.num_movable; ++i) {
        inflated += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
      }
      inflated_ratio = area ? inflated / area : 1;
    }
    if (warm) {
      for (const auto &p : state->align_pairs) {
        require(
            p.a >= 0 && p.b >= 0 && p.a < d.num_movable &&
                p.b < d.num_movable && p.a != p.b &&
                (p.axis == Axis::x || p.axis == Axis::y) &&
                !d.node_align_class.empty() && d.node_align_class[p.a] >= 0 &&
                d.node_align_class[p.a] == d.node_align_class[p.b] &&
                !d.node_macro.empty() && d.node_macro[p.a] && d.node_macro[p.b],
            "invalid GP continuation alignment pairs");
      }
      x = state->x;
      y = state->y;
      prev_x = state->prev_x;
      prev_y = state->prev_y;
      old_gx = state->old_gx;
      old_gy = state->old_gy;
      old_ex = state->old_ex;
      old_ey = state->old_ey;
      density_weight = state->density_weight;
      t = state->momentum;
      alpha = state->step;
      // External compaction may have moved physical nodes since the checkpoint.
      for (int i = 0; i < original; ++i) {
        if (i >= d.num_movable) {
          // Refreshed fixed pseudo-terminals must retain the caller's exact
          // coordinates. Translating old state by (new - old) can lose an ulp.
          x[i] = prev_x[i] = old_ex[i] = d.node_x[i];
          y[i] = prev_y[i] = old_ey[i] = d.node_y[i];
          continue;
        }
        const double dx = d.node_x[i] - x[i], dy = d.node_y[i] - y[i];
        x[i] += dx;
        y[i] += dy;
        prev_x[i] += dx;
        prev_y[i] += dy;
        old_ex[i] += dx;
        old_ey[i] += dy;
      }
    }
    // Telemetry has its own lifetime and never enables magnet force/stopping
    // rules. Rigid placement uses its existing shared runtime instead.
    std::optional<detail::MagnetRuntime> step_recorder;
    if (o.macro_step_observer && !o.magnet_runtime) {
      step_recorder.emplace(o);
    }
    if (step_recorder) {
      step_recorder->observe(model, x, y, offset);
    }
    if (o.magnet_runtime) {
      o.magnet_runtime->observe(model, x, y, offset);
    }
    if (o.macro_motion_observer) {
      o.macro_motion_observer(model, x, y, offset);
    }
    const auto save_state = [&](int iteration) {
      GpState s;
      s.pin_timing = d.pin_timing;
      s.density_reservations = d.density_reservations;
      s.reservation_accounting = reservation_accounting;
      s.convergence = convergence;
      s.signature = signature;
      s.iteration = iteration;
      s.x = x;
      s.y = y;
      s.prev_x = prev_x;
      s.prev_y = prev_y;
      s.old_gx = old_gx;
      s.old_gy = old_gy;
      s.old_ex = old_ex;
      s.old_ey = old_ey;
      s.density_weight = density_weight;
      s.density_initialized = density_initialized;
      s.previous_overflow = previous_overflow;
      s.previous_hpwl = previous_hpwl;
      s.momentum = t;
      s.step = alpha;
      s.inflation = inflation;
      s.routing_utilization = routing_utilization;
      s.route_model = active_route;
      s.route_hold = route_hold;
      s.route_uroute_hold = uroute_hold;
      s.route_probe_deadline = probe_deadline;
      s.route_probed = r.route_probed;
      s.route_switch_iteration = r.route_switch_iteration;
      s.route_switch_rudy = r.route_switch_rudy;
      s.route_switch_rsmt = r.route_switch_rsmt;
      s.route_uroute_switch_iteration = r.route_uroute_switch_iteration;
      s.route_switch_uroute_rsmt = r.route_switch_uroute_rsmt;
      s.route_switch_uroute = r.route_switch_uroute;
      s.route_probe_overflow = r.route_probe_overflow;
      s.align_pairs = align_pairs;
      return s;
    };
    struct RetainedCheckpoint {
      Design design, model;
      GpState state;
      bool restart = false, pin_timing_seen = false;
      double overflow = 0;
      std::vector<double> sx, sy, ox, oy, ratio;
      double effective_filler_area = 0;
    };
    std::optional<RetainedCheckpoint> retained_checkpoint;
    int restored_selection = 0;
    GpTrialTuning tuning = o.trial_tuning.value_or(
        GpTrialTuning{1, o.vanished_memory ? o.vanished_memory_local_floor
                                           : o.local_density_weight});
    int search_round = 0;
    bool tuned = o.trial_tuning.has_value();
    for (int iter = 0; iter < o.max_iterations; ++iter) {
      const int global_iter = offset + iter;
      bool retain_checkpoint = false;
      restored_selection = 0;
      if (o.gp_search.enabled && global_iter >= o.checkpoint_start &&
          global_iter >= o.gp_search.start &&
          (global_iter - o.gp_search.start) % o.gp_search.interval == 0 &&
          (o.gp_search.rounds == 0 || search_round < o.gp_search.rounds) &&
          iter + o.gp_search.steps + 20 <= o.max_iterations) {
        Design parent = model;
        std::copy_n(x.begin(), original, parent.node_x.begin());
        std::copy_n(y.begin(), original, parent.node_y.begin());
        GpTrialSample sample;
        sample.iteration = global_iter;
        sample.steps = o.gp_search.steps;
        sample.round = search_round++;
        sample.adjusting_fillers = sample.round % 2 == 0;
        sample.settings.push_back(tuning); // unchanged control is always first
        for (int direction : {-1, 1}) {
          auto candidate = tuning;
          double &value = sample.adjusting_fillers ? candidate.filler_charge
                                                   : candidate.local_weight;
          value = std::clamp(value + direction * (sample.adjusting_fillers
                                                      ? o.gp_search.filler_step
                                                      : o.gp_search.local_step),
                             0., 1.);
          if (candidate.filler_charge != tuning.filler_charge ||
              candidate.local_weight != tuning.local_weight) {
            sample.settings.push_back(candidate);
          }
        }
        const auto score = [&](const Design &physical) {
          GpTrialScore q;
          const auto route = measure_route(physical, RouteModel::rsmt, false);
          q.route_overflow = route.overflow;
          q.route_unroutable_nets = route.unroutable_nets;
          q.route_blocked_pins = route.blocked_pins;
          q.route_measured_nets = route.measured_nets;
          q.route_congestion_overflow =
              std::max(0., route.overflow - double(route.unroutable_nets +
                                                   route.blocked_pins) /
                                                (aux_x * aux_y));
          q.route_peak_utilization = route.max_utilization;
          q.route_detour_wirelength = route.detour_wirelength;
          const auto [wns, tns] = o.trial_timing(physical);
          q.wns = wns;
          q.tns = tns;
          auto unweighted = physical;
          std::fill(unweighted.net_weights.begin(),
                    unweighted.net_weights.end(), 1.);
          q.hpwl = hpwl(unweighted);
          q.physical_overflow =
              overflow(physical, o.bins_x, o.bins_y, o.target_density,
                       o.scale_macro_density);
          return q;
        };
        sample.before = score(parent);
        auto checkpoint = save_state(global_iter);
        // Every candidate, including unchanged control, starts with the same
        // objective-boundary momentum restart. Retain filler positions, density
        // normalization, inflation, weights and the global smoothing schedule.
        checkpoint.momentum = 1;
        checkpoint.step = 0;
        checkpoint.prev_x = checkpoint.x;
        checkpoint.prev_y = checkpoint.y;
        checkpoint.old_ex = checkpoint.x;
        checkpoint.old_ey = checkpoint.y;
        std::fill(checkpoint.old_gx.begin(), checkpoint.old_gx.end(), 0.);
        std::fill(checkpoint.old_gy.begin(), checkpoint.old_gy.end(), 0.);
        struct Trial {
          Design physical;
          GpState state;
          Result result;
        };
        std::vector<std::future<Trial>> futures;
        // Small thread budgets must not launch one worker per candidate in
        // addition to each candidate's own native worker team.
        const bool parallel_trials =
            o.num_threads >= int(sample.settings.size());
        for (const auto &setting : sample.settings) {
          auto child = o;
          child.gp_search.enabled = false;
          child.trial_tuning = setting;
          child.max_iterations = sample.steps;
          child.legalize = child.detailed_place = false;
          child.before_legalize = {};
          child.after_macro_legalize = {};
          child.checkpoint = {};
          child.checkpoint_interval = 0;
          child.checkpoint_observer = {};
          child.observer = {};
          child.progress_observer = {};
          child.trial_timing = {};
          child.trial_observer = {};
          child.num_threads = parallel_trials
                                  ? o.num_threads / int(sample.settings.size())
                                  : o.num_threads;
          futures.push_back(std::async(
              parallel_trials ? std::launch::async : std::launch::deferred,
              [parent, checkpoint, child]() mutable {
                Trial trial{std::move(parent), std::move(checkpoint), {}};
                trial.result = place(trial.physical, child, &trial.state);
                return trial;
              }));
        }
        std::vector<Trial> trials;
        for (auto &f : futures) {
          trials.push_back(f.get());
        }
        // STA/host graph access is serialized; losing trials never refresh the
        // shared timing weights, observer stream or placement state.
        for (const auto &trial : trials) {
          require(trial.state.iteration == global_iter + sample.steps,
                  "unequal GP trial horizon");
          sample.scores.push_back(score(trial.physical));
        }
        sample.selected =
            select_gp_trial(sample.scores, o.gp_search, &sample.reason);
        if (o.trial_observer) {
          o.trial_observer(sample);
        }
        auto &winner = trials[sample.selected];
        tuning = sample.settings[sample.selected];
        tuned = true;
        auto &s = winner.state;
        x = std::move(s.x);
        y = std::move(s.y);
        prev_x = std::move(s.prev_x);
        prev_y = std::move(s.prev_y);
        old_gx = std::move(s.old_gx);
        old_gy = std::move(s.old_gy);
        old_ex = std::move(s.old_ex);
        old_ey = std::move(s.old_ey);
        density_weight = s.density_weight;
        density_initialized = s.density_initialized;
        previous_overflow = s.previous_overflow;
        previous_hpwl = s.previous_hpwl;
        t = s.momentum;
        alpha = s.step;
        convergence = s.convergence;
        inflation = std::move(s.inflation);
        routing_utilization = std::move(s.routing_utilization);
        active_route = s.route_model;
        route_hold = s.route_hold;
        uroute_hold = s.route_uroute_hold;
        probe_deadline = s.route_probe_deadline;
        r.route_probed = s.route_probed;
        r.route_switch_iteration = s.route_switch_iteration;
        r.route_switch_rudy = s.route_switch_rudy;
        r.route_switch_rsmt = s.route_switch_rsmt;
        r.route_uroute_switch_iteration = s.route_uroute_switch_iteration;
        r.route_switch_uroute_rsmt = s.route_switch_uroute_rsmt;
        r.route_switch_uroute = s.route_switch_uroute;
        r.route_probe_overflow = s.route_probe_overflow;
        align_pairs = std::move(s.align_pairs);
        model = std::move(winner.physical);
        double inflated = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          inflated += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
        }
        inflated_ratio = area ? inflated / area : 1;
        route_overflow = sample.scores[sample.selected].route_overflow;
        r.overflow = sample.scores[sample.selected].physical_overflow;
        restart = false;
        last_checkpoint = global_iter + sample.steps;
        convergence.previous_hpwl = -1;
        convergence.previous_inflation = -1;
        convergence.stable_samples = 0;
        best_feasible_x.clear();
        best_feasible_y.clear();
        best_feasible_hpwl = std::numeric_limits<double>::infinity();
        for (const auto &h : winner.result.history) {
          r.history.push_back(h);
          if (o.progress_observer) {
            o.progress_observer(h);
          }
        }
        iter += sample.steps - 1;
        r.iterations = iter + 1;
        if (o.observer) {
          o.observer(model, PlacementStage::global_placement,
                     offset + r.iterations);
        }
        continue;
      }
      const double next_scale =
          o.vanished_memory
              ? vanished_memory_scale(global_iter + 1,
                                      o.vanished_memory_initial_scale,
                                      o.vanished_memory_growth_iterations,
                                      o.vanished_memory_start)
              : 1;
      const double local_weight =
          tuned ? tuning.local_weight
          : o.vanished_memory
              ? vanished_memory_local_weight(
                    o.local_density_weight, global_iter + 1,
                    o.vanished_memory_local_floor, o.vanished_memory_start)
              : o.local_density_weight;
      const bool next_transparent =
          o.vanished_memory
              ? next_scale == 0
              : macro_transparent(o.macro_transparency, global_iter);
      const bool schedule_boundary =
          o.vanished_memory &&
          (global_iter + 1 ==
               o.vanished_memory_start + o.vanished_memory_growth_iterations ||
           global_iter + 1 == o.vanished_memory_start + 1 ||
           (global_iter >= o.vanished_memory_start + 1 &&
            local_weight == o.vanished_memory_local_floor &&
            vanished_memory_local_weight(o.local_density_weight, global_iter,
                                         o.vanished_memory_local_floor,
                                         o.vanished_memory_start) >
                o.vanished_memory_local_floor));
      exclusion_scale = next_scale;
      if (next_transparent != transparent || schedule_boundary) {
        transparent = next_transparent;
        t = 1;
        prev_x = x;
        prev_y = y;
        restart = true;
        last_checkpoint = global_iter;
        convergence.previous_hpwl = -1;
        convergence.previous_inflation = -1;
        convergence.stable_samples = 0;
        previous_hpwl = -1;
        best_feasible_x.clear();
        best_feasible_y.clear();
        best_feasible_hpwl = std::numeric_limits<double>::infinity();
        // Inflation from a blocked routing model must not persist through a
        // transparent phase (or vice versa). Re-estimate using the new model.
        std::fill(inflation.begin(), inflation.end(), 1);
        inflated_ratio = 1;
        routing_utilization.clear();
      }
      const bool opacity_settled =
          o.vanished_memory ? global_iter + 1 >= memory_settled_iteration(o)
                            : (!cycling || global_iter + 1 >=
                                               transparency_end +
                                                   o.macro_transparency.period);
      double progress = std::min(
          1., double(global_iter) /
                  std::max(1, (o.schedule_iterations ? o.schedule_iterations
                                                     : o.max_iterations) -
                                  1));
      double gamma = std::max(l.site_width,
                              4 * std::max(dx, dy) *
                                  std::pow(o.wirelength_gamma_end, progress));
      if (o.wirelength_gamma_start != 1) {
        gamma = std::max(
            l.site_width,
            4 * std::max(dx, dy) * o.wirelength_gamma_start *
                std::pow(o.wirelength_gamma_end / o.wirelength_gamma_start,
                         progress));
      }
      double nt = (1 + std::sqrt(1 + 4 * t * t)) / 2, momentum = (t - 1) / nt;
      for (int i = 0; i < n; ++i) {
        if (can_move(i)) {
          d.node_x[i] = x[i] + momentum * (x[i] - prev_x[i]);
          d.node_y[i] = y[i] + momentum * (y[i] - prev_y[i]);
          project(i, d.node_x[i], d.node_y[i]);
        }
      }
      auto wg = detail::wirelength(d, gamma, o.num_threads);
      double pin_timing_cost = 0;
      if (d.pin_timing && d.pin_timing->enabled()) {
        const auto timing = detail::pin_timing_force(d, gamma);
        pin_timing_cost = timing.value;
        wg.value += timing.value;
        for (int i = 0; i < n; ++i) {
          wg.x[i] += timing.x[i];
          wg.y[i] += timing.y[i];
        }
      }
      const double lambda =
          o.order_weight * std::pow(o.order_growth, global_iter);
      require(finite(lambda), "relative-order schedule overflow");
      if (lambda) {
        for (const auto &e : d.order_edges) {
          const bool horizontal = e.axis == Axis::x;
          const auto &pos = horizontal ? d.node_x : d.node_y;
          const auto &size = horizontal ? d.node_size_x : d.node_size_y;
          const double violation = std::max(0., pos[e.from] + size[e.from] +
                                                    e.separation - pos[e.to]);
          auto &g = horizontal ? wg.x : wg.y;
          const double force = 2 * lambda * e.strength.weight * violation;
          g[e.from] += force;
          g[e.to] -= force;
        }
      }
      if (o.routability &&
          (global_iter % 20 == 0 ||
           (o.vanished_memory && global_iter >= o.vanished_memory_start &&
            global_iter < o.vanished_memory_start +
                              o.vanished_memory_growth_iterations)) &&
          (!o.improved_routability ||
           (previous_overflow >= 0 && previous_overflow < .3))) {
        Design physical = model;
        std::copy_n(d.node_x.begin(), original, physical.node_x.begin());
        std::copy_n(d.node_y.begin(), original, physical.node_y.begin());
        const auto route = estimate_route(physical, global_iter + 1);
        route_overflow = route.overflow;
        if (macro_routing) {
          // Learned utilization only drives inflation and convergence. Macro
          // channel forces retain their physical RSMT demand measurement.
          std::optional<RudyResult> rsmt_pressure;
          if (active_route == RouteModel::uroute)
            rsmt_pressure = measure_route(physical, RouteModel::rsmt, false);
          const auto &pressure = rsmt_pressure ? *rsmt_pressure : route;
          routing_utilization.resize(pressure.horizontal.size());
          for (size_t k = 0; k < routing_utilization.size(); ++k) {
            routing_utilization[k] =
                std::max({pressure.horizontal[k], pressure.vertical[k],
                          pressure.pin_utilization[k]});
          }
        }
        auto proposed = inflation;
        double added = 0, current = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          current += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
        }
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_inflatable.empty() && d.node_inflatable[i]) {
            const int bx =
                std::clamp(int((d.node_x[i] + d.node_size_x[i] / 2 - l.xl) /
                               (width / aux_x)),
                           0, aux_x - 1);
            const int by =
                std::clamp(int((d.node_y[i] + d.node_size_y[i] / 2 - l.yl) /
                               (height / aux_y)),
                           0, aux_y - 1);
            const double demand =
                std::clamp(std::max({route.pin_utilization[bx * aux_y + by],
                                     route.horizontal[bx * aux_y + by],
                                     route.vertical[bx * aux_y + by]}),
                           1., 3.);
            proposed[i] = o.improved_routability
                              ? std::max(inflation[i], demand)
                              : demand;
            added += (proposed[i] - inflation[i]) * d.node_size_x[i] *
                     d.node_size_y[i];
          }
        }
        const double budget = o.improved_routability && added > 0
                                  ? std::min(1., .05 * current / added)
                                  : 1.;
        double inflated = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          inflation[i] += budget * (proposed[i] - inflation[i]);
          inflated += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
        }
        inflated_ratio = area ? inflated / area : 1;
      }
      const double target =
          o.target_density / (o.improved_routability ? inflated_ratio : 1);
      double filler_scale = 1;
      double excluded_macro_area = transparent ? 0 : fixed_macro_area;
      if (o.vanished_memory && exclusion_scale > 0 && exclusion_scale < 1) {
        excluded_macro_area = 0;
        for (int i = d.num_movable; i < original; ++i) {
          if (d.node_macro.empty() || !d.node_macro[i]) {
            continue;
          }
          const double w = d.node_size_x[i] * exclusion_scale;
          const double h = d.node_size_y[i] * exclusion_scale;
          const double mx = d.node_x[i] + (d.node_size_x[i] - w) / 2;
          const double my = d.node_y[i] + (d.node_size_y[i] - h) / 2;
          excluded_macro_area +=
              std::max(0., std::min(l.xh, mx + w) - std::max(l.xl, mx)) *
              std::max(0., std::min(l.yh, my + h) - std::max(l.yl, my));
        }
      }
      if ((o.gp_search.enabled || o.trial_tuning || o.vanished_memory ||
           (o.improved_routability && o.routability_charge_compensation) ||
           (o.macro_transparency.always ||
            (cycling && global_iter >= o.macro_transparency.warmup))) &&
          n > original) {
        double filler_area = 0;
        for (int i = original; i < n; ++i) {
          filler_area += d.node_size_x[i] * d.node_size_y[i];
        }
        // Keep total charge consistent with the effective capacity. A uniform
        // target shift is the discarded DC mode and cannot change forces.
        const double real_charge =
            area * inflated_ratio -
            (o.scale_macro_density ? (1 - target) * macro_area : 0);
        const double freed_macro_area =
            o.vanished_memory ? fixed_macro_area - excluded_macro_area
                              : (transparent ? fixed_macro_area : 0);
        const double budget =
            o.filler_cell_space
                ? (1 - o.target_density) * cell_space * o.filler_area_fraction
                : std::max(0., target * (width * height - fixed_area +
                                         freed_macro_area) -
                                   real_charge) *
                      o.filler_area_fraction;
        filler_scale = filler_area ? budget / filler_area : 0;
        for (int i = original; i < n; ++i) {
          ratio[i] = filler_scale * d.node_size_x[i] * d.node_size_y[i] /
                     (sx[i] * sy[i]);
        }
        for (int i = d.num_movable; i < original; ++i) {
          ratio[i] = sx[i] * sy[i] > 0 ? target * d.node_size_x[i] *
                                             d.node_size_y[i] / (sx[i] * sy[i])
                                       : 0;
        }
      }
      if (tuned) {
        filler_scale *= tuning.filler_charge;
        for (int i = original; i < n; ++i) {
          ratio[i] *= tuning.filler_charge;
        }
      }
      const auto apply_reservation_fillers = [&](double current_target) {
        require(o.add_fillers,
                "density reservations require enabled filler accounting");
        const auto charge = reservation_charge(d, inflation);
        const double real_charge =
            area * inflated_ratio -
            (o.scale_macro_density ? (1 - current_target) * macro_area : 0);
        const double freed = o.vanished_memory
                                 ? fixed_macro_area - excluded_macro_area
                                 : (transparent ? fixed_macro_area : 0);
        const double budget =
            o.filler_cell_space
                ? (1 - o.target_density) * cell_space * o.filler_area_fraction
                : std::max(0., current_target *
                                       (width * height - fixed_area + freed) -
                                   (real_charge - charge.materialized_area)) *
                      o.filler_area_fraction;
        const double effective = detail::reservation_filler_area(
            budget * (tuned ? tuning.filler_charge : 1), filler_rounding_loss,
            charge);
        require(r.filler_area > 0 || effective == 0,
                "density reservations require available filler particles");
        filler_scale = r.filler_area > 0 ? effective / r.filler_area : 0;
        for (int i = original; i < n; ++i)
          ratio[i] = filler_scale * d.node_size_x[i] * d.node_size_y[i] /
                     (sx[i] * sy[i]);
        r.effective_filler_area = effective;
      };
      if (reservation_accounting)
        apply_reservation_fillers(target);

      r.effective_filler_area = r.filler_area * filler_scale;
      if (o.macro_transparency.always ||
          (cycling && global_iter >= o.macro_transparency.warmup)) {
        for (int i = d.num_movable; i < original; ++i) {
          const bool hidden =
              transparent && !d.node_macro.empty() && d.node_macro[i];
          ratio[i] = !hidden && sx[i] * sy[i] > 0
                         ? target * d.node_size_x[i] * d.node_size_y[i] /
                               (sx[i] * sy[i])
                         : 0;
        }
      }
      if (o.routability) {
        for (int i = 0; i < d.num_movable; ++i) {
          const double factor = std::sqrt(inflation[i]);
          sx[i] = std::max(d.node_size_x[i] * factor,
                           sharp(i) ? 0. : std::sqrt(2.) * dx);
          sy[i] = std::max(d.node_size_y[i] * factor,
                           sharp(i) ? 0. : std::sqrt(2.) * dy);
          ox[i] = (d.node_size_x[i] - sx[i]) / 2;
          oy[i] = (d.node_size_y[i] - sy[i]) / 2;
          ratio[i] = d.node_size_x[i] * d.node_size_y[i] * inflation[i] /
                     (sx[i] * sy[i]);
          if (o.scale_macro_density && !d.node_macro.empty() &&
              d.node_macro[i]) {
            ratio[i] *=
                o.routability_charge_compensation ? target : o.target_density;
          }
        }
      }
      if (o.vanished_memory) {
        for (int i = d.num_movable; i < original; ++i) {
          if (d.node_macro.empty() || !d.node_macro[i]) {
            continue;
          }
          sx[i] = d.node_size_x[i] * exclusion_scale;
          sy[i] = d.node_size_y[i] * exclusion_scale;
          ox[i] = (d.node_size_x[i] - sx[i]) / 2;
          oy[i] = (d.node_size_y[i] - sy[i]) / 2;
          ratio[i] = exclusion_scale > 0 ? target : 0;
        }
      }
      detail::density_map(d, o.bins_x, o.bins_y, o.num_threads, sx, sy, ox, oy,
                          ratio, map);
      if (reservation_accounting) {
        if (reservation_map_dirty) {
          detail::reservation_density_map(d, o.bins_x, o.bins_y,
                                          o.exact_cell_density, o.num_threads,
                                          reserved_map);
          reservation_map_dirty = false;
        }
        for (size_t k = 0; k < map.size(); ++k)
          map[k] += reserved_map[k];
      }
      for (auto &v : map) {
        v = v / (dx * dy) - target;
        if (o.density_upper_bound) {
          v = std::max(0., v);
        }
      }
      const auto poisson_start = std::chrono::steady_clock::now();
      solver.solve(map, phi, fx, fy);
      r.poisson_seconds += std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - poisson_start)
                               .count();
      auto dg = detail::density_force(d, o.bins_x, o.bins_y, o.num_threads, sx,
                                      sy, ox, oy, ratio, fx, fy);
      if ((!o.density_upper_bound && !global_iter) ||
          (o.density_upper_bound && !density_initialized)) {
        double wsum = 0, dsum = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          wsum += std::abs(wg.x[i]) + std::abs(wg.y[i]);
          dsum += std::abs(dg.x[i]) + std::abs(dg.y[i]);
        }
        if (!o.density_upper_bound || dsum > 1e-12) {
          density_weight *= std::max(wsum, 1.) / std::max(dsum, 1e-12);
          density_initialized = true;
        }
      }
      std::vector<double> cell_density;
      // Calibrate against ordinary, already preconditioned motion before
      // adding independent local-density/routing/alignment objectives.
      const auto macro_pull_forces =
          detail::macro_pull_forces(d, macro_pulls, wg, dg, density_weight);
      double ordinary_macro_force = 0;
      if (local_density &&
          (local_weight > 0 || o.macro_cell_density_weight > 0 ||
           o.macro_density_pull_weight > 0)) {
        const auto begin = std::chrono::steady_clock::now();
        cell_density = local_density->cell_density(d);
        r.local_density_grid_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          begin)
                .count();
        if (o.macro_density_pull_weight > 0 ||
            o.macro_cell_density_weight > 0) {
          for (int i = 0; i < d.num_movable; ++i) {
            if (!d.node_macro.empty() && d.node_macro[i]) {
              ordinary_macro_force +=
                  std::abs(wg.x[i] + density_weight * dg.x[i]) +
                  std::abs(wg.y[i] + density_weight * dg.y[i]);
            }
          }
        }
      }
      if (local_density && local_weight > 0) {
        const auto begin = std::chrono::steady_clock::now();
        const auto lg = local_density->gradient(d, false, cell_density,
                                                transparent, exclusion_scale);
        double base = 0, force = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (d.node_macro.empty() || !d.node_macro[i]) {
            base += std::abs(wg.x[i] + density_weight * dg.x[i]) +
                    std::abs(wg.y[i] + density_weight * dg.y[i]);
            force += std::abs(lg.x[i]) + std::abs(lg.y[i]);
          }
        }
        // Calibrate to a physical reference at target density, not the
        // current field. An initially sparse design must not turn tiny noise
        // into wirelength-sized force. Cap aggregate pressure on dense cases.
        double cell_area = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (d.node_macro.empty() || !d.node_macro[i]) {
            cell_area += d.node_size_x[i] * d.node_size_y[i];
          }
        }
        const double pitch = std::sqrt(width * height) / o.local_density_bins;
        const double reference =
            cell_area * local_target * pitch / (2 * std::acos(-1.));
        const double scale =
            local_weight * base / std::max({reference, force, 1e-12});
        double applied_force = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          // Calibrate before attenuation; otherwise a uniformly critical
          // region would have its protection normalized away.
          const double cell_scale =
              spread_mobility.empty() ? scale : scale * spread_mobility[i];
          wg.x[i] += cell_scale * lg.x[i];
          wg.y[i] += cell_scale * lg.y[i];
          if (!spread_mobility.empty()) {
            applied_force +=
                cell_scale * (std::abs(lg.x[i]) + std::abs(lg.y[i]));
          }
        }
        r.local_density_force =
            std::max(r.local_density_force,
                     spread_mobility.empty() ? scale * force : applied_force);
        r.local_density_seconds += std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - begin)
                                       .count();
      }
      if (macro_cell_density) {
        const auto begin = std::chrono::steady_clock::now();
        const auto mg = macro_cell_density->gradient(d, true, cell_density);
        double force = 0, macro_area = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_macro.empty() && d.node_macro[i]) {
            force += std::abs(mg.x[i]) + std::abs(mg.y[i]);
            macro_area += d.node_size_x[i] * d.node_size_y[i];
          }
        }
        const double pitch = std::sqrt(width * height) / o.local_density_bins;
        const double reference =
            macro_area * o.target_density * pitch / (2 * std::acos(-1.));
        // Calibrate once per netlist, before other optional forces. Freezing
        // the coefficient lets weak cell contributions accumulate linearly;
        // renormalizing by the current field would cap away that accumulation.
        // A physical floor avoids amplifying an initially tiny sparse field.
        if (macro_cell_scale < 0 && ordinary_macro_force > 0) {
          macro_cell_scale = o.macro_cell_density_weight *
                             ordinary_macro_force /
                             std::max({reference, force, 1e-12});
          r.macro_cell_density_coefficient = macro_cell_scale;
        }
        const double scale = std::max(0., macro_cell_scale);
        r.macro_cell_density_force =
            std::max(r.macro_cell_density_force, scale * force);
        for (int i = 0; i < d.num_movable; ++i) {
          wg.x[i] += scale * mg.x[i];
          wg.y[i] += scale * mg.y[i];
        }
        r.macro_cell_density_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          begin)
                .count();
      }
      if (r.centroid.started &&
          (o.centroid_distance_attraction || parent_centroids) &&
          r.centroid.attraction_calibration_iteration < 0) {
        // Calibrate once against ordinary WL+density, before extra penalties.
        // At one circle radius the spring has the requested relative strength;
        // later displacement must increase it, not normalize it away.
        const auto reference =
            detail::centroid_attraction_reference(d, o.target_density);
        double base = 0, norm = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (reference[i] == 0) {
            continue;
          }
          base += std::pow(wg.x[i] + density_weight * dg.x[i], 2) +
                  std::pow(wg.y[i] + density_weight * dg.y[i], 2);
          norm += reference[i] * reference[i];
        }
        if (base > 0 && norm > 0) {
          r.centroid.attraction_coefficient =
              o.centroid_weight * std::sqrt(base / norm);
          r.centroid.attraction_calibration_iteration = global_iter;
        }
      }
      for (const auto &force : macro_pull_forces) {
        wg.x[force.node] += force.x;
        wg.y[force.node] += force.y;
      }
      if (local_density && o.macro_density_pull_weight > 0) {
        const auto begin = std::chrono::steady_clock::now();
        auto pg = local_density->density_pull(d, cell_density,
                                              o.macro_density_pull_range);
        double norm = 0, reference = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_macro.empty() && d.node_macro[i]) {
            norm += std::abs(pg.x[i]) + std::abs(pg.y[i]);
            reference +=
                d.node_size_x[i] * d.node_size_y[i] * o.target_density /
                std::max(1e-12, std::min(d.node_size_x[i], d.node_size_y[i]));
          }
        }
        const double scale = (centroids && !r.centroid.started)
                                 ? 0
                                 : o.macro_density_pull_weight *
                                       ordinary_macro_force /
                                       std::max({reference, norm, 1e-12});
        const double anchor_ramp =
            d.centroid_targets.empty()
                ? std::pow(std::clamp(double(global_iter + 1 -
                                             r.centroid.start_iteration) /
                                          o.centroid_ramp,
                                      0., 1.),
                           2)
                : 1;
        std::vector<double> anchor;
        if (r.centroid.attraction_coefficient > 0) {
          anchor = detail::centroid_attraction_reference(d, o.target_density);
        }
        double applied = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          double gx = scale * pg.x[i], gy = scale * pg.y[i];
          // Compare with the spring at one circle radius, not its current
          // force (zero at the centroid). This keeps the slope weaker without
          // disabling small shifts around an exactly centered macro.
          const double magnitude = std::hypot(gx, gy);
          if (!anchor.empty() && anchor[i] > 0 && magnitude > 0) {
            const double anchor_force =
                r.centroid.attraction_coefficient * anchor[i] * anchor_ramp;
            const double cap = o.macro_density_pull_weight * anchor_force;
            const double fraction = std::min(1., cap / magnitude);
            gx *= fraction;
            gy *= fraction;
            if (anchor_force > 0) {
              r.macro_density_pull_anchor_ratio =
                  std::max(r.macro_density_pull_anchor_ratio,
                           std::hypot(gx, gy) / anchor_force);
            }
          }
          wg.x[i] += gx;
          wg.y[i] += gy;
          applied += std::abs(gx) + std::abs(gy);
        }
        r.macro_density_pull_force =
            std::max(r.macro_density_pull_force, applied);
        r.macro_density_pull_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          begin)
                .count();
      }
      macro_routing_force = 0;
      macro_routing_channels = 0;
      if (macro_routing && !routing_utilization.empty()) {
        const auto routing_start = std::chrono::steady_clock::now();
        const auto pressure = detail::macro_routing_pressure(
            d, aux_x, aux_y, routing_utilization, o.routability_macro_range);
        double base = 0;
        int macros = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_macro[i]) {
            continue;
          }
          base += std::abs(wg.x[i] + density_weight * dg.x[i]) +
                  std::abs(wg.y[i] + density_weight * dg.y[i]);
          ++macros;
        }
        // Keep pressure magnitude: normalization by pressure itself would make
        // an empty channel as strong as a congested one. The unit-force floor
        // also lets an otherwise unconnected macro move out of a routing lane.
        const double scale =
            o.routability_macro_weight * std::max(1., base / macros);
        for (int i = 0; i < d.num_movable; ++i) {
          wg.x[i] += scale * pressure.x[i];
          wg.y[i] += scale * pressure.y[i];
          macro_routing_force +=
              scale * (std::abs(pressure.x[i]) + std::abs(pressure.y[i]));
        }
        r.macro_routing_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          routing_start)
                .count();
        macro_routing_channels = pressure.channels;
        macro_routing_force_peak =
            std::max(macro_routing_force_peak, macro_routing_force);
        macro_routing_channels_peak =
            std::max(macro_routing_channels_peak, macro_routing_channels);
      }
      if (fixed_macro_weight > 0 && !transparent) {
        const auto og = detail::fixed_macro_overlap(
            d, o.bins_x, o.bins_y, fixed_macro_weight * progress * progress,
            exclusion_scale);
        for (int i = 0; i < d.num_movable; ++i) {
          wg.x[i] += og.x[i];
          wg.y[i] += og.y[i];
        }
      }
      if (o.macro_align_weight > 0) {
        const auto align_start = std::chrono::steady_clock::now();
        if (o.macro_align_global || global_iter % 20 == 0) {
          align_pairs = alignment_pairs(d);
        }
        align_ramp = detail::macro_align_ramp(o, global_iter);
        double mean_x = 0, mean_y = 0;
        int macros = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_align_class.empty() && d.node_align_class[i] >= 0 &&
              !d.node_macro.empty() && d.node_macro[i]) {
            mean_x += std::abs(wg.x[i] + density_weight * dg.x[i]);
            mean_y += std::abs(wg.y[i] + density_weight * dg.y[i]);
            ++macros;
          }
        }
        const double scale =
            macros ? o.macro_align_weight * align_ramp / macros : 0;
        const auto ag = detail::macro_alignment(
            d, align_pairs, scale * mean_x, scale * mean_y, o.macro_align_sigma,
            &align_residual, o.macro_align_global);
        double align_force = 0;
        for (int i = 0; i < original; ++i) {
          wg.x[i] += ag.x[i];
          wg.y[i] += ag.y[i];
          align_force += std::abs(ag.x[i]) + std::abs(ag.y[i]);
        }
        align_force_peak = std::max(align_force_peak, align_force);
        r.macro_alignment_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          align_start)
                .count();
      }
      if (o.cohesion_weight > 0 && progress > 0 && !d.cohesion_groups.empty()) {
        const auto cg = detail::cohesion(d, o.cohesion_weight * progress);
        for (int i = 0; i < original; ++i) {
          wg.x[i] += cg.x[i];
          wg.y[i] += cg.y[i];
        }
      }
      if (o.macro_overlap_weight > 0) {
        const auto mg = detail::macro_overlap(
            d, o.macro_overlap_weight *
                   (1 + o.macro_align_overlap_boost * align_ramp));
        for (int i = 0; i < original; ++i) {
          wg.x[i] += mg.x[i];
          wg.y[i] += mg.y[i];
        }
      }
      if (r.centroid.started) {
        const auto begin = std::chrono::steady_clock::now();
        const auto cg = detail::centroid_attraction(d, 1);
        detail::Gradient rg;
        if (o.centroid_repulsion_weight > 0) {
          rg = detail::centroid_repulsion(d, 1, o.target_density);
        }
        double base = 0, attraction = 0, repulsion = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          const bool repels =
              o.centroid_repulsion_weight > 0 && (rg.x[i] != 0 || rg.y[i] != 0);
          if (cg.x[i] != 0 || cg.y[i] != 0 || repels) {
            if (repels) {
              repulsion += rg.x[i] * rg.x[i] + rg.y[i] * rg.y[i];
            }
            base += std::pow(wg.x[i] + density_weight * dg.x[i], 2) +
                    std::pow(wg.y[i] + density_weight * dg.y[i], 2);
            attraction += cg.x[i] * cg.x[i] + cg.y[i] * cg.y[i];
          }
        }
        const double ramp =
            std::clamp(double(global_iter + 1 - r.centroid.start_iteration) /
                           o.centroid_ramp,
                       0., 1.);
        // Preserve the existing attraction normalization. A perfectly compact
        // group still needs repulsion, so use its norm only in that zero case.
        const double norm = attraction > 0 ? attraction : repulsion;
        const double scale =
            norm > 0 ? o.centroid_weight * ramp * ramp * std::sqrt(base / norm)
                     : 0;
        detail::Gradient alignment_gradient;
        if (align_centroids) {
          const auto alignment_start = std::chrono::steady_clock::now();
          centroid_align_ramp = detail::macro_align_ramp(
              centroid_alignment_schedule, global_iter);
          const auto [fx, fy] =
              detail::centroid_force_reference(d, wg, dg, density_weight);
          const double weight = o.centroid_align_weight * centroid_align_ramp;
          alignment_gradient = detail::centroid_alignment(
              d, weight * fx, weight * fy, o.centroid_align_sigma,
              o.target_density, &centroid_align_pairs,
              &centroid_align_residual);
          r.centroid.alignment_seconds +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                            alignment_start)
                  .count();
        }
        detail::Gradient parent_gradient;
        if (parent_centroids) {
          const auto parent_start = std::chrono::steady_clock::now();
          parent_gradient = detail::centroid_parent_attraction(
              d,
              r.centroid.attraction_coefficient * ramp * ramp *
                  o.centroid_parent_weight,
              o.target_density, o.centroid_parent_outside);
          r.centroid.parent_seconds +=
              std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                            parent_start)
                  .count();
        }
        detail::Gradient distance_gradient;
        if (o.centroid_distance_attraction) {
          // Frozen targets take effect immediately in the second GP. Only the
          // moving-centroid first phase ramps from scratch.
          const double spring_ramp =
              d.centroid_targets.empty() ? ramp * ramp : 1;
          distance_gradient = detail::centroid_circle_attraction(
              d, r.centroid.attraction_coefficient * spring_ramp,
              o.target_density, o.centroid_attraction_outside,
              o.centroid_anchor_macros_only);
        }
        for (int i = 0; i < original; ++i) {
          const bool anchor_active =
              !o.centroid_anchor_macros_only ||
              (i < d.num_movable && !d.node_macro.empty() && d.node_macro[i]);
          if (anchor_active) {
            wg.x[i] += o.centroid_distance_attraction ? distance_gradient.x[i]
                                                      : scale * cg.x[i];
            wg.y[i] += o.centroid_distance_attraction ? distance_gradient.y[i]
                                                      : scale * cg.y[i];
          }
          // Keep circle repulsion on its existing independent schedule/scale.
          if (o.centroid_repulsion_weight > 0) {
            wg.x[i] += scale * centroid_repulsion_weight * rg.x[i];
            wg.y[i] += scale * centroid_repulsion_weight * rg.y[i];
          }
        }
        if (align_centroids) {
          double force = 0;
          for (int i = 0; i < original; ++i) {
            wg.x[i] += alignment_gradient.x[i];
            wg.y[i] += alignment_gradient.y[i];
            if (i < d.num_movable) {
              force += std::abs(alignment_gradient.x[i]) +
                       std::abs(alignment_gradient.y[i]);
            }
          }
          centroid_align_force = std::max(centroid_align_force, force);
        }
        if (parent_centroids) {
          double force = 0;
          for (int i = 0; i < original; ++i) {
            wg.x[i] += parent_gradient.x[i];
            wg.y[i] += parent_gradient.y[i];
            if (i < d.num_movable) {
              force += std::abs(parent_gradient.x[i]) +
                       std::abs(parent_gradient.y[i]);
            }
          }
          centroid_parent_force = std::max(centroid_parent_force, force);
        }
        r.centroid.force_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          begin)
                .count();
        r.centroid.effective_repulsion_weight = centroid_repulsion_weight;
        ++r.centroid.iterations;
      }
      double max_grad = 0, ss = 0, yy = 0;
      std::vector<double> gx(n), gy(n);
      for (int i = 0; i < n; ++i) {
        if (can_move(i)) {
          const double pre = placement_preconditioner(d, i, density_weight);
          double cell_density_weight = density_weight;
          if (o.routability && !spread_mobility.empty() && inflation[i] > 1) {
            // Keep the full inflated source in the field so neighbors still
            // see the hotspot. Attenuate only this receiver's extra charge:
            // unit-area response + mobility * virtual inflation response,
            // both integrated over the same inflated footprint.
            cell_density_weight *=
                (1 + (inflation[i] - 1) * spread_mobility[i]) / inflation[i];
          }
          gx[i] = (wg.x[i] + cell_density_weight * dg.x[i]) / pre;
          gy[i] = (wg.y[i] + cell_density_weight * dg.y[i]) / pre;
          if (!finite(gx[i]) || !finite(gy[i])) {
            throw std::runtime_error("nonfinite placement gradient");
          }
          max_grad =
              std::max({max_grad, std::abs(gx[i]) / dx, std::abs(gy[i]) / dy});
          ss += std::pow(d.node_x[i] - old_ex[i], 2) +
                std::pow(d.node_y[i] - old_ey[i], 2);
          yy += std::pow(gx[i] - old_gx[i], 2) + std::pow(gy[i] - old_gy[i], 2);
        }
      }
      double cap = o.max_step / std::max(max_grad, 1e-12);
      alpha = global_iter && !restart && yy > 1e-20
                  ? std::clamp(std::sqrt(ss / yy), cap * 0.01, cap)
                  : cap;
      restart = false;
      old_ex = d.node_x;
      old_ey = d.node_y;
      old_gx = gx;
      old_gy = gy;
      prev_x = x;
      prev_y = y;
      double macro_momentum_work = 0;
      for (int i = 0; i < n; ++i) {
        if (can_move(i)) {
          x[i] = d.node_x[i] - alpha * gx[i];
          y[i] = d.node_y[i] - alpha * gy[i];
          project(i, x[i], y[i]);
          if (o.magnet_runtime && i < original && !d.node_macro.empty() &&
              d.node_macro[i]) {
            macro_momentum_work +=
                gx[i] * (x[i] - prev_x[i]) + gy[i] * (y[i] - prev_y[i]);
          }
        }
      }
      // Adaptive gradient restart: channel pressure changes with the macro
      // geometry and routing map. Retaining near-unit Nesterov momentum when
      // its proposed macro displacement goes uphill causes repeated channel
      // overshoot. Drop extrapolation on the next step without removing any
      // routing force or changing the step/curvature estimate.
      t = o.magnet_runtime && macro_momentum_work > 0 ? 1 : nt;
      r.iterations = iter + 1;
      // Decrement after this step, before installing newly requested holds.
      // A hold changes only motion, never net forces or density accounting.
      for (auto &remaining : movable_holds) {
        if (remaining > 0) {
          --remaining;
        }
      }
      if (movable_updates_active &&
          (global_iter + 1) % o.movable_interval == 0) {
        Design snapshot = model;
        std::copy_n(x.begin(), original, snapshot.node_x.begin());
        std::copy_n(y.begin(), original, snapshot.node_y.begin());
        const auto updates = o.movable_update(snapshot, global_iter + 1);
        std::vector<bool> seen(original);
        bool changed = false;
        for (const auto &position : updates) {
          const int i = position.node;
          require(i >= 0 && i < d.num_movable &&
                      (d.node_macro.empty() || !d.node_macro[i]) && !seen[i] &&
                      finite(position.x) && finite(position.y) &&
                      position.fix_iterations >= 0,
                  "invalid movable standard-cell update");
          double px = position.x, py = position.y;
          project(i, px, py);
          require(std::abs(px - position.x) <= 1e-6 &&
                      std::abs(py - position.y) <= 1e-6,
                  "movable update outside legal origin bounds");
          seen[i] = true;
          changed |= x[i] != px || y[i] != py;
          changed |= position.fix_iterations > 0;
        }
        for (const auto &position : updates) {
          const int i = position.node;
          x[i] = d.node_x[i] = model.node_x[i] = position.x;
          y[i] = d.node_y[i] = model.node_y[i] = position.y;
          movable_holds[i] = position.fix_iterations;
        }
        if (changed) {
          t = 1;
          prev_x = x;
          prev_y = y;
          restart = true;
          last_checkpoint = global_iter + 1;
          convergence.previous_hpwl = -1;
          convergence.previous_inflation = -1;
          convergence.stable_samples = 0;
          previous_hpwl = -1;
          best_feasible_x.clear();
          best_feasible_y.clear();
          best_feasible_hpwl = std::numeric_limits<double>::infinity();
        }
      }
      if (o.terminal_update && o.terminal_interval > 0 &&
          (global_iter + 1) % o.terminal_interval == 0) {
        Design snapshot = model;
        std::copy_n(x.begin(), original, snapshot.node_x.begin());
        std::copy_n(y.begin(), original, snapshot.node_y.begin());
        const auto updates = o.terminal_update(snapshot, global_iter + 1);
        std::vector<bool> seen(original);
        bool changed = false;
        for (const auto &pin : updates) {
          require(pin.node >= d.num_movable && pin.node < original &&
                      d.node_size_x[pin.node] == 0 &&
                      d.node_size_y[pin.node] == 0 && !seen[pin.node] &&
                      finite(pin.x) && finite(pin.y),
                  "invalid fixed zero-area terminal update");
          seen[pin.node] = true;
          changed |= x[pin.node] != pin.x || y[pin.node] != pin.y;
        }
        for (const auto &pin : updates) {
          x[pin.node] = d.node_x[pin.node] = model.node_x[pin.node] = pin.x;
          y[pin.node] = d.node_y[pin.node] = model.node_y[pin.node] = pin.y;
        }
        if (changed) {
          if (o.magnet_runtime) {
            o.magnet_runtime->event(global_iter + 1, MagnetEvent::terminal);
          }
          t = 1;
          prev_x = x;
          prev_y = y;
          restart = true;
          last_checkpoint = global_iter + 1;
          convergence.previous_hpwl = -1;
          convergence.previous_inflation = -1;
          convergence.stable_samples = 0;
          previous_hpwl = -1;
          best_feasible_x.clear();
          best_feasible_y.clear();
          best_feasible_hpwl = std::numeric_limits<double>::infinity();
        }
      }
      if (o.macro_pull_update && o.macro_pull_interval > 0 &&
          (global_iter + 1) % o.macro_pull_interval == 0) {
        Design snapshot = model;
        std::copy_n(x.begin(), original, snapshot.node_x.begin());
        std::copy_n(y.begin(), original, snapshot.node_y.begin());
        std::optional<std::vector<MacroPull>> replacement;
        if (!o.macro_pull_after_density ||
            overflow(snapshot, o.bins_x, o.bins_y, o.target_density,
                     o.scale_macro_density,
                     o.macro_transparency.always) <= o.stop_overflow) {
          replacement = o.macro_pull_update(snapshot, global_iter + 1);
        }
        if (replacement) {
          std::vector<bool> seen(original);
          for (const auto &pull : *replacement) {
            require(pull.node >= 0 && pull.node < d.num_movable &&
                        !d.node_macro.empty() && d.node_macro[pull.node] &&
                        !seen[pull.node] && finite(pull.x) && finite(pull.y) &&
                        finite(pull.weight) && pull.weight >= 0 &&
                        pull.weight <= 1,
                    "invalid virtual macro pull");
            seen[pull.node] = true;
          }
          macro_pulls = std::move(*replacement);
          if (o.magnet_runtime) {
            o.magnet_runtime->event(global_iter + 1, MagnetEvent::pull);
          }
          macro_pull_changed = global_iter + 1;
          t = 1;
          prev_x = x;
          prev_y = y;
          restart = true;
          convergence.previous_hpwl = -1;
          convergence.previous_inflation = -1;
          convergence.stable_samples = 0;
          previous_hpwl = -1;
          best_feasible_x.clear();
          best_feasible_y.clear();
          best_feasible_hpwl = std::numeric_limits<double>::infinity();
        }
      }
      const bool forced_checkpoint =
          density_gated_checkpoints &&
          (iter + 1 == o.max_iterations ||
           (convergence.deadline && global_iter + 1 == convergence.deadline));
      if (o.checkpoint && o.checkpoint_interval > 0 &&
          global_iter + 1 >= o.checkpoint_start &&
          ((global_iter + 1) % o.checkpoint_interval == 0 ||
           forced_checkpoint ||
           (o.checkpoint_final && iter + 1 == o.max_iterations))) {
        Design snapshot = model;
        std::copy_n(x.begin(), original, snapshot.node_x.begin());
        std::copy_n(y.begin(), original, snapshot.node_y.begin());
        bool invoke = true;
        if (density_gated_checkpoints) {
          const double measured =
              overflow(snapshot, o.bins_x, o.bins_y, o.target_density,
                       o.scale_macro_density, o.macro_transparency.always,
                       exclusion_scale);
          invoke = forced_checkpoint || measured <= o.stop_overflow;
          r.checkpoint_history.push_back(
              {global_iter + 1, measured, invoke, forced_checkpoint});
          if (o.checkpoint_observer) {
            o.checkpoint_observer(r.checkpoint_history.back());
          }
        }
        GpUpdate update;
        if (invoke) {
          checkpoint_invoked = true;
          update = o.checkpoint(snapshot, global_iter + 1);
          if (o.magnet_runtime) {
            o.magnet_runtime->event(global_iter + 1, MagnetEvent::timing);
          }
        }
        require(update.action == GpCheckpointAction::update ||
                    o.checkpoint_rollback,
                "GP checkpoint transaction requires checkpoint_rollback");
        if (o.checkpoint_rollback) {
          require(
              !update.replacement && update.previous_node.empty() &&
                  update.node_size_x.empty() && update.node_size_y.empty() &&
                  update.pin_offset_x.empty() && update.pin_offset_y.empty(),
              "GP rollback requires unchanged topology and footprint");
        }
        if (update.action == GpCheckpointAction::restore) {
          require(retained_checkpoint.has_value(),
                  "GP rollback has no retained checkpoint");
          require(!update.stop_movable_updates && !update.timing_violation &&
                      update.net_weights.empty() && !update.pin_timing &&
                      !update.density_reservations,
                  "GP restore cannot carry objective updates");
          const auto &saved = *retained_checkpoint;
          const auto &s = saved.state;
          d = saved.design;
          model = saved.model;
          reservation_accounting = s.reservation_accounting;
          reservation_map_dirty = true;
          sx = saved.sx;
          sy = saved.sy;
          ox = saved.ox;
          oy = saved.oy;
          ratio = saved.ratio;
          r.effective_filler_area = saved.effective_filler_area;
          x = s.x;
          y = s.y;
          prev_x = s.prev_x;
          prev_y = s.prev_y;
          old_gx = s.old_gx;
          old_gy = s.old_gy;
          old_ex = s.old_ex;
          old_ey = s.old_ey;
          density_weight = s.density_weight;
          density_initialized = s.density_initialized;
          previous_overflow = s.previous_overflow;
          previous_hpwl = s.previous_hpwl;
          t = s.momentum;
          alpha = s.step;
          convergence = s.convergence;
          inflation = s.inflation;
          routing_utilization = s.routing_utilization;
          active_route = s.route_model;
          route_hold = s.route_hold;
          probe_deadline = s.route_probe_deadline;
          r.route_probed = s.route_probed;
          r.route_switch_iteration = s.route_switch_iteration;
          r.route_switch_rudy = s.route_switch_rudy;
          r.route_switch_rsmt = s.route_switch_rsmt;
          align_pairs = s.align_pairs;
          restart = saved.restart;
          pin_timing_seen = saved.pin_timing_seen;
          spread_mobility = make_spread_mobility();
          r.overflow = saved.overflow;
          r.objective_converged = r.converged = false;
          restored_selection = s.iteration;
          r.rollback_history.push_back({global_iter + 1, s.iteration});
          // Policy clocks and the iteration cap are monotonic. Suppress
          // convergence on this rejected observation without changing the
          // restored optimizer, fillers, weights or density normalization.
          last_checkpoint = global_iter + 1;
          continue;
        }
        retain_checkpoint =
            invoke && update.action == GpCheckpointAction::retain;
        if (update.stop_movable_updates) {
          movable_updates_active = false;
        }
        if (update.timing_violation) {
          require(finite(*update.timing_violation) &&
                      *update.timing_violation >= 0,
                  "invalid checkpoint timing violation");
          convergence.timing_observed = true;
          convergence.timing_satisfied = *update.timing_violation == 0;
        }
        bool objective_changed = update.replacement.has_value();
        if (update.replacement) {
          require(!o.gp_search.enabled,
                  "GP A/B search requires fixed topology checkpoints");
          require(!movable_updates_active,
                  "retire movable updates before topology checkpoints");
          macro_pulls.clear();
          macro_pull_changed = -1;
          require(update.net_weights.empty() && update.node_size_x.empty() &&
                      update.node_size_y.empty() &&
                      update.pin_offset_x.empty() &&
                      update.pin_offset_y.empty() && !update.pin_timing &&
                      !update.density_reservations,
                  "GP replacement cannot be combined with field updates");
          Design candidate = std::move(*update.replacement);
          candidate.pin_timing.reset();
          candidate.validate();
          const auto &a = d.layout;
          const auto &b = candidate.layout;
          require(a.xl == b.xl && a.yl == b.yl && a.xh == b.xh &&
                      a.yh == b.yh && a.site_width == b.site_width &&
                      a.row_height == b.row_height,
                  "GP replacement changed layout");
          const int size = count(candidate.node_x.size());
          require(update.previous_node.size() == candidate.node_x.size(),
                  "GP replacement missing node mapping");
          require(candidate.centroid_groups.size() == d.centroid_groups.size(),
                  "GP replacement changed centroid candidates");
          for (size_t g = 0; g < d.centroid_groups.size(); ++g) {
            require(candidate.centroid_groups[g].name ==
                            d.centroid_groups[g].name &&
                        candidate.centroid_groups[g].parent ==
                            d.centroid_groups[g].parent,
                    "GP replacement changed centroid identity");
          }
          std::vector<bool> seen(original);
          std::vector<double> next_inflation(size, 1);
          const auto is_macro = [](const Design &v, int i) {
            return !v.node_macro.empty() && v.node_macro[i];
          };
          const auto pin_geometry = [](const Design &v, int i) {
            std::vector<std::pair<double, double>> pins;
            for (int k = v.flat_node2pin_start[i];
                 k < v.flat_node2pin_start[i + 1]; ++k) {
              const int p = v.flat_node2pin[k];
              pins.emplace_back(v.pin_offset_x[p], v.pin_offset_y[p]);
            }
            std::sort(pins.begin(), pins.end());
            return pins;
          };
          for (int i = 0; i < size; ++i) {
            const int old = update.previous_node[i];
            require(old >= -1 && old < original,
                    "invalid GP replacement node mapping");
            require(candidate.node_size_x[i] <= width &&
                        candidate.node_size_y[i] <= height,
                    "GP replacement cell does not fit layout");
            if (old < 0) {
              require(i < candidate.num_movable && !is_macro(candidate, i) &&
                          !candidate.node_macro.empty(),
                      "GP replacement inserted a macro or fixed node");
              continue;
            }
            require(candidate.placement_mask(i) == d.placement_mask(old),
                    "GP replacement changed placement-mask semantics");
            require(!seen[old], "duplicate GP replacement node mapping");
            seen[old] = true;
            require((i < candidate.num_movable) == (old < d.num_movable) &&
                        is_macro(candidate, i) == is_macro(d, old),
                    "GP replacement changed node mobility/type");
            if (old >= d.num_movable || is_macro(d, old)) {
              require(candidate.node_size_x[i] == d.node_size_x[old] &&
                          candidate.node_size_y[i] == d.node_size_y[old] &&
                          pin_geometry(candidate, i) == pin_geometry(d, old),
                      "GP replacement changed macro/fixed geometry");
              if (old >= d.num_movable) {
                require(candidate.node_x[i] == x[old] &&
                            candidate.node_y[i] == y[old],
                        "GP replacement moved a fixed node");
              }
            }
            candidate.node_x[i] = x[old];
            candidate.node_y[i] = y[old];
            next_inflation[i] = inflation[old];
          }
          for (int i = 0; i < original; ++i) {
            require(seen[i] || (i < d.num_movable && !is_macro(d, i) &&
                                !d.node_macro.empty()),
                    "GP replacement removed a macro or fixed node");
          }
          for (int i = original; i < n; ++i) {
            candidate.node_x.push_back(x[i]);
            candidate.node_y.push_back(y[i]);
            candidate.node_size_x.push_back(d.node_size_x[i]);
            candidate.node_size_y.push_back(d.node_size_y[i]);
            next_inflation.push_back(inflation[i]);
          }
          candidate.build_maps();
          d = std::move(candidate);
          original = size;
          movable_holds.assign(d.num_movable, 0);
          n = count(d.node_x.size());
          inflation = std::move(next_inflation);
          x = d.node_x;
          y = d.node_y;
          sx.resize(n);
          sy.resize(n);
          ox.resize(n);
          oy.resize(n);
          ratio.resize(n);
          old_gx.assign(n, 0);
          old_gy.assign(n, 0);
          old_ex = x;
          old_ey = y;
          align_pairs = alignment_pairs(d);
          // A saved feasible point belongs to the old netlist and dimensions.
          best_feasible_x.clear();
          best_feasible_y.clear();
          best_feasible_hpwl = std::numeric_limits<double>::infinity();
          // Recompute density-step scaling for the new graph. Memory geometry
          // remains fixed; the next iteration reapplies its current schedule.
          density_initialized = false;
          alpha = 0;
        } else {
          Design candidate = model;
          auto assign = [](auto &to, const auto &from) {
            if (!from.empty()) {
              require(from.size() == to.size(),
                      "GP checkpoint array size changed");
              to = from;
            }
          };
          assign(candidate.net_weights, update.net_weights);
          if (update.pin_timing)
            candidate.pin_timing = update.pin_timing;
          if (update.density_reservations)
            candidate.density_reservations = *update.density_reservations;
          assign(candidate.node_size_x, update.node_size_x);
          assign(candidate.node_size_y, update.node_size_y);
          assign(candidate.pin_offset_x, update.pin_offset_x);
          assign(candidate.pin_offset_y, update.pin_offset_y);
          candidate.validate();
          for (int i = 0; i < original; ++i) {
            const bool resizable =
                i < d.num_movable && !d.node_macro.empty() && !d.node_macro[i];
            if (!resizable) {
              require(candidate.node_size_x[i] == d.node_size_x[i] &&
                          candidate.node_size_y[i] == d.node_size_y[i],
                      "GP checkpoint resized a macro or fixed node");
              for (int k = d.flat_node2pin_start[i];
                   k < d.flat_node2pin_start[i + 1]; ++k) {
                const int pin = d.flat_node2pin[k];
                require(candidate.pin_offset_x[pin] == d.pin_offset_x[pin] &&
                            candidate.pin_offset_y[pin] == d.pin_offset_y[pin],
                        "GP checkpoint changed a macro or fixed pin");
              }
            }
            require(candidate.node_size_x[i] <= width &&
                        candidate.node_size_y[i] <= height,
                    "GP checkpoint cell does not fit layout");
          }
          const bool pin_timing_changed = candidate.pin_timing != d.pin_timing;
          const bool reservations_changed =
              candidate.density_reservations != d.density_reservations;
          const bool changed = reservations_changed ||
                               candidate.net_weights != d.net_weights ||
                               pin_timing_changed ||
                               candidate.node_size_x != model.node_size_x ||
                               candidate.node_size_y != model.node_size_y ||
                               candidate.pin_offset_x != d.pin_offset_x ||
                               candidate.pin_offset_y != d.pin_offset_y;
          if (changed) {
            d.net_weights = candidate.net_weights;
            d.pin_timing = candidate.pin_timing;
            d.density_reservations = candidate.density_reservations;
            pin_timing_seen =
                pin_timing_seen || (d.pin_timing && d.pin_timing->enabled());
            if (pin_timing_changed || reservations_changed) {
              best_feasible_x.clear();
              best_feasible_y.clear();
              best_feasible_hpwl = std::numeric_limits<double>::infinity();
            }
            std::copy_n(candidate.node_size_x.begin(), original,
                        d.node_size_x.begin());
            std::copy_n(candidate.node_size_y.begin(), original,
                        d.node_size_y.begin());
            d.pin_offset_x = candidate.pin_offset_x;
            d.pin_offset_y = candidate.pin_offset_y;
          }
          objective_changed = changed;
        }
        if (objective_changed) {
          reservation_accounting =
              reservation_accounting || !d.density_reservations.empty();
          reservation_map_dirty = true;
          spread_mobility = make_spread_mobility();
          area = 0;
          double inflated = 0;
          for (int i = 0; i < original; ++i) {
            const double a = d.node_size_x[i] * d.node_size_y[i];
            if (i < d.num_movable) {
              area += a;
              inflated += a * inflation[i];
            }
            if (!d.centroid_node_area.empty() &&
                (d.node_macro.empty() || !d.node_macro[i])) {
              d.centroid_node_area[i] = physical_area(i);
            }
            if (active(i)) {
              project(i, x[i], y[i]);
            }
          }
          inflated_ratio = area ? inflated / area : 1;
          const double capacity =
              o.target_density / (o.improved_routability ? inflated_ratio : 1);
          for (int i = 0; i < n; ++i) {
            const double factor = o.routability ? std::sqrt(inflation[i]) : 1;
            sx[i] = std::max(d.node_size_x[i] * factor,
                             active(i) && !sharp(i) ? std::sqrt(2.) * dx : 0.);
            sy[i] = std::max(d.node_size_y[i] * factor,
                             active(i) && !sharp(i) ? std::sqrt(2.) * dy : 0.);
            ox[i] = (d.node_size_x[i] - sx[i]) / 2;
            oy[i] = (d.node_size_y[i] - sy[i]) / 2;
            ratio[i] = sx[i] * sy[i] != 0
                           ? d.node_size_x[i] * d.node_size_y[i] *
                                 inflation[i] / (sx[i] * sy[i])
                           : 0;
            // Fixed obstacles keep target_density, as in the no-checkpoint
            // path; only scaled movable macros follow the compensation.
            if (!active(i)) {
              ratio[i] *= o.target_density;
            } else if (o.scale_macro_density && i < original &&
                       !d.node_macro.empty() && d.node_macro[i]) {
              ratio[i] *= o.routability_charge_compensation ? capacity
                                                            : o.target_density;
            }
          }
          for (auto &g : d.cohesion_groups) {
            if (!g.areas.empty()) {
              for (size_t k = 0; k < g.nodes.size(); ++k) {
                const int i = g.nodes[k];
                if (d.node_macro.empty() || !d.node_macro[i]) {
                  g.areas[k] = physical_area(i);
                }
              }
            }
          }
          if (reservation_accounting)
            apply_reservation_fillers(capacity);
          model = d;
          model.node_x.assign(x.begin(), x.begin() + original);
          model.node_y.assign(y.begin(), y.begin() + original);
          model.node_size_x.resize(original);
          model.node_size_y.resize(original);
          model.build_maps();
          if (local_density) {
            local_density = std::make_unique<detail::LocalDensity>(
                model, o.local_density_bins, o.local_density_range,
                local_target, o.num_threads);
          }
          rebuild_macro_cell_density();
          t = 1;
          prev_x = x;
          prev_y = y;
          restart = true;
          last_checkpoint = global_iter + 1;
          // New weights/netlist make successive objective values incomparable.
          // Do not restart the hard continuation deadline after an ECO.
          convergence.previous_hpwl = -1;
          convergence.previous_inflation = -1;
          convergence.stable_samples = 0;
          previous_overflow = -1;
          previous_hpwl = hpwl(model);
        }
      }
      if (step_recorder) {
        step_recorder->observe(model, x, y, global_iter + 1);
      }
      if (o.magnet_runtime) {
        o.magnet_runtime->observe(model, x, y, global_iter + 1);
      }
      if (o.macro_motion_observer) {
        o.macro_motion_observer(model, x, y, global_iter + 1);
      }
      if (o.observer && r.iterations % o.observer_interval == 0) {
        Design physical = model;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        o.observer(physical, PlacementStage::global_placement,
                   offset + r.iterations);
      }
      if ((global_iter + 1) % 20 == 0) {
        Design physical = model;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        r.overflow = overflow(physical, o.bins_x, o.bins_y, o.target_density,
                              o.scale_macro_density,
                              o.macro_transparency.always, exclusion_scale);
        const double current_hpwl = hpwl(physical);
        const double selection_cost =
            current_hpwl +
            (fixed_macro_weight > 0
                 ? detail::fixed_macro_overlap(physical, o.bins_x, o.bins_y,
                                               fixed_macro_weight)
                       .value
                 : 0);
        if (keep_feasible && !pin_timing_seen && opacity_settled &&
            global_iter >= 39 && inflated_ratio == 1 &&
            (!fixed_macro_weight || progress >= 1) &&
            r.overflow <= o.stop_overflow &&
            selection_cost < best_feasible_hpwl) {
          best_feasible_hpwl = selection_cost;
          best_feasible_x = x;
          best_feasible_y = y;
          best_feasible_iteration = global_iter + 1;
        }
        const auto overlaps = detail::macro_overlap_counts(physical);
        if (o.macro_align_weight > 0) {
          // Report the current legal-geometry candidates, not stale pairs from
          // before this iteration's move. The next force refreshes
          // independently.
          if (o.macro_align_global) {
            align_pairs = alignment_pairs(physical);
          }
          detail::macro_alignment(physical, align_pairs, 0, 0,
                                  o.macro_align_sigma, &align_residual,
                                  o.macro_align_global);
        }
        r.history.push_back(
            {global_iter + 1, r.overflow, current_hpwl, density_weight,
             route_overflow, inflated_ratio, overlaps.first, overlaps.second,
             overflow(physical, o.bins_x, o.bins_y, o.target_density), target,
             filler_scale, align_pairs.size(), align_residual, align_ramp});
        r.history.back().macro_opacity = transparent ? 0 : 1;
        r.history.back().macro_exclusion_scale =
            transparent ? 0 : exclusion_scale;
        r.history.back().local_density_weight = local_weight;
        r.history.back().pin_timing_cost = pin_timing_cost;
        const auto sampled_charge = reservation_charge(d, inflation);
        r.history.back().reserved_area = sampled_charge.virtual_area;
        r.history.back().materialized_reserved_area =
            sampled_charge.materialized_area;
        r.history.back().effective_filler_area = r.effective_filler_area;
        r.history.back().filler_count = r.filler_count;
        r.history.back().filler_area = r.filler_area;
        r.history.back().filler_size_x = r.filler_size_x;
        r.history.back().filler_size_y = r.filler_size_y;
        if (align_centroids) {
          detail::centroid_alignment(physical, 0, 0, o.centroid_align_sigma,
                                     o.target_density, &centroid_align_pairs,
                                     &centroid_align_residual);
          auto &sample = r.history.back();
          sample.centroid_align_pairs = centroid_align_pairs;
          sample.centroid_align_residual = centroid_align_residual;
          sample.centroid_align_ramp = centroid_align_ramp;
          sample.centroid_align_force = centroid_align_force;
          centroid_align_force = 0;
        }
        if (parent_centroids) {
          detail::centroid_parent_attraction(
              physical, 0, o.target_density, o.centroid_parent_outside,
              &centroid_parent_edges, &centroid_parent_distance);
          auto &sample = r.history.back();
          sample.centroid_parent_edges = centroid_parent_edges;
          sample.centroid_parent_max_distance = centroid_parent_distance;
          sample.centroid_parent_force = centroid_parent_force;
          centroid_parent_force = 0;
        }
        r.history.back().align_force = align_force_peak;
        align_force_peak = 0;
        r.history.back().macro_routing_force = macro_routing_force_peak;
        r.history.back().macro_routing_channels = macro_routing_channels_peak;
        macro_routing_force_peak = 0;
        macro_routing_channels_peak = 0;
        if (separate_circles) {
          const auto circles =
              detail::centroid_circle_overlaps(physical, o.target_density);
          auto &sample = r.history.back();
          sample.centroid_overlaps = circles.overlaps;
          sample.centroid_max_overlap = circles.maximum;
          sample.centroid_repulsion_weight = centroid_repulsion_weight;
          r.centroid.circle_overlaps = circles.overlaps;
          r.centroid.max_circle_overlap = circles.maximum;
          r.centroid.circles_converged =
              circles.maximum <= o.centroid_stop_overlap;
          r.centroid.effective_repulsion_weight = centroid_repulsion_weight;
        }
        const bool settled =
            opacity_settled && global_iter + 1 >= o.min_iterations &&
            r.overflow <= o.stop_overflow &&
            (!o.macro_align_weight || global_iter >= o.macro_align_full) &&
            (!align_centroids || global_iter >= o.centroid_align_full) &&
            (!fixed_macro_weight || progress >= 1);
        const bool geometry_feasible =
            settled && (!centroids || global_iter + 1 >= o.centroid_ramp) &&
            (!separate_circles || r.centroid.circles_converged);
        const bool macro_pulls_ready =
            macro_pull_changed < 0 || global_iter + 1 - macro_pull_changed >=
                                          o.macro_pull_settle_iterations;
        const bool geometry_ready =
            geometry_feasible && global_iter + 1 > last_checkpoint &&
            std::none_of(movable_holds.begin(), movable_holds.end(),
                         [](int hold) { return hold > 0; }) &&
            macro_pulls_ready &&
            (!density_gated_checkpoints || checkpoint_invoked);
        if (o.post_density_factor) {
          const int iteration = global_iter + 1;
          if (geometry_feasible && !convergence.density_iteration) {
            convergence.density_iteration = iteration;
            convergence.deadline = static_cast<int>(std::min<std::int64_t>(
                std::numeric_limits<int>::max(),
                std::int64_t(iteration) * (1 + o.post_density_factor)));
            convergence.previous_hpwl = -1;
            convergence.previous_inflation = -1;
            convergence.stable_samples = 0;
          }
          const auto close = [&](double previous, double current) {
            return previous >= 0 &&
                   std::abs(current - previous) <=
                       o.convergence_tolerance *
                           std::max(
                               {std::abs(previous), std::abs(current), 1e-12});
          };
          const bool stable =
              settled && close(convergence.previous_hpwl, current_hpwl) &&
              close(convergence.previous_inflation, inflated_ratio);
          convergence.stable_samples =
              stable ? convergence.stable_samples + 1 : 0;
          convergence.previous_hpwl = current_hpwl;
          convergence.previous_inflation = inflated_ratio;
        }
        const bool route_candidate =
            o.post_density_factor
                ? convergence.stable_samples >= 2
                : (geometry_ready &&
                   (o.route.probe || active_route == RouteModel::uroute ||
                    learned_ladder));
        if (route_candidate && o.routability) {
          auto route = estimate_route(physical, global_iter + 1);
          const auto probe_model =
              o.route.probe_model == RouteModel::uroute && !o.uroute
                  ? RouteModel::rsmt
                  : o.route.probe_model;
          if (o.route.probe && !r.route_probed && geometry_ready &&
              convergence.timing_satisfied &&
              active_route != RouteModel::uroute &&
              (probe_model == RouteModel::uroute ||
               active_route == RouteModel::rudy)) {
            r.route_probed = true;
            if (o.route_probe_observer) {
              if (o.route_physical_snapshot) {
                const auto expanded = o.route_physical_snapshot(physical);
                o.route_probe_observer(expanded, aux_x, aux_y);
              } else
                o.route_probe_observer(physical, aux_x, aux_y);
            }
            auto probe = measure_route(physical, probe_model);
            r.route_probe_overflow = probe.overflow;
            if (probe_model == RouteModel::uroute) {
              r.route_switch_uroute_rsmt =
                  active_route == RouteModel::rsmt
                      ? route.overflow
                      : measure_route(physical, RouteModel::rsmt).overflow;
              r.route_switch_uroute = probe.overflow;
            } else {
              r.route_switch_rudy = route.overflow;
              r.route_switch_rsmt = probe.overflow;
            }
            if (probe.overflow > o.route.target_overflow) {
              active_route = r.route_model = probe_model;
              if (probe_model == RouteModel::uroute)
                r.route_uroute_switch_iteration = global_iter + 1;
              else
                r.route_switch_iteration = global_iter + 1;
              probe_deadline =
                  std::min(o.max_iterations + offset,
                           global_iter + 1 + o.route.probe_iterations);
              convergence.stable_samples = 0;
              route = std::move(probe);
            }
          }
          route_peak_utilization = route.max_utilization;
          route_overflow = route.overflow;
          r.history.back().route_overflow = route_overflow;
        }
        // A qualifying RSMT candidate must collect its hold window rather
        // than exit before the learned model has a chance to measure it.
        const bool uroute_pending =
            learned_ladder && active_route == RouteModel::rsmt &&
            route_overflow <= o.route.uroute_switch_overflow;
        if (o.post_density_factor || probe_deadline ||
            active_route == RouteModel::uroute || learned_ladder) {
          r.objective_converged =
              (!density_gated_checkpoints || checkpoint_invoked) &&
              !uroute_pending && route_candidate &&
              convergence.stable_samples >= (o.post_density_factor ? 2 : 0) &&
              (!o.routability || route_overflow <= o.route.target_overflow ||
               (r.route_probed && r.route_probe_overflow >= 0 &&
                r.route_probe_overflow <= o.route.target_overflow)) &&
              convergence.timing_satisfied;
          r.continuation_limit =
              probe_deadline ? global_iter + 1 >= probe_deadline
                             : (convergence.deadline &&
                                global_iter + 1 >= convergence.deadline);
          r.density_iteration = convergence.density_iteration;
          r.continuation_deadline =
              probe_deadline ? probe_deadline : convergence.deadline;
          r.timing_observed = convergence.timing_observed;
          r.timing_satisfied = convergence.timing_satisfied;
        }
        route_history(r.history.back());
        if (o.magnet_runtime) {
          o.magnet_runtime->sample(physical, r.history.back());
        }
        if (o.progress_observer) {
          o.progress_observer(r.history.back());
        }
        if (o.magnet_runtime && !o.magnet_runtime->pending.empty()) {
          break;
        }
        const bool magnet_ready =
            !o.magnet_runtime || o.magnet_runtime->ready(global_iter + 1);
        if (!o.trial_tuning && magnet_ready && geometry_ready &&
            ((!o.post_density_factor && !probe_deadline &&
              active_route != RouteModel::uroute && !learned_ladder) ||
             r.objective_converged)) {
          r.converged = true;
          break;
        }
        if (!o.trial_tuning && r.continuation_limit && macro_pulls_ready &&
            magnet_ready) {
          break;
        }
        if (separate_circles && global_iter + 1 >= o.centroid_ramp &&
            !r.centroid.circles_converged) {
          const double improvement =
              (previous_circle_overlap - r.centroid.max_circle_overlap) /
              std::max(previous_circle_overlap, 1e-12);
          if (improvement < .1) {
            centroid_repulsion_weight =
                std::min(centroid_repulsion_weight * 2,
                         o.centroid_repulsion_weight * 10000);
          }
          previous_circle_overlap = r.centroid.max_circle_overlap;
        }
        double growth = 1.15;
        if (o.adaptive_density) {
          const double improvement =
              previous_overflow > 0
                  ? (previous_overflow - r.overflow) / previous_overflow
                  : 1;
          growth = improvement < .05 ? 1.4 : 1.05;
          if (previous_hpwl > 0 && current_hpwl > previous_hpwl * 1.2) {
            growth = std::min(growth, 1.05);
          }
        }
        previous_overflow = r.overflow;
        previous_hpwl = current_hpwl;
        // If only circle separation is pending, do not keep escalating an
        // already-satisfied density constraint. Resume if density rises again.
        if ((!separate_circles && !o.density_upper_bound &&
             !o.post_density_factor) ||
            r.overflow > o.stop_overflow) {
          density_weight = std::min(density_weight * growth, 1e12);
        }
      }
      if (retain_checkpoint) {
        retained_checkpoint = RetainedCheckpoint{d,
                                                 model,
                                                 save_state(global_iter + 1),
                                                 restart,
                                                 pin_timing_seen,
                                                 r.overflow,
                                                 sx,
                                                 sy,
                                                 ox,
                                                 oy,
                                                 ratio,
                                                 r.effective_filler_area};
      }
    }
    r.centroid.completed = r.centroid.started;
    if (o.magnet_runtime) {
      o.magnet_runtime->density_initialized = density_initialized;
      o.magnet_runtime->routing =
          detail::MagnetRuntime::RouteProgress{active_route,
                                               probe_deadline,
                                               r.route_switch_iteration,
                                               r.route_switch_rudy,
                                               r.route_switch_rsmt,
                                               r.route_probed,
                                               r.route_uroute_switch_iteration,
                                               r.route_switch_uroute_rsmt,
                                               r.route_switch_uroute,
                                               r.route_probe_overflow};
      o.magnet_runtime->convergence = convergence;
      o.magnet_runtime->checkpoint_invoked = checkpoint_invoked;
      o.magnet_runtime->macro_pulls = std::move(macro_pulls);
      o.magnet_runtime->macro_pull_changed = macro_pull_changed;
    }
    if (density_gated_checkpoints || o.magnet_runtime) {
      // A forced cap need not coincide with a 20-iteration history sample.
      r.timing_observed = convergence.timing_observed;
      r.timing_satisfied = convergence.timing_satisfied;
      r.density_iteration = convergence.density_iteration;
      r.continuation_deadline =
          probe_deadline ? probe_deadline : convergence.deadline;
      if (!convergence.timing_satisfied) {
        r.objective_converged = false;
      }
    }
    if (state) {
      next.pin_timing = d.pin_timing;
      next.density_reservations = d.density_reservations;
      next.reservation_accounting = reservation_accounting;
      next.convergence = convergence;
      next.signature = signature;
      next.iteration = offset + r.iterations;
      next.x = x;
      next.y = y;
      next.prev_x = prev_x;
      next.prev_y = prev_y;
      next.old_gx = old_gx;
      next.old_gy = old_gy;
      next.old_ex = old_ex;
      next.old_ey = old_ey;
      next.density_weight = density_weight;
      next.density_initialized = density_initialized;
      next.previous_overflow = previous_overflow;
      next.previous_hpwl = previous_hpwl;
      next.momentum = t;
      next.step = alpha;
      next.inflation = inflation;
      next.routing_utilization = routing_utilization;
      next.route_model = active_route;
      next.route_hold = route_hold;
      next.route_uroute_hold = uroute_hold;
      next.route_probe_deadline = probe_deadline;
      next.route_probed = r.route_probed;
      next.route_switch_iteration = r.route_switch_iteration;
      next.route_switch_rudy = r.route_switch_rudy;
      next.route_switch_rsmt = r.route_switch_rsmt;
      next.route_uroute_switch_iteration = r.route_uroute_switch_iteration;
      next.route_switch_uroute_rsmt = r.route_switch_uroute_rsmt;
      next.route_switch_uroute = r.route_switch_uroute;
      next.route_probe_overflow = r.route_probe_overflow;
      next.align_pairs = align_pairs;
    }
    r.selected_iteration =
        restored_selection ? restored_selection : offset + r.iterations;
    if (keep_feasible && !pin_timing_seen && inflated_ratio == 1 &&
        !best_feasible_x.empty()) {
      x = std::move(best_feasible_x);
      y = std::move(best_feasible_y);
      r.selected_iteration = best_feasible_iteration;
    }
    d.node_x = std::move(x);
    d.node_y = std::move(y);
    d.node_x.resize(original);
    d.node_y.resize(original);
    d.node_size_x.resize(original);
    d.node_size_y.resize(original);
    d.build_maps();
    if (o.observer && (r.iterations % o.observer_interval != 0 ||
                       r.selected_iteration != offset + r.iterations)) {
      o.observer(d, PlacementStage::global_placement, offset + r.iterations);
    }
  }
  r.global_placement_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - gp_start)
          .count();
  r.overflow = overflow(d, o.bins_x, o.bins_y, o.target_density,
                        o.scale_macro_density, o.macro_transparency.always);
  if (separate_circles) {
    const auto circles = detail::centroid_circle_overlaps(d, o.target_density);
    r.centroid.circle_overlaps = circles.overlaps;
    r.centroid.max_circle_overlap = circles.maximum;
    r.centroid.circles_converged = circles.maximum <= o.centroid_stop_overlap;
  }
  // The checkpoint and objective criteria only exist once the GP loop ran;
  // otherwise only the overflow test applies.
  const bool gp_ran = d.num_movable && o.max_iterations;
  r.converged =
      (!gp_ran || !o.checkpoint_after_density || !o.checkpoint ||
       !o.checkpoint_interval ||
       (o.magnet_runtime && o.magnet_runtime->checkpoint_invoked) ||
       std::any_of(r.checkpoint_history.begin(), r.checkpoint_history.end(),
                   [](const auto &sample) { return sample.invoked; })) &&
      (!gp_ran ||
       (!o.post_density_factor && r.route_model != RouteModel::uroute &&
        !learned_ladder) ||
       r.objective_converged) &&
      r.overflow <= o.stop_overflow &&
      (!separate_circles || r.centroid.circles_converged) &&
      (!align_centroids || r.iterations > o.centroid_align_full) &&
      (!o.macro_align_weight ||
       (warm ? state->iteration : o.iteration_offset) + r.iterations >
           o.macro_align_full) &&
      (!o.magnet_runtime ||
       o.magnet_runtime->ready(o.iteration_offset + r.iterations));
  if (o.legalize && d.num_movable) {
    if (o.before_legalize) {
      o.before_legalize(d);
      d.validate();
    }
    if (o.macro_snap) {
      const auto before = d;
      const auto snap =
          o.macro_snap_regular
              ? detail::macro_regular_snap(d, o.macro_snap_options,
                                           o.macro_snap_weight)
              : detail::macro_align_snap(d, o.macro_snap_tolerance,
                                         o.macro_snap_weight);
      AlignmentReport times;
      detail::legalize(d, aux_x, aux_y, &snap, &times, o.after_macro_legalize,
                       o.macro_only_legalize);
      r.alignment = detail::macro_alignment_report(before, d, snap);
      r.alignment.macro_legalization_seconds = times.macro_legalization_seconds;
      r.alignment.cell_legalization_seconds = times.cell_legalization_seconds;
    } else {
      detail::legalize(d, aux_x, aux_y, nullptr, &r.alignment,
                       o.after_macro_legalize, o.macro_only_legalize);
    }
    std::string legality_reason;
    const bool legal =
        o.macro_only_legalize
            ? macros_are_legal(d, &legality_reason)
            : is_legal(d, &legality_reason) && detail::kernel_is_legal(d);
    if (!legal) {
      throw std::runtime_error("DREAMPlace legalization failed: insufficient "
                               "row space or unsupported geometry: " +
                               legality_reason);
    }
    if (o.observer) {
      o.observer(d, PlacementStage::legalization,
                 o.iteration_offset + r.iterations);
    }
    if (o.detailed_place && !d.net_weights.empty()) {
      const auto detail_start = std::chrono::steady_clock::now();
      detail::detailed_place(d, aux_x, aux_y, o.num_threads,
                             o.macro_snap || bool(o.after_macro_legalize),
                             o.detailed_iterations, &r.detailed_history);
      r.detailed_placement_seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        detail_start)
              .count();
      if (o.observer) {
        o.observer(d, PlacementStage::detailed_placement,
                   o.iteration_offset + r.iterations);
      }
    }
  }
  r.legal = is_legal(d);
  r.hpwl = hpwl(d);
  const auto final_charge = reservation_charge(d);
  r.reserved_area = final_charge.virtual_area;
  r.materialized_reserved_area = final_charge.materialized_area;
  r.physical_overflow = reservation_accounting
                            ? overflow(d, o.bins_x, o.bins_y, o.target_density,
                                       o.scale_macro_density, false, 1, false)
                            : r.overflow;
  for (int i = d.num_movable; i < original; ++i) {
    if (d.node_x[i] != model.node_x[i] || d.node_y[i] != model.node_y[i]) {
      throw std::runtime_error("kernel moved a fixed node");
    }
  }
  if (o.observer) {
    o.observer(d, PlacementStage::finished, o.iteration_offset + r.iterations);
  }
  if (state && next.signature) {
    *state = std::move(next);
  }
  if ((o.checkpoint && o.checkpoint_interval > 0) || o.before_legalize ||
      o.after_macro_legalize) {
    design = std::move(d);
    return r;
  }
  design.node_x = std::move(d.node_x);
  design.node_y = std::move(d.node_y);
  return r;
}
} // namespace DreamPlace::native
