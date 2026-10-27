// BSD-3-Clause. See LICENSE. Checks deliberately remain active under -c opt.
#include "dreamplace/native/placement.hpp"
#include "dreamplace/native/rudy.hpp"

#include <array>
#include <barrier>
#include <cmath>
#include <future>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

#include "dreamplace/native/spectral.hpp"
#include "utility/src/detailed_place_db.h"
#include "utility/src/utils.h"

#include "abacus_legalize/src/abacus_legalize_cpu.h"
#include "independent_set_matching/src/bin2node_3d_map.h"

namespace DreamPlace {
template <typename T> int globalSwapCPULauncher(DetailedPlaceDB<T>, int);
}

using namespace DreamPlace::native;
namespace {
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
void near(double a, double b, double eps, const char *message) {
  check(std::abs(a - b) <= eps, message);
}
template <typename F> void rejects(F f) {
  bool rejected = false;
  try {
    f();
  } catch (const std::exception &) {
    rejected = true;
  }
  check(rejected, "expected rejection");
}
void density_plateau_hands_off_without_claiming_convergence() {
  Design initial;
  initial.layout = {0, 0, 2, 2, 1, 2};
  initial.num_movable = 1;
  initial.node_x = initial.node_y = {0};
  initial.node_size_x = initial.node_size_y = {2};
  initial.flat_net2pin_start = {0};
  initial.build_maps();
  Options options;
  options.bins_x = options.bins_y = 8;
  options.max_iterations = 400;
  options.min_iterations = 200;
  options.target_density = .5;
  options.stop_overflow = .1;
  options.post_density_factor = 4;
  options.adaptive_density = true;
  options.add_fillers = false;
  options.detailed_place = false;
  // One legal cell fills the entire core: no movement can satisfy 50% density.
  auto cells = initial;
  const auto stalled = place(cells, options);
  check(stalled.density_stalled && !stalled.converged && stalled.legal,
        "density plateau did not hand off to legal placement honestly");
  check(stalled.iterations >= options.min_iterations &&
            stalled.iterations < options.max_iterations,
        "density plateau ignored its minimum or exhausted its maximum");
  auto fixed_budget = initial;
  options.min_iterations = options.max_iterations;
  const auto full = place(fixed_budget, options);
  check(full.iterations == options.max_iterations,
        "density plateau shortened an explicit minimum iteration budget");
  auto macros = initial;
  macros.node_macro = {1};
  options.min_iterations = 200;
  const auto macro_result = place(macros, options);
  check(!macro_result.density_stalled &&
            macro_result.iterations == options.max_iterations,
        "cell plateau policy interrupted a macro solve");
}

void setup_timing_does_not_veto_convergence() {
  Design initial;
  initial.layout = {0, 0, 100, 20, 1, 1};
  initial.num_movable = 2;
  initial.node_x = {40, 50, 0, 100};
  initial.node_y = {10, 10, 10, 10};
  initial.node_size_x = initial.node_size_y = {1, 1, 0, 0};
  initial.pin2node = {2, 0, 0, 1, 1, 3};
  initial.pin_offset_x = initial.pin_offset_y = {0, 0, 0, 0, 0, 0};
  initial.flat_net2pin = {0, 1, 2, 3, 4, 5};
  initial.flat_net2pin_start = {0, 2, 4, 6};
  initial.build_maps();
  Options options;
  options.bins_x = options.bins_y = 8;
  options.max_iterations = 400;
  options.min_iterations = 20;
  options.stop_overflow = .1;
  options.post_density_factor = 4;
  options.convergence_tolerance = 1;
  options.checkpoint_interval = 20;
  options.legalize = options.detailed_place = false;
  auto positive = initial;
  options.checkpoint = [](const Design &, int) {
    GpUpdate update;
    update.timing_violation = 0;
    return update;
  };
  const auto reference = place(positive, options);
  auto negative = initial;
  options.checkpoint = [](const Design &, int) {
    GpUpdate update;
    update.timing_violation = 1000;
    return update;
  };
  const auto failing_timing = place(negative, options);
  check(reference.converged && failing_timing.converged,
        "setup timing vetoed feasible placement");
  check(reference.iterations < options.max_iterations &&
            reference.iterations == failing_timing.iterations,
        "setup timing extended the GP iteration budget");
  check(positive.node_x == negative.node_x &&
            positive.node_y == negative.node_y,
        "timing-only observation changed placement");
  check(failing_timing.timing_observed && !failing_timing.timing_satisfied,
        "convergence erased a failing timing observation");
  auto seeded = initial;
  auto seed_options = options;
  seed_options.max_iterations = seed_options.min_iterations = 1;
  seed_options.density_weight_seed = 3;
  const auto seed_result = place(seeded, seed_options);
  near(seed_result.density_weight, 3, 1e-12,
       "refinement renormalized its density seed");
  seed_options.density_weight_seed = 0;
  rejects([&] { place(seeded, seed_options); });
  // Milestones coalesce when crossed together, and never repeat just because
  // the sample clock ticks or the iteration cap is reached.
  options.min_iterations = options.max_iterations = 100;
  options.checkpoint_overflow_milestones = {1, .5};
  options.checkpoint_final = true;
  unsigned evaluations = 0;
  options.checkpoint = [&](const Design &, int) {
    ++evaluations;
    return GpUpdate{};
  };
  auto milestones = initial;
  place(milestones, options);
  check(evaluations == 1, "crossed density milestones repeated timing");
  options.checkpoint_overflow_milestones = {0};
  options.target_density = .00001;
  evaluations = 0;
  auto infeasible = initial;
  place(infeasible, options);
  check(evaluations == 0,
        "iteration cap forced timing on an overflowing placement");
  options.checkpoint_overflow_milestones = {.1, .2};
  rejects([&] { place(infeasible, options); });
}
void placement_masks() {
  // A full-height placement-only strip must not sever an estimated signal.
  Design d;
  d.layout = {0, 0, 20, 4, 1, 1};
  d.node_x = {9, 0, 20};
  d.node_y = {0, 2, 2};
  d.node_size_x = {2, 0, 0};
  d.node_size_y = {4, 0, 0};
  d.node_macro = {1, 0, 0};
  d.pin2node = {1, 2};
  d.pin_offset_x = d.pin_offset_y = {0, 0};
  d.flat_net2pin = {0, 1};
  d.flat_net2pin_start = {0, 2};
  d.net_weights = {1};
  d.build_maps();
  RudyOptions route;
  route.bins_x = 10;
  route.bins_y = 2;
  route.macro_blockages = route.pin_access = true;
  route.model = RudyModel::rsmt;
  check(rudy(d, route).unroutable_nets == 1,
        "real macro wall must block routing");
  for (auto kind : {PlacementMask::hard, PlacementMask::soft}) {
    d.node_placement_mask = {kind, PlacementMask::none, PlacementMask::none};
    const auto result = rudy(d, route);
    check(result.unroutable_nets == 0 && result.blocked_pins == 0,
          "placement-only wall must allow routing");
    for (double f : result.blocked_fraction)
      check(f == 0, "mask reduced routing capacity");
  }
  auto bad = d;
  bad.node_placement_mask[1] = PlacementMask::hard;
  rejects([&] { bad.validate(); }); // A pin-bearing object cannot be a mask.
  bad = d;
  bad.node_placement_mask.resize(1);
  rejects([&] { bad.validate(); });
  // A repair buffer has space only in the soft strip between two hard masks.
  Design cells;
  cells.layout = {0, 0, 20, 4, 1, 1};
  cells.num_movable = 1;
  cells.node_x = {9, 0, 8, 12};
  cells.node_y = {1, 0, 0, 0};
  cells.node_size_x = {2, 8, 4, 8};
  cells.node_size_y = {1, 4, 4, 4};
  cells.node_macro = {0, 1, 1, 1};
  cells.node_placement_mask = {PlacementMask::none, PlacementMask::hard,
                               PlacementMask::soft, PlacementMask::hard};
  cells.build_maps();
  check(is_legal(cells), "buffer in a soft channel should be legal");
  check(detail::kernel_is_legal(cells), "row checker must release soft masks");
  const auto sizes = cells.node_size_x;
  Options options;
  options.max_iterations = 0;
  options.routability = false;
  options.bins_x = options.bins_y = 2;
  check(place(cells, options).legal,
        "buffer legalization in a soft channel failed");
  check(cells.node_x[0] >= 8 && cells.node_x[0] + 2 <= 12,
        "buffer crossed a hard placement exclusion");
  check(cells.node_size_x == sizes, "legalization changed mask geometry");
  auto occupied = cells;
  occupied.node_placement_mask[2] = PlacementMask::none;
  check(!is_legal(occupied), "real macro overlap was allowed");
  occupied = cells;
  occupied.node_placement_mask[2] = PlacementMask::hard;
  check(!is_legal(occupied), "hard mask overlap was allowed");
  // Both masks retain identical GP density charge despite different legality.
  near(density_overflow(cells, 2, 2, .5), density_overflow(occupied, 2, 2, .5),
       1e-12, "soft mask lost its GP density charge");
}

Design design() {
  Design d;
  d.layout = {10, 20, 50, 30, 1, 1};
  d.num_movable = 6;
  d.node_x = {30, 30, 30, 30, 30, 30, 25, 10, 50};
  d.node_y = {25, 25, 25, 25, 25, 25, 24, 25, 25};
  d.node_size_x = {2, 2, 2, 2, 2, 2, 4, 0, 0};
  d.node_size_y = {1, 1, 1, 1, 1, 1, 2, 0, 0};
  d.pin2node = {7, 0, 0, 1, 2, 1, 3, 2, 3, 4, 4, 5, 5, 8};
  d.pin_offset_x.assign(14, 0.5);
  d.pin_offset_y.assign(14, 0.25);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  d.flat_net2pin_start = {0, 2, 5, 8, 10, 12, 14};
  d.build_maps();
  return d;
}
void spectral() {
  const int nx = 7, ny = 6;
  const double w = 14, h = 9;
  PoissonSolver s(nx, ny, w, h);
  std::vector<double> rho(nx * ny, 1), phi, fx, fy;
  s.solve(rho, phi, fx, fy);
  for (int i = 0; i < nx * ny; ++i) {
    near(phi[i], 0, 1e-12, "constant charge potential");
    near(fx[i], 0, 1e-12, "constant field");
  }
  double wx = std::numbers::pi * 2 / w, wy = std::numbers::pi * 3 / h;
  for (int x = 0; x < nx; ++x)
    for (int y = 0; y < ny; ++y)
      rho[x * ny + y] = std::cos(std::numbers::pi * 2 * (x + 0.5) / nx) *
                        std::cos(std::numbers::pi * 3 * (y + 0.5) / ny);
  s.solve(rho, phi, fx, fy);
  for (int x = 0; x < nx; ++x)
    for (int y = 0; y < ny; ++y) {
      double ax = std::numbers::pi * 2 * (x + 0.5) / nx,
             ay = std::numbers::pi * 3 * (y + 0.5) / ny;
      near(phi[x * ny + y], rho[x * ny + y] / (wx * wx + wy * wy), 1e-12,
           "Poisson eigenmode");
      near(fx[x * ny + y],
           wx * std::sin(ax) * std::cos(ay) / (wx * wx + wy * wy), 1e-12,
           "Poisson x field");
      near(fy[x * ny + y],
           wy * std::cos(ax) * std::sin(ay) / (wx * wx + wy * wy), 1e-12,
           "Poisson y field");
    }
}
void gradients() {
  auto d = design();
  d.net_weights = {0, 2, 0.5, 4, 1, 3};
  d.node_x[1] += 2.4;
  d.node_y[3] -= 0.7;
  auto g = detail::wirelength(d, 1.7, 2);
  for (int i = 0; i < d.num_movable; ++i)
    for (int axis = 0; axis < 2; ++axis) {
      auto &pos = axis ? d.node_y : d.node_x;
      double old = pos[i], eps = 1e-5;
      pos[i] = old + eps;
      double plus = detail::wirelength(d, 1.7, 1).value;
      pos[i] = old - eps;
      double minus = detail::wirelength(d, 1.7, 1).value;
      pos[i] = old;
      near((plus - minus) / (2 * eps), axis ? g.y[i] : g.x[i], 2e-7,
           "weighted WA gradient");
    }
  // Large absolute coordinates must not cause exp() overflow.
  for (auto &x : d.node_x)
    x += 1e6;
  check(std::isfinite(detail::wirelength(d, 0.1, 1).value), "stable WA kernel");
}
void native_flow() {
  auto d = design(), before = d;
  check(!is_legal(d), "overlapping input should be illegal");
  Options o;
  o.bins_x = 8;
  o.bins_y = 8;
  o.max_iterations = 80;
  o.num_threads = 1;
  auto r = place(d, o);
  check(r.legal && is_legal(d), "native GP/LG/DP must legalize");
  check(r.iterations > 0 && r.iterations <= o.max_iterations,
        "iteration accounting");
  check(std::isfinite(r.hpwl) && std::isfinite(r.overflow),
        "finite placement metrics");
  check(d.node_x.size() == before.node_x.size() &&
            d.pin2node == before.pin2node,
        "fillers must not escape");
  for (size_t i = d.num_movable; i < d.node_x.size(); ++i) {
    check(d.node_x[i] == before.node_x[i] && d.node_y[i] == before.node_y[i],
          "fixed obstacle/terminal moved");
  }
  auto again = before;
  auto rr = place(again, o);
  check(d.node_x == again.node_x && d.node_y == again.node_y &&
            rr.hpwl == r.hpwl,
        "deterministic placement");
  auto observed = before;
  Options traced = o;
  traced.observer_interval =
      7; // Deliberately not aligned with 20-step convergence checks.
  std::vector<PlacementStage> stages;
  Design last;
  int previous_iteration = 0;
  traced.observer = [&](const Design &frame, PlacementStage stage,
                        int iteration) {
    frame.validate();
    check(frame.node_x.size() == before.node_x.size() &&
              frame.pin2node == before.pin2node,
          "observer must see physical nodes/connectivity without fillers");
    for (size_t i = frame.num_movable; i < frame.node_x.size(); ++i)
      check(frame.node_x[i] == before.node_x[i] &&
                frame.node_y[i] == before.node_y[i],
            "observer fixed nodes moved");
    check(iteration >= previous_iteration, "observer iteration order");
    previous_iteration = iteration;
    if (stage == PlacementStage::initial)
      check(frame.node_x == before.node_x && frame.node_y == before.node_y &&
                iteration == 0,
            "observer initial state");
    stages.push_back(stage);
    last = frame;
  };
  auto tr = place(observed, traced);
  check(observed.node_x == d.node_x && observed.node_y == d.node_y &&
            tr.hpwl == r.hpwl && tr.iterations == r.iterations,
        "observation must not change placement");
  check(stages.front() == PlacementStage::initial &&
            stages.back() == PlacementStage::finished &&
            stages[stages.size() - 2] == PlacementStage::detailed_placement &&
            stages[stages.size() - 3] == PlacementStage::legalization &&
            stages[stages.size() - 4] == PlacementStage::global_placement,
        "observer phase sequence");
  check(last.node_x == observed.node_x && last.node_y == observed.node_y,
        "observer final state");
  traced.observer = [](const Design &, PlacementStage stage, int) {
    if (stage == PlacementStage::finished)
      throw std::runtime_error("observer abort");
  };
  auto aborted = before;
  rejects([&] { place(aborted, traced); });
  check(aborted.node_x == before.node_x && aborted.node_y == before.node_y,
        "observer abort must not write back");
  traced.observer_interval = 0;
  rejects([&] { place(aborted, traced); });
  auto parallel = before;
  o.num_threads = 2;
  check(place(parallel, o).legal, "parallel native flow");
  auto gp = before;
  o.legalize = false;
  o.detailed_place = false;
  o.add_fillers = false;
  o.max_iterations = 30;
  auto gr = place(gp, o);
  check(gr.iterations == 30, "GP-only bounded iteration limit");
  check(gp.node_x != before.node_x || gp.node_y != before.node_y,
        "GP should move cells");
  auto lg = before;
  o.legalize = true;
  o.detailed_place = true;
  o.max_iterations = 0;
  check(place(lg, o).legal, "legalization-only path");
}
void run_global_swap(Design &d, int bx, int by) {
  auto ix = d.node_x, iy = d.node_y;
  std::vector<int> regions(d.node_x.size(), 0);
  DreamPlace::DetailedPlaceDB<double> db{};
  db.init_x = ix.data();
  db.init_y = iy.data();
  db.x = d.node_x.data();
  db.y = d.node_y.data();
  db.node_size_x = d.node_size_x.data();
  db.node_size_y = d.node_size_y.data();
  db.node2fence_region_map = regions.data();
  db.flat_net2pin_map = d.flat_net2pin.data();
  db.flat_net2pin_start_map = d.flat_net2pin_start.data();
  db.flat_node2pin_map = d.flat_node2pin.data();
  db.flat_node2pin_start_map = d.flat_node2pin_start.data();
  db.pin2node_map = d.pin2node.data();
  db.pin2net_map = d.pin2net.data();
  db.pin_offset_x = d.pin_offset_x.data();
  db.pin_offset_y = d.pin_offset_y.data();
  db.net_mask = d.net_mask.data();
  db.net_weights = d.net_weights.data();
  db.xl = d.layout.xl;
  db.yl = d.layout.yl;
  db.xh = d.layout.xh;
  db.yh = d.layout.yh;
  db.site_width = d.layout.site_width;
  db.row_height = d.layout.row_height;
  db.num_bins_x = bx;
  db.num_bins_y = by;
  db.bin_size_x = (db.xh - db.xl) / bx;
  db.bin_size_y = (db.yh - db.yl) / by;
  db.num_sites_x = std::round((db.xh - db.xl) / db.site_width);
  db.num_sites_y = std::round((db.yh - db.yl) / db.row_height);
  db.num_nodes = d.node_x.size();
  db.num_movable_nodes = d.num_movable;
  db.num_pins = d.pin2node.size();
  db.num_nets = d.net_weights.size();
  DreamPlace::globalSwapCPULauncher(db, 2);
}

void offgrid_macro_swap() {
  Design d;
  d.layout = {0, 0, 16, 2, 1, 1};
  d.num_movable = 2;
  d.node_x = {10, 10, 7.2, 10, 15};
  d.node_y = {0, 1, 1, 1, 0};
  d.node_size_x = {3, 1, 2.3, 0, 0};
  d.node_size_y = {1, 1, 1, 0, 0};
  d.pin2node = {0, 3, 1, 4};
  d.pin_offset_x.assign(4, 0);
  d.pin_offset_y.assign(4, 0);
  d.flat_net2pin = {0, 1, 2, 3};
  d.flat_net2pin_start = {0, 2, 4};
  d.net_weights = {1, 1};
  d.build_maps();
  check(is_legal(d), "off-grid obstacle fixture starts legal");
  const double before = hpwl(d);
  run_global_swap(d, 4, 2);
  check(is_legal(d), "global swap clamped a cell to an off-grid macro edge");
  check(hpwl(d) < before, "legal cross-row swap must improve wirelength");
  near(d.node_x[2], 7.2, 0, "fixed macro must not move");
}

void same_bin_swap() {
  Design d;
  d.layout = {0, 0, 12, 1, 1, 1};
  d.num_movable = 2;
  d.node_x = {2, 4, 5, 1};
  d.node_y.assign(4, 0);
  d.node_size_x = {1, 1, 0, 0};
  d.node_size_y = {1, 1, 0, 0};
  d.pin2node = {0, 2, 1, 3};
  d.pin_offset_x.assign(4, 0);
  d.pin_offset_y.assign(4, 0);
  d.flat_net2pin = {0, 1, 2, 3};
  d.flat_net2pin_start = {0, 2, 4};
  d.net_weights = {1, 1};
  d.build_maps();
  const auto before = hpwl(d);
  run_global_swap(d, 2, 1);
  check(is_legal(d), "same-bin swap must preserve legality");
  check(hpwl(d) < before, "global swap skipped cells contained in one bin");
}

void weighted_detailed_placement() {
  Design d;
  d.layout = {0, 0, 12, 1, 1, 1};
  d.num_movable = 2;
  d.node_x = {1, 9, 0, 12};
  d.node_y = {0, 0, 0, 0};
  d.node_size_x = {1, 1, 0, 0};
  d.node_size_y = {1, 1, 0, 0};
  d.pin2node = {0, 3, 0, 2, 1, 2, 1, 3};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.net_weights = {10, 1, 10, 1};
  d.build_maps();
  const auto x = d.node_x, y = d.node_y;
  const double before = hpwl(d);
  auto explicit_default = d;
  std::vector<DetailedSample> samples;
  detail::detailed_place(explicit_default, 4, 1, 1, false, 2, &samples);
  check(samples.size() == 3, "one diagnostic per detailed-placement kernel");
  for (const auto &sample : samples)
    check(!sample.accepted || (sample.legal && sample.macros_preserved &&
                               sample.after <= sample.before + 1e-8),
          "accepted detailed moves must preserve legality and objective");
  auto extended = d;
  detail::detailed_place(extended, 4, 1, 1, false, 10);
  check(is_legal(extended) && hpwl(extended) <= before,
        "extended weighted detail must preserve legality and objective");
  detail::detailed_place(d, 4, 1, 1);
  check(d.node_x == explicit_default.node_x &&
            d.node_y == explicit_default.node_y,
        "default detailed iteration limit changed");
  check(is_legal(d), "weighted detailed placement must remain legal");
  check(hpwl(d) < before,
        "weighted detailed placement must optimize active timing weights");
  check(d.node_x[0] > d.node_x[1],
        "critical connections must beat the unweighted tie");
  check(d.node_x[2] == x[2] && d.node_x[3] == x[3] && d.node_y == y,
        "weighted detail must preserve fixed terminal locations");
  auto overlapping = d;
  overlapping.node_x[1] = overlapping.node_x[0];
  std::string reason;
  check(!is_legal(overlapping, &reason) &&
            reason.find("overlap") != std::string::npos,
        "legality diagnostic must identify the colliding nodes");
}
// A DBU-scale die with large timing weights must not overflow the int cost
// matrices of independent-set matching (large_number ~ perimeter*set_size).
void weighted_detailed_placement_wide_die() {
  const double w = 8e6;
  Design d;
  // Matching only solves sets of >= 3 equal-size cells sharing a bin, and
  // (upstream quirk) needs two such sets to run at all. Cells 0-3 are the
  // 1-wide set that must swap; cells 4-6 are a 2-wide filler set. One cell per
  // row keeps k_reorder from doing the swap first.
  d.layout = {0, 0, w, 7, 1, 1};
  d.num_movable = 7;
  d.node_x = {1, 3, 5, 7, 1, 4, 7, 0, w};
  d.node_y = {0, 1, 2, 3, 4, 5, 6, 0, 0};
  d.node_size_x = {1, 1, 1, 1, 2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 1, 1, 1, 1, 0, 0};
  // net0: cell0-right (heavy), net1: cell3-left (heavy), net2: cell1-left,
  // net3: cell2-right
  d.pin2node = {0, 8, 3, 7, 1, 7, 2, 8};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.net_weights = {1000, 1000, 1,
                   1}; // weight * die span exceeds INT_MAX unscaled
  d.build_maps();
  const double before = hpwl(d);
  std::vector<DetailedSample> samples;
  detail::detailed_place(d, 16, 1, 1, false, 2, &samples);
  for (const auto &sample : samples)
    check(!sample.accepted || (sample.legal && sample.after <= sample.before),
          "wide-die weighted detail must stay legal and not worsen HPWL");
  check(is_legal(d) && hpwl(d) <= before,
        "wide-die weighted detail must preserve legality and objective");
  // The int-overflow guard is only observable if matching actually ran and won.
  const DetailedSample *matching = nullptr;
  for (const auto &sample : samples)
    if (sample.kernel == "independent_set_matching")
      matching = &sample;
  check(matching && matching->accepted && matching->after < matching->before,
        "wide-die independent-set matching must be accepted and improve HPWL");
  check(hpwl(d) < before && d.node_x[0] > d.node_x[3],
        "wide-die critical connections must beat the unweighted tie");
}
void placer_instances() {
  Options initial;
  initial.bins_x = initial.bins_y = 8;
  initial.max_iterations = 17;
  initial.legalize = false;
  initial.detailed_place = false;
  Placer global(initial), legalizer;
  legalizer.options().max_iterations = 0;
  // Construction owns a configuration copy, not a reference to the caller's.
  initial.max_iterations = 2;
  const Placer &configured = global;
  auto d = design();
  check(configured.place(d).iterations == 17,
        "placer instance must own its options");
  check(configured.options().max_iterations == 17 &&
            legalizer.options().legalize &&
            legalizer.options().max_iterations == 0,
        "separate placer instances must have independent configuration");
  check(legalizer.place(d).legal && is_legal(d),
        "instance must support legalizing a caller-owned placement");
  // Reconfigure and reuse one object; failed runs must preserve caller data.
  global.options().max_iterations = -1;
  auto before = d;
  rejects([&] { global.place(d); });
  check(d.node_x == before.node_x && d.node_y == before.node_y,
        "invalid instance options must not change design");
  global.options() = legalizer.options();
  check(global.place(d).legal,
        "placer instance must allow reconfiguration and reuse");
  check(legalizer.options().max_iterations == 0,
        "reconfiguring another instance must not change its source");
}
void concurrent_instances() {
  struct Run {
    Design design;
    Result result;
  };
  constexpr size_t count = 4;
  std::array<Design, count> inputs;
  std::array<Options, count> options;
  std::array<Run, count> reference;
  for (size_t i = 0; i < count; ++i) {
    inputs[i] = design();
    auto &d = inputs[i];
    double offset = 100 * i;
    d.layout.xl += offset;
    d.layout.xh += offset;
    for (auto &x : d.node_x)
      x += offset;
    options[i].bins_x = options[i].bins_y = 8 + 2 * i;
    options[i].max_iterations = 80 + 20 * i;
  }
  // Include the LEMON macro LP alongside standard-cell GP/LG/DP jobs.
  inputs.back().node_size_x[0] = 6;
  inputs.back().node_size_y[0] = 3;
  inputs.back().node_size_x[1] = 5;
  inputs.back().node_size_y[1] = 2;
  options.back().max_iterations = 0;
  options.back().detailed_place = false;
  for (size_t i = 0; i < count; ++i) {
    auto d = inputs[i];
    auto result = Placer(options[i]).place(d);
    check(result.legal, "serial concurrency reference must be legal");
    reference[i] = {std::move(d), result};
  }
  // First check exact isolation with one worker per instance. Then exercise
  // concurrent OpenMP teams, where reduction ordering may change coordinates.
  for (int threads : {1, 2}) {
    std::barrier start(static_cast<std::ptrdiff_t>(count));
    std::array<std::future<Run>, count> jobs;
    for (size_t i = 0; i < count; ++i)
      jobs[i] = std::async(std::launch::async, [&, i] {
        Placer placer(options[i]);
        placer.options().num_threads = threads;
        placer.options().observer = [&](const Design &frame,
                                        PlacementStage stage, int) {
          if (stage == PlacementStage::initial)
            start.arrive_and_wait();
          check(frame.node_x.size() == inputs[i].node_x.size(),
                "parallel snapshot leaked another run's nodes");
          for (size_t j = frame.num_movable; j < frame.node_x.size(); ++j)
            check(frame.node_x[j] == inputs[i].node_x[j] &&
                      frame.node_y[j] == inputs[i].node_y[j],
                  "parallel run moved another instance's fixed node");
        };
        auto d = inputs[i];
        auto result = placer.place(d);
        return Run{std::move(d), result};
      });
    for (size_t i = 0; i < count; ++i) {
      auto run = jobs[i].get();
      check(run.result.legal && is_legal(run.design) &&
                std::isfinite(run.result.hpwl),
            "parallel instances must produce legal finite placements");
      if (threads == 1)
        check(run.design.node_x == reference[i].design.node_x &&
                  run.design.node_y == reference[i].design.node_y &&
                  run.result.hpwl == reference[i].result.hpwl &&
                  run.result.iterations == reference[i].result.iterations,
              "independent parallel calls must match their serial references");
    }
  }
}
// Floating-point atomic accumulation into the density map makes the result
// depend on thread arrival order. Repeated identical runs must agree exactly
// for a fixed worker count, including when OpenMP is enabled (--config=openmp).
void threaded_determinism() {
  Design d;
  d.layout = {0, 0, 64, 32, 1, 1};
  const int cells = 300, nets = 380;
  d.num_movable = cells;
  unsigned state = 12345;
  const auto next = [&](unsigned bound) {
    state = state * 1664525u + 1013904223u;
    return (state >> 8) % bound;
  };
  for (int i = 0; i < cells; ++i) {
    d.node_x.push_back(28 + next(8));
    d.node_y.push_back(14 + next(4));
    d.node_size_x.push_back(2);
    d.node_size_y.push_back(1);
  }
  d.flat_net2pin_start = {0};
  for (int n = 0; n < nets; ++n) {
    const int degree = 2 + next(3);
    for (int k = 0; k < degree; ++k) {
      d.flat_net2pin.push_back(d.pin2node.size());
      d.pin2node.push_back(next(cells));
      d.pin_offset_x.push_back(next(3) * .5);
      d.pin_offset_y.push_back(.5);
    }
    d.flat_net2pin_start.push_back(d.flat_net2pin.size());
  }
  // A cell may repeat inside a net; build_maps accepts repeated nodes.
  d.build_maps();
  Options o;
  o.bins_x = o.bins_y = 16;
  o.max_iterations = 120;
  o.num_threads = 4;
  o.legalize = false;
  o.detailed_place = false;
  auto first = d;
  place(first, o);
  for (int run = 0; run < 6; ++run) {
    auto again = d;
    place(again, o);
    check(again.node_x == first.node_x && again.node_y == first.node_y,
          "multithreaded global placement is nondeterministic");
  }
  // The legalizer and detailed-placement kernels also take the worker count.
  o.legalize = o.detailed_place = true;
  auto full = d;
  place(full, o);
  for (int run = 0; run < 6; ++run) {
    auto again = d;
    place(again, o);
    check(again.node_x == full.node_x && again.node_y == full.node_y,
          "multithreaded legalization/detailed placement is nondeterministic");
  }
}
void hpwl_case() {
  Design d;
  d.layout = {0, 0, 20, 10, 1, 1};
  d.num_movable = 1;
  d.node_x = {2, 10};
  d.node_y = {3, 8};
  d.node_size_x = {1, 0};
  d.node_size_y = {1, 0};
  d.pin2node = {0, 1};
  d.pin_offset_x = {1, 0};
  d.pin_offset_y = {0, 0};
  d.flat_net2pin = {1, 0};
  d.flat_net2pin_start = {0, 2};
  d.net_weights = {2};
  d.build_maps();
  near(hpwl(d), 24, 0, "weighted HPWL including pin offsets");
  d.net_mask[0] = 0;
  near(hpwl(d), 0, 0, "masked net");
  // Real pins with no net, duplicate pins, nonmonotonic CSR and bad maps fail
  // before kernels.
  auto bad = d;
  bad.flat_net2pin = {1, 1};
  rejects([&] { bad.build_maps(); });
  bad = d;
  bad.pin2node[0] = 20;
  rejects([&] { bad.build_maps(); });
  bad = d;
  bad.pin2net[0] = 9;
  rejects([&] { bad.validate(); });
  bad = d;
  bad.node_x[0] = std::numeric_limits<double>::quiet_NaN();
  rejects([&] { place(bad); });
  bad = d;
  bad.flat_node2pin_start.back() = 1;
  rejects([&] { place(bad); });
  auto original = d;
  Options o;
  o.target_density = 0;
  rejects([&] { place(d, o); });
  check(d.node_x == original.node_x, "invalid options mutated design");
  o = Options{};
  o.detailed_iterations = 0;
  rejects([&] { place(d, o); });
  o = Options{};
  o.wirelength_gamma_end = 0;
  rejects([&] { place(d, o); });
}
void macros() {
  Design d;
  d.layout = {0, 0, 40, 20, 1, 1};
  d.num_movable = 4;
  d.node_x = {10, 10, 10, 10, 10};
  d.node_y = {5, 5, 5, 5, 5};
  d.node_size_x = {8, 6, 2, 2, 4};
  d.node_size_y = {6, 7, 1, 1, 8};
  d.build_maps();
  auto original = d;
  Options o;
  o.max_iterations = 0;
  o.bins_x = 8;
  o.bins_y = 8;
  o.detailed_place = false;
  check(place(d, o).legal && is_legal(d), "Hannan/LP macro legalization");
  check(d.node_x.back() == original.node_x.back() &&
            d.node_y.back() == original.node_y.back(),
        "macro moved fixed obstacle");
  auto shifted = original;
  shifted.layout = {-50, -40, -10, -20, 1, 1};
  for (auto &x : shifted.node_x)
    x -= 50;
  for (auto &y : shifted.node_y)
    y -= 40;
  check(place(shifted, o).legal,
        "macro legalization with negative core origin");
}
template <typename T> void zero_area_fixed_terminals() {
  const T x[] = {T(3.5), T(5.5), T(10), T(7), T(8)};
  const T y[] = {T(0), T(.5), T(0), T(.5), T(.5)};
  const T w[] = {T(2), T(0), T(1), T(0), T(1)};
  const T h[] = {T(1), T(0), T(1), T(1), T(0)};
  std::vector<std::vector<int>> bins(2);
  DreamPlace::distributeMovableAndFixedCells2BinsCPU(
      x, y, w, h, T(20), T(1), T(0), T(0), T(20), T(2), T(1), 1, 2, 5, 1, bins);
  check(bins[0] == std::vector<int>({0, 2}) && bins[1].empty(),
        "zero-area fixed terminals must not become Abacus row obstacles");
  DreamPlace::sortNodesInRow(x, y, w, h, 1, bins[0]);
  DreamPlace::DetailedPlaceDB<T> db{};
  db.num_nodes = 5;
  db.num_movable_nodes = 1;
  db.num_sites_y = 2;
  db.yl = 0;
  db.row_height = 1;
  db.node_size_x = w;
  db.node_size_y = h;
  std::vector<std::vector<int>> rows(2);
  db.make_row2node_map(x, y, rows, 1);
  check(rows[0] == std::vector<int>({0, 2}) && rows[1].empty(),
        "zero-area fixed terminals must not truncate detailed-placement row "
        "space");
}
void fractional_io_detailed_placement() {
  auto d = design();
  d.node_x.back() = 35.5;
  d.node_y.back() = 25.5;
  Options o;
  o.max_iterations = 0;
  o.add_fillers = false;
  const auto result = place(d, o);
  check(result.legal && is_legal(d),
        "fractional IO terminal broke detailed placement");
  near(d.node_x.back(), 35.5, 0, "detailed placement moved fixed IO x");
  near(d.node_y.back(), 25.5, 0, "detailed placement moved fixed IO y");
}
void continued_fixed_terminals() {
  auto d = design();
  Options o;
  o.max_iterations = 2;
  o.bins_x = o.bins_y = 8;
  o.add_fillers = o.legalize = o.detailed_place = false;
  GpState state;
  place(d, o, &state);
  // H2D refreshes fixed pseudo-terminals between continued GP rounds. The
  // new coordinate cannot be recovered exactly by old + (new - old).
  d.node_x.back() = 10.001;
  d.node_y.back() = 20.003;
  const auto before = d;
  place(d, o, &state);
  for (size_t i = d.num_movable; i < d.node_x.size(); ++i)
    check(d.node_x[i] == before.node_x[i] && d.node_y[i] == before.node_y[i] &&
              state.x[i] == before.node_x[i] && state.y[i] == before.node_y[i],
          "continued GP moved a refreshed fixed terminal");
}
void mixed_sizes() {
  // Integer DBU, with macro and cell sizes not divisible by site/row pitch.
  Design original;
  original.layout = {0, 0, 400, 200, 10, 10};
  original.num_movable = 5;
  original.node_x = {100, 100, 100, 100, 100, 100, 0};
  original.node_y = {50, 50, 50, 50, 50, 50, 0};
  original.node_size_x = {83, 67, 13, 17, 20, 40, 0};
  original.node_size_y = {63, 77, 23, 10, 10, 80, 0};
  original.pin2node = {6, 0, 0, 1, 1, 2, 2, 3, 3, 4};
  original.pin_offset_x.assign(10, 1);
  original.pin_offset_y.assign(10, 2);
  original.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  original.flat_net2pin_start = {0, 2, 4, 6, 8, 10};
  original.build_maps();
  for (int iterations : {0, 20}) {
    for (bool detailed : {false, true}) {
      auto d = original;
      Options o;
      o.max_iterations = iterations;
      o.bins_x = o.bins_y = 8;
      o.add_fillers = false;
      o.detailed_place = detailed;
      o.observer = [&](const Design &frame, PlacementStage, int) {
        check(frame.node_size_x == original.node_size_x &&
                  frame.node_size_y == original.node_size_y &&
                  frame.pin_offset_x == original.pin_offset_x &&
                  frame.pin_offset_y == original.pin_offset_y,
              "mixed-size solver changed physical sizes or pin offsets");
      };
      const auto result = place(d, o);
      check(result.legal && is_legal(d) && detail::kernel_is_legal(d),
            "nonintegral mixed-size legalization");
      check(d.node_size_x == original.node_size_x &&
                d.node_size_y == original.node_size_y,
            "mixed-size writeback changed master dimensions");
      for (size_t i = d.num_movable; i < d.node_x.size(); ++i)
        check(d.node_x[i] == original.node_x[i] &&
                  d.node_y[i] == original.node_y[i],
              "mixed-size placement moved a fixed object");
    }
  }
}
void relocated_macro_obstacles() {
  Design macro;
  macro.layout = {0, 0, 120, 80, 10, 10};
  macro.num_movable = 1;
  macro.node_x = {40, 40};
  macro.node_y = {0, 0};
  macro.node_size_x = {31, 40};
  macro.node_size_y = {61, 80};
  macro.build_maps();
  auto preview = macro;
  detail::legalize(preview, 4, 4);
  check(is_legal(preview) && preview.node_x[0] != macro.node_x[0],
        "macro fixture must relocate");
  // Put five row cells inside the macro's new position. Greedy legalization
  // must reserve that location, rather than the macro's obsolete GP anchor.
  Design d;
  d.layout = macro.layout;
  d.num_movable = 6;
  d.node_x = {40};
  d.node_y = {0};
  d.node_size_x = {31};
  d.node_size_y = {61};
  for (int i = 0; i < 5; ++i) {
    d.node_x.push_back(preview.node_x[0]);
    d.node_y.push_back(10 * i);
    d.node_size_x.push_back(10);
    d.node_size_y.push_back(10);
  }
  d.node_x.push_back(40);
  d.node_y.push_back(0);
  d.node_size_x.push_back(40);
  d.node_size_y.push_back(80);
  d.build_maps();
  Options o;
  o.max_iterations = 0;
  o.bins_x = o.bins_y = 4;
  o.add_fillers = false;
  o.detailed_place = false;
  check(place(d, o).legal && is_legal(d),
        "row cells must avoid the relocated macro");
  check(d.node_x.back() == 40 && d.node_y.back() == 0, "fixed obstacle moved");
}
void failures_and_empty() {
  Design empty;
  empty.layout = {0, 0, 20, 10, 1, 1};
  empty.build_maps();
  check(place(empty).legal, "empty design");
  Design d;
  d.layout = {0, 0, 2, 1, 1, 1};
  d.num_movable = 2;
  d.node_x = {0, 0};
  d.node_y = {0, 0};
  d.node_size_x = {2, 2};
  d.node_size_y = {1, 1};
  d.build_maps();
  auto before = d;
  Options o;
  o.max_iterations = 0;
  rejects([&] { place(d, o); });
  check(d.node_x == before.node_x && d.node_y == before.node_y,
        "failed legalization mutated design");
  d = design();
  d.node_x[0] = 10;
  d.node_y[0] = 20;
  d.node_x[1] = 12;
  d.node_y[1] = 20;
  // Independent legality checker distinguishes touching edges from overlap.
  for (int i = 2; i < 6; ++i) {
    d.node_x[i] = 10 + 2 * i;
    d.node_y[i] = 20;
  }
  check(is_legal(d), "touching edges are legal");
  d.node_x[1] -= 0.5;
  check(!is_legal(d), "off-site/overlap missed");
}
void halos() {
  Design d;
  d.layout = {0, 0, 80, 40, 1, 1};
  d.num_movable = 4;
  d.node_x = {25, 25, 25, 25, 60, 0};
  d.node_y = {20, 20, 20, 20, 25, 20};
  d.node_size_x = {10, 8.3, 2, 2, 8, 0};
  d.node_size_y = {6, 5.7, 1, 1, 5, 0};
  d.node_halo_x = {3, 1.2, 0, 0, 2, 0};
  d.node_halo_y = {2, 3, 0, 0, 1, 0};
  d.pin2node = {0, 1, 1, 2, 2, 3, 3, 4, 0, 5};
  d.pin_offset_x.assign(10, .5);
  d.pin_offset_y.assign(10, .25);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  d.flat_net2pin_start = {0, 2, 4, 6, 8, 10};
  d.build_maps();
  const auto original = d;
  Options o;
  o.bins_x = o.bins_y = 8;
  o.max_iterations = 80;
  int callbacks = 0;
  o.observer = [&](const Design &physical, PlacementStage, int) {
    check(physical.node_size_x == original.node_size_x &&
              physical.node_size_y == original.node_size_y,
          "observer exposed padded dimensions");
    check(physical.pin_offset_x == original.pin_offset_x &&
              physical.pin_offset_y == original.pin_offset_y,
          "observer exposed padded pin offsets");
    ++callbacks;
  };
  auto result = place(d, o);
  check(result.legal && is_legal(d) && callbacks > 0, "halo placement failed");
  check(d.node_size_x == original.node_size_x &&
            d.node_size_y == original.node_size_y &&
            d.pin_offset_x == original.pin_offset_x &&
            d.pin_offset_y == original.pin_offset_y,
        "halo changed physical geometry");
  near(hpwl(d), result.hpwl, 1e-9, "halo HPWL used occupancy pins");
  check(d.node_x[4] == original.node_x[4] && d.node_y[4] == original.node_y[4],
        "halo moved fixed macro");
  for (int i = 0; i < 4; ++i) {
    const auto hx = std::ceil(d.node_halo_x[i]),
               hy = std::ceil(d.node_halo_y[i]);
    check(d.node_x[i] >= hx && d.node_y[i] >= hy &&
              d.node_x[i] + d.node_size_x[i] + hx <= 80 &&
              d.node_y[i] + d.node_size_y[i] + hy <= 40,
          "halo outside core");
    for (int j = i + 1; j < 5; ++j) {
      const double x = hx + std::ceil(d.node_halo_x[j]),
                   y = hy + std::ceil(d.node_halo_y[j]);
      check(d.node_x[i] + d.node_size_x[i] + x <= d.node_x[j] + 1e-6 ||
                d.node_x[j] + d.node_size_x[j] + x <= d.node_x[i] + 1e-6 ||
                d.node_y[i] + d.node_size_y[i] + y <= d.node_y[j] + 1e-6 ||
                d.node_y[j] + d.node_size_y[j] + y <= d.node_y[i] + 1e-6,
            "missing halo clearance");
    }
  }
  auto invalid = d;
  invalid.node_x[2] = d.node_x[0] + d.node_size_x[0];
  invalid.node_y[2] = d.node_y[0];
  check(!is_legal(invalid), "legal checker ignored halo-only overlap");
  invalid = d;
  invalid.node_halo_x.pop_back();
  rejects([&] { invalid.validate(); });
  invalid = d;
  invalid.node_halo_y[0] = -1;
  rejects([&] { invalid.validate(); });
  invalid = d;
  invalid.node_halo_x[5] = 1;
  rejects([&] { invalid.validate(); });
  invalid = d;
  invalid.node_halo_x[0] = 100;
  const auto before = invalid;
  rejects([&] { place(invalid, o); });
  check(invalid.node_x == before.node_x && invalid.node_y == before.node_y,
        "failed halo placement mutated input");
}

void transparent_macros() {
  MacroTransparency schedule{5, 100, 20};
  for (int i = 0; i < 1000; ++i) {
    const bool expected = i >= 100 && i < 300 && (i / 20) % 2 == 1;
    check(macro_transparent(schedule, i) == expected,
          "wrong transparency phase boundary");
    check(!macro_transparent({}, i), "default transparency must be disabled");
  }
  Options o;
  o.bins_x = o.bins_y = 8;
  o.num_threads = 1;
  o.min_iterations = o.max_iterations = 120;
  o.legalize = o.detailed_place = false;
  o.add_fillers = true;
  o.filler_max_count = 64;
  o.target_density = .5;
  o.fixed_macro_overlap_weight = 1;
  auto base = design();
  base.node_macro.assign(base.node_x.size(), 0);
  base.node_macro[6] = 1;
  auto plain = base, cycle = base;
  const auto baseline = place(plain, o);
  o.macro_transparency = {2, 20, 20};
  const auto result = place(cycle, o);
  check(result.iterations == 120,
        "transparency stopped before opaque settling");
  check(result.filler_count == baseline.filler_count &&
            result.filler_area == baseline.filler_area,
        "transparency changed the filler population");
  check(result.history.size() == 6, "missing transparency phase samples");
  for (const auto &s : result.history) {
    const bool transparent = s.iteration == 40 || s.iteration == 80;
    near(s.macro_opacity, transparent ? 0 : 1, 0,
         "wrong reported macro opacity");
    const double effective_fixed = transparent ? 0 : 8;
    const double expected_scale =
        (.5 * (400 - effective_fixed) - 12) / result.filler_area;
    near(s.filler_charge_scale, expected_scale, 1e-12,
         "inconsistent transparent filler charge");
    if (s.iteration == 20) {
      near(s.hpwl, baseline.history[0].hpwl, 0,
           "transparency changed warmup HPWL");
      near(s.overflow, baseline.history[0].overflow, 0,
           "transparency changed warmup density");
    }
  }
  check(plain.node_x != cycle.node_x || plain.node_y != cycle.node_y,
        "transparency did not affect movement");
  for (int i = base.num_movable; i < int(base.node_x.size()); ++i) {
    near(cycle.node_x[i], base.node_x[i], 0, "transparency moved a fixed node");
    near(cycle.node_y[i], base.node_y[i], 0, "transparency moved a fixed node");
  }
  check(cycle.node_size_x == base.node_size_x &&
            cycle.pin_offset_x == base.pin_offset_x &&
            cycle.flat_net2pin == base.flat_net2pin,
        "transparency changed physical geometry/connectivity");
  // The final legalizer receives fully opaque original physical geometry.
  auto legal = base;
  o.legalize = true;
  check(place(legal, o).legal && is_legal(legal),
        "transparent-phase overlap escaped legalization");
  o.legalize = false;
  auto invalid = base;
  o.max_iterations = 100;
  rejects([&] { place(invalid, o); });
  check(invalid.node_x == base.node_x,
        "invalid transparency schedule mutated input");
  o.max_iterations = 120;
  GpState state;
  rejects([&] { place(invalid, o, &state); });
  invalid.node_macro[0] = 1;
  rejects([&] { place(invalid, o); });

  // Permanent transparency must behave like zero-area fixed terminals with
  // identical pin coordinates, including density feedback and the GP forces.
  o.macro_transparency = {};
  o.macro_transparency.always = true;
  o.add_fillers = false;
  auto always = base;
  std::vector<double> physical_overflows;
  o.observer_interval = 20;
  o.observer = [&](const Design &snapshot, PlacementStage stage, int) {
    if (stage == PlacementStage::global_placement) {
      physical_overflows.push_back(density_overflow(snapshot, 8, 8, .5));
    }
  };
  const auto all = place(always, o);
  for (const auto &sample : all.history) {
    near(sample.macro_opacity, 0, 0, "permanent transparency became opaque");
    near(sample.raw_overflow, physical_overflows[sample.iteration / 20 - 1],
         1e-12, "raw overflow lost physical macro occupancy");
  }
  auto zero_area = base;
  zero_area.node_size_x[6] = zero_area.node_size_y[6] = 0;
  o.macro_transparency.always = false;
  o.fixed_macro_overlap_weight = 0;
  o.observer = {};
  const auto reference = place(zero_area, o);
  check(always.node_x == zero_area.node_x && always.node_y == zero_area.node_y,
        "transparent macro differs from fixed pins without a body");
  near(all.overflow, reference.overflow, 0,
       "transparent overflow still subtracts macro area");
  check(always.node_size_x == base.node_size_x &&
            always.node_size_y == base.node_size_y &&
            always.pin_offset_x == base.pin_offset_x &&
            always.pin_offset_y == base.pin_offset_y &&
            always.flat_net2pin == base.flat_net2pin,
        "permanent transparency changed macro geometry or pins");
  for (int i = base.num_movable; i < int(base.node_x.size()); ++i) {
    near(always.node_x[i], base.node_x[i], 0,
         "permanent transparency moved fixed x");
    near(always.node_y[i], base.node_y[i], 0,
         "permanent transparency moved fixed y");
  }
  o.macro_transparency.always = true;
  o.add_fillers = true;
  auto filled = base;
  const auto fill = place(filled, o);
  near(fill.history.front().filler_charge_scale,
       (.5 * 400 - 12) / fill.filler_area, 1e-12,
       "permanent transparency did not compensate fillers from step one");
  o.legalize = true;
  auto final_legal = base;
  check(place(final_legal, o).legal && is_legal(final_legal),
        "permanent transparent GP changed final physical legality");
  o.macro_transparency.cycles = 1;
  rejects([&] { place(invalid, o); });
}

void vanished_memories() {
  near(vanished_memory_scale(1), .6, 0, "initial memory size");
  near(vanished_memory_scale(70), .6, 0, "hold initial size through 70");
  near(vanished_memory_scale(71), .616, 1e-12, "memory growth starts at 71");
  near(vanished_memory_scale(75), .68, 1e-12, "memory width at 75");
  near(vanished_memory_scale(94), .984, 1e-12, "last partial memory step");
  near(vanished_memory_scale(95), 1, 0, "memory full at 95");
  near(vanished_memory_scale(1000), 1, 0, "memory stays full");
  near(vanished_memory_scale(100, 0, 100, 100), 0, 0,
       "explicit vanished start");
  near(vanished_memory_scale(101, 0, 100, 100), .01, 0,
       "explicit original ramp");
  near(vanished_memory_scale(100, .8, 100, 100), .8, 0,
       "partial initial memory");
  near(vanished_memory_scale(150, .8, 100, 100), .9, 1e-12,
       "standard partial ramp");
  near(vanished_memory_scale(200, .5, 200, 100), .75, 0, "slow partial ramp");
  near(vanished_memory_scale(300, .5, 200, 100), 1, 0, "slow ramp ends at 300");
  near(vanished_memory_local_weight(.67, 70), .67, 0, "local holds through 70");
  near(vanished_memory_local_weight(.67, 71), .62, 1e-12,
       "local starts falling at 71");
  near(vanished_memory_local_weight(.67, 83), .02, 1e-12,
       "local ramp before zero");
  near(vanished_memory_local_weight(.67, 84), 0, 0, "local reaches zero at 84");
  near(vanished_memory_local_weight(.67, 81, .1), .12, 1e-12,
       "local ramp before retained floor");
  near(vanished_memory_local_weight(.67, 82, .1), .1, 0,
       "local reaches retained floor at 82");
  near(vanished_memory_local_weight(.67, 1000, .1), .1, 0,
       "local floor persists after memory growth");
  auto base = design();
  base.node_macro.assign(base.node_x.size(), 0);
  base.node_macro[6] = 1;
  Options o;
  o.vanished_memory = true;
  o.min_iterations = o.max_iterations = 340;
  o.bins_x = o.bins_y = 8;
  o.local_density_bins = 16;
  o.local_density_weight = .67;
  o.local_density_target = .05;
  o.local_density_jitter = 0;
  o.num_threads = 1;
  o.target_density = .5;
  o.add_fillers = true;
  o.filler_max_count = 64;
  o.fixed_macro_overlap_weight = 1;
  o.legalize = true;
  o.detailed_place = false;
  std::vector<double> expected;
  o.observer_interval = 20;
  o.observer = [&](const Design &d, PlacementStage stage, int it) {
    check(d.node_size_x == base.node_size_x &&
              d.node_size_y == base.node_size_y &&
              d.pin_offset_x == base.pin_offset_x &&
              d.pin_offset_y == base.pin_offset_y,
          "memory growth leaked into physical geometry or pins");
    for (int i = base.num_movable; i < int(base.node_x.size()); ++i) {
      near(d.node_x[i], base.node_x[i], 0, "memory growth moved fixed x");
      near(d.node_y[i], base.node_y[i], 0, "memory growth moved fixed y");
    }
    if (stage == PlacementStage::global_placement) {
      auto body = d;
      const double scale = vanished_memory_scale(
          it, o.vanished_memory_initial_scale,
          o.vanished_memory_growth_iterations, o.vanished_memory_start);
      body.node_x[6] += body.node_size_x[6] * (1 - scale) / 2;
      body.node_y[6] += body.node_size_y[6] * (1 - scale) / 2;
      body.node_size_x[6] *= scale;
      body.node_size_y[6] *= scale;
      expected.push_back(density_overflow(body, 8, 8, .5));
    }
  };
  for (auto [initial, duration] : {std::pair{.6, 25},
                                   {0., 100},
                                   {.5, 100},
                                   {.8, 100},
                                   {.5, 200},
                                   {.8, 200}}) {
    o.vanished_memory_initial_scale = initial;
    o.vanished_memory_growth_iterations = duration;
    for (double floor : {0., .1}) {
      o.vanished_memory_local_floor = floor;
      expected.clear();
      auto d = base;
      const auto result = place(d, o);
      check(result.legal && is_legal(d), "vanished-memory final legality");
      check(result.iterations == 340 && result.history.size() == 17,
            "schedule ended early");
      for (size_t i = 0; i < result.history.size(); ++i) {
        const auto &h = result.history[i];
        const double scale =
            vanished_memory_scale(h.iteration, initial, duration);
        near(h.macro_exclusion_scale, scale, 0, "wrong memory extent report");
        near(h.local_density_weight,
             vanished_memory_local_weight(.67, h.iteration, floor), 0,
             "wrong local report");
        near(h.overflow, expected[i], 1e-12,
             "density still uses full memory geometry");
        near(h.filler_charge_scale,
             (.5 * (400 - 8 * scale * scale) - 12) / result.filler_area, 1e-12,
             "filler charge disagrees with grown memory area");
      }
    }
  }
  o.vanished_memory_initial_scale = .6;
  o.vanished_memory_growth_iterations = 25;
  o.vanished_memory_local_floor = 0;
  o.min_iterations = o.max_iterations = 400;
  for (int start : {50, 70, 100, 150, 200, 250, 300}) {
    o.vanished_memory_start = start;
    near(vanished_memory_scale(start, .6, 25, start), .6, 0,
         "hold shifted initial size");
    near(vanished_memory_scale(start + 5, .6, 25, start), .68, 1e-12,
         "shifted growth");
    near(vanished_memory_scale(start + 25, .6, 25, start), 1, 0,
         "shifted full size");
    near(vanished_memory_local_weight(.67, start, 0, start), .67, 0,
         "hold shifted local");
    near(vanished_memory_local_weight(.67, start + 1, 0, start), .62, 1e-12,
         "shifted local decay");
    near(vanished_memory_local_weight(.67, start + 14, 0, start), 0, 0,
         "shifted local zero");
    auto d = base;
    expected.clear();
    const auto result = place(d, o);
    check(result.legal && is_legal(d) && result.iterations == 400,
          "shifted ramp legality/iterations");
    check(result.history.size() == 20, "shifted ramp history");
    for (size_t i = 0; i < result.history.size(); ++i) {
      const auto &h = result.history[i];
      near(h.macro_exclusion_scale,
           vanished_memory_scale(h.iteration, .6, 25, start), 0,
           "shifted extent report");
      near(h.local_density_weight,
           vanished_memory_local_weight(.67, h.iteration, 0, start), 0,
           "shifted local report");
      near(h.overflow, expected[i], 1e-12, "shifted density footprint");
    }
  }
  o.observer = {};
  auto invalid = base;
  o.vanished_memory_start = 300;
  o.vanished_memory_growth_iterations = 1;
  o.gp_search.enabled = true;
  o.gp_search.start = 321;
  // Full memory at 301 is not enough: local pressure settles only at 334.
  rejects([&] { place(invalid, o); });
  o.gp_search.enabled = false;
  o.vanished_memory_growth_iterations = 25;
  o.vanished_memory_start = -1;
  rejects([&] { place(invalid, o); });
  o.vanished_memory_start = 1000001;
  rejects([&] { place(invalid, o); });
  o.vanished_memory_start = 300;
  o.max_iterations = 340;
  rejects(
      [&] { place(invalid, o); }); // Requires at least 345 including settling.
  o.vanished_memory_start = 100;
  o.vanished_memory_growth_iterations = 200;
  o.min_iterations = 340;
  o.max_iterations = 300;
  rejects([&] {
    place(invalid, o);
  }); // 200-step growth needs opaque settling past 300
  o.max_iterations = 340;
  o.vanished_memory_initial_scale = 1.1;
  rejects([&] { place(invalid, o); });
  o.vanished_memory_initial_scale = 0;
  o.vanished_memory_growth_iterations = 0;
  rejects([&] { place(invalid, o); });
  o.vanished_memory_growth_iterations = 100;
  o.macro_transparency.always = true;
  auto bad = base;
  rejects([&] { place(bad, o); });
  o.macro_transparency.always = false;
  o.max_iterations = 200;
  rejects([&] { place(bad, o); });
  o.max_iterations = 340;
  o.vanished_memory_local_floor = -.1;
  rejects([&] { place(bad, o); });
  o.vanished_memory_local_floor = .68;
  rejects([&] { place(bad, o); });
  o.vanished_memory_local_floor = .1;
  o.vanished_memory = false;
  rejects([&] { place(bad, o); });
  check(bad.node_x == base.node_x, "invalid growth schedule changed input");
}

void moving_cells() {
  Options o;
  o.max_iterations = o.min_iterations = 4;
  o.bins_x = o.bins_y = 8;
  o.legalize = o.detailed_place = false;
  o.num_threads = 1;
  auto plain = design(), off = plain;
  place(plain, o);
  o.movable_update = [](const Design &, int) -> std::vector<MovablePosition> {
    throw std::runtime_error("disabled movable callback was called");
  };
  place(off, o);
  check(plain.node_x == off.node_x && plain.node_y == off.node_y,
        "disabled movable update changed GP");
  o.movable_interval = 4;
  for (bool haloed : {false, true}) {
    auto d = design();
    d.node_macro.assign(d.node_x.size(), 0);
    d.node_macro[6] = 1;
    if (haloed) {
      d.node_halo_x.assign(d.node_x.size(), 0);
      d.node_halo_y.assign(d.node_x.size(), 0);
      d.node_halo_x[0] = d.node_halo_y[0] = 1;
    }
    const auto before = d;
    int calls = 0;
    o.add_fillers = true;
    o.target_density = .8;
    o.filler_max_count = 64;
    o.movable_update = [&](const Design &snapshot, int iteration) {
      ++calls;
      check(iteration == 4 && snapshot.node_x.size() == before.node_x.size(),
            "movable callback iteration or filler isolation");
      check(snapshot.node_size_x == before.node_size_x &&
                snapshot.pin_offset_x == before.pin_offset_x,
            "movable callback must see physical halo geometry");
      return std::vector<MovablePosition>{{0, 12, 22}};
    };
    place(d, o);
    check(calls == 1, "missing movable callback");
    near(d.node_x[0], 12, 0, "jumped x overwritten by stale iterate");
    near(d.node_y[0], 22, 0, "jumped y overwritten by stale iterate");
    for (int i = d.num_movable; i < int(d.node_x.size()); ++i) {
      near(d.node_x[i], before.node_x[i], 0, "jump moved fixed x");
      near(d.node_y[i], before.node_y[i], 0, "jump moved fixed y");
    }
    for (auto invalid : std::vector<std::vector<MovablePosition>>{
             {{6, 12, 22}},
             {{0, 12, 22}, {0, 13, 22}},
             {{-1, 12, 22}},
             {{0, -100, 22}},
             {{0, std::numeric_limits<double>::quiet_NaN(), 22}}}) {
      auto candidate = before;
      o.movable_update = [invalid](const Design &, int) { return invalid; };
      rejects([&] { place(candidate, o); });
      check(candidate.node_x == before.node_x &&
                candidate.node_y == before.node_y,
            "invalid jump mutated caller design");
    }
  }
  auto d = design();
  o.movable_interval = 1;
  o.observer_interval = 1;
  int held_samples = 0, released_samples = 0;
  o.movable_update = [](const Design &, int iteration) {
    return iteration == 1 ? std::vector<MovablePosition>{{0, 12, 22, 2}}
                          : std::vector<MovablePosition>{};
  };
  o.observer = [&](const Design &snapshot, PlacementStage stage,
                   int iteration) {
    if (stage != PlacementStage::global_placement)
      return;
    if (iteration == 2 || iteration == 3) {
      near(snapshot.node_x[0], 12, 0, "hold expired too early in x");
      near(snapshot.node_y[0], 22, 0, "hold expired too early in y");
      ++held_samples;
    }
    if (iteration == 4) {
      check(snapshot.node_x[0] != 12 || snapshot.node_y[0] != 22,
            "cell not released after the requested two iterations");
      ++released_samples;
    }
  };
  place(d, o);
  check(held_samples == 2 && released_samples == 1,
        "missing hold/release observations");
  GpState state;
  rejects([&] { place(d, o, &state); });
}

void moving_terminals() {
  auto base = design();
  Options o;
  o.max_iterations = 12;
  o.bins_x = o.bins_y = 8;
  o.add_fillers = false;
  o.legalize = o.detailed_place = false;
  o.num_threads = 1;
  auto plain = base, off = base;
  place(plain, o);
  o.terminal_update = [&](const Design &,
                          int) -> std::vector<TerminalPosition> {
    throw std::runtime_error("disabled terminal callback was called");
  };
  place(off, o);
  check(off.node_x == plain.node_x && off.node_y == plain.node_y,
        "disabled terminal adjustment changes GP");
  for (bool haloed : {false, true}) {
    auto d = base;
    d.node_macro.assign(d.node_x.size(), 0);
    if (haloed) {
      d.node_macro[6] = 1;
      d.node_halo_x.assign(d.node_x.size(), 0);
      d.node_halo_y.assign(d.node_x.size(), 0);
      d.node_halo_x[6] = d.node_halo_y[6] = 1;
    }
    int callbacks = 0, observed = 0;
    double latest_x = d.node_x[7];
    o.terminal_interval = 2;
    o.terminal_update = [&](const Design &snapshot, int iteration) {
      check(snapshot.node_size_x == base.node_size_x &&
                snapshot.pin_offset_x == base.pin_offset_x,
            "terminal callback must use physical geometry with halo");
      near(snapshot.node_x[7], latest_x, 0,
           "terminal movement lost between iterations");
      latest_x = iteration % 4 == 0 ? 11 : 12;
      ++callbacks;
      return std::vector<TerminalPosition>{{7, latest_x, 25}};
    };
    o.observer = [&](const Design &snapshot, PlacementStage stage, int) {
      if (stage == PlacementStage::finished) {
        near(snapshot.node_x[7], latest_x, 0,
             "observer used stale IO positions");
        ++observed;
      }
    };
    o.checkpoint_interval = 2;
    o.checkpoint = [&](const Design &snapshot, int) {
      near(snapshot.node_x[7], latest_x, 0,
           "STA checkpoint used stale IO positions");
      return GpUpdate{};
    };
    const auto result = place(d, o);
    o.checkpoint = {};
    check(callbacks > 0 && observed == 1, "terminal callback/observer missing");
    near(d.node_x[7], latest_x, 0, "moved IO not written back");
    near(hpwl(d), result.hpwl, 1e-9, "reported HPWL uses stale IO positions");
    near(d.node_x[6], base.node_x[6], 0, "moving terminals moved fixed macro");
    o.observer = {};
    o.terminal_update = [&](const Design &, int) {
      return std::vector<TerminalPosition>{{6, 11, 25}};
    };
    auto invalid = base;
    rejects([&] { place(invalid, o); });
    check(invalid.node_x == base.node_x,
          "invalid terminal update mutated input");
  }
  o.terminal_update = [](const Design &, int) {
    return std::vector<TerminalPosition>{};
  };
  o.add_fillers = true;
  o.filler_max_count = 64;
  const auto physical_count = base.node_x.size();
  o.terminal_update = [&](const Design &snapshot, int) {
    snapshot.validate();
    check(snapshot.node_x.size() == physical_count,
          "virtual fillers escaped to IO callback");
    return std::vector<TerminalPosition>{};
  };
  const auto filled = place(base, o);
  check(filled.filler_count > 0 && filled.filler_count <= 64 &&
            base.node_x.size() == physical_count,
        "virtual filler budget or export");
  o.add_fillers = false;
  GpState state;
  rejects([&] { place(base, o, &state); });
}

// Fixed nodes only need finite coordinates; one far above the core used to turn
// an out-of-range double into a negative row index inside is_legal().
void far_fixed_nodes() {
  for (const double far : {1e12, -1e12, 3e9}) {
    for (const auto &[far_x, far_y] :
         {std::pair{far, far}, std::pair{far, 5.}, std::pair{5., far}}) {
      Design d;
      d.layout = {0, 0, 40, 20, 1, 1};
      d.num_movable = 4;
      d.node_x = {2, 5, 8, 20, far_x, 30};
      d.node_y = {2, 6, 2, 9, far_y, 5};
      d.node_size_x = {2, 3, 2, 4, 6, 0};
      d.node_size_y = {1, 1, 1, 1, 3, 0};
      d.pin2node = {0, 1, 2, 3, 4, 5};
      d.pin_offset_x.assign(6, 0.5);
      d.pin_offset_y.assign(6, 0.5);
      d.flat_net2pin = {0, 1, 2, 3, 4, 5};
      d.flat_net2pin_start = {0, 3, 6};
      d.build_maps();
      const auto before = d;
      (void)is_legal(d); // must not index outside the row table
      Options o;
      o.max_iterations = 30;
      o.bins_x = o.bins_y = 8;
      const auto r = place(d, o);
      check(r.legal && is_legal(d), "far fixed node broke legalization");
      check(d.node_x[4] == before.node_x[4] && d.node_y[4] == before.node_y[4],
            "far fixed node moved");
    }
  }
}
void power2_grid() {
  // ISM sizes its bin grid with ceil_power2(bins / 8); fewer than eight bins
  // gives zero, whose log2 used to be cast to int and used as a shift count.
  using DreamPlace::ceil_power2;
  const int expected[] = {1, 1, 2, 4, 4, 8, 8, 8, 8, 16};
  for (int v = 0; v < 10; ++v)
    check(ceil_power2(v) == expected[v], "ceil_power2 mismatch");
  check(ceil_power2(-3) == 1, "ceil_power2 of a negative count");
}
void timing_spreading() {
  // Symmetric crowded cells with coincident fixed terminals: initial WL force
  // is zero even with different net weights. Movement tests the spreading
  // policy itself, rather than the existing timing-weighted WL attraction.
  Design seed;
  seed.layout = {0, 0, 64, 64, 1, 1};
  seed.num_movable = 2;
  seed.node_x = {29.5, 33.5, 30, 34};
  seed.node_y = {31.5, 31.5, 32, 32};
  seed.node_size_x = {1, 1, 0, 0};
  seed.node_size_y = {1, 1, 0, 0};
  seed.node_macro = {0, 0, 0, 0};
  seed.node_inflatable = {1, 1, 0, 0};
  seed.pin2node = {0, 2, 1, 3};
  seed.pin_offset_x = seed.pin_offset_y = {.5, 0, .5, 0};
  seed.flat_net2pin = {0, 1, 2, 3};
  seed.flat_net2pin_start = {0, 2, 4};
  seed.net_weights = {5, 1};
  seed.build_maps();
  Options o;
  o.bins_x = o.bins_y = 8;
  o.max_iterations = 1;
  o.add_fillers = o.legalize = o.detailed_place = false;
  o.local_density_jitter = 0;
  o.local_density_bins = 64;
  o.target_density = .01;
  auto plain = seed, protected_cells = seed;
  place(plain, o);
  o.timing_spread_weight = 1;
  place(protected_cells, o);
  check(plain.node_x == protected_cells.node_x &&
            plain.node_y == protected_cells.node_y,
        "timing spreading changed ordinary physical-density force");
  o.local_density_weight = 20;
  protected_cells = seed;
  const auto protected_result = place(protected_cells, o);
  const double critical_move = seed.node_x[0] - protected_cells.node_x[0];
  const double slack_move = protected_cells.node_x[1] - seed.node_x[1];
  check(critical_move > 0 && critical_move < .5 * slack_move,
        "critical cell must receive less local spreading than its slack "
        "neighbor");
  check(protected_result.timing_spread_cells == 1 &&
            protected_result.timing_spread_min_mobility == .2,
        "timing spreading census");
  check(protected_cells.node_x[2] == seed.node_x[2] &&
            protected_cells.node_y[3] == seed.node_y[3] &&
            protected_cells.node_size_x == seed.node_size_x,
        "timing spreading changed fixed terminals or physical sizes");
  auto uniform = seed;
  uniform.net_weights = {1, 1};
  auto uniform_on = uniform;
  place(uniform_on, o);
  o.timing_spread_weight = 0;
  place(uniform, o);
  check(uniform.node_x == uniform_on.node_x &&
            uniform.node_y == uniform_on.node_y,
        "unit weights must leave spreading exactly unchanged");
  auto masked = seed;
  masked.net_mask[0] = 0;
  o.timing_spread_weight = 1;
  check(place(masked, o).timing_spread_cells == 0,
        "masked nets must not protect cells from spreading");
  // A timing checkpoint changes which cell is critical. The cached per-cell
  // mobility must follow the new weights, even without topology replacement.
  auto refreshed = seed;
  o.max_iterations = 2;
  o.checkpoint_interval = 1;
  o.checkpoint = [](const Design &, int) {
    GpUpdate update;
    update.net_weights = {1, 2};
    return update;
  };
  const auto refreshed_result = place(refreshed, o);
  check(refreshed_result.timing_spread_cells == 1 &&
            refreshed_result.timing_spread_min_mobility == .5,
        "checkpoint did not refresh spreading priority");
  o.checkpoint = {};
  // A fixed-terminal crossing net produces the same demand around both
  // movable cells without adding WL force to them. Protection must retain
  // outward movement while allowing the slack cell to move farther.
  o.local_density_weight = 0;
  o.max_iterations = 1;
  o.routability = true;
  o.route_capacity_h = o.route_capacity_v = 1e-6;
  auto routed = seed;
  routed.node_x.insert(routed.node_x.end(), {20, 44});
  routed.node_y.insert(routed.node_y.end(), {28, 36});
  routed.node_size_x.insert(routed.node_size_x.end(), {0, 0});
  routed.node_size_y.insert(routed.node_size_y.end(), {0, 0});
  routed.node_macro.insert(routed.node_macro.end(), {0, 0});
  routed.node_inflatable.insert(routed.node_inflatable.end(), {0, 0});
  routed.pin2node.insert(routed.pin2node.end(), {4, 5});
  routed.pin_offset_x.insert(routed.pin_offset_x.end(), {0, 0});
  routed.pin_offset_y.insert(routed.pin_offset_y.end(), {0, 0});
  routed.flat_net2pin.insert(routed.flat_net2pin.end(), {4, 5});
  routed.flat_net2pin_start.push_back(6);
  routed.net_weights.push_back(1);
  routed.net_mask.push_back(1);
  routed.build_maps();
  GpState state;
  place(routed, o, &state);
  check(state.inflation[0] > 1 && state.inflation[1] > 1,
        "routing spreading fixture did not create congestion");
  check(seed.node_x[0] - routed.node_x[0] > 0 &&
            seed.node_x[0] - routed.node_x[0] <
                routed.node_x[1] - seed.node_x[1],
        "routing receiver protection must favor the slack cell");
  auto changed = o;
  changed.timing_spread_weight = .5;
  rejects([&] { place(routed, changed, &state); });
  for (double value : {-1., 1.01, std::numeric_limits<double>::quiet_NaN()}) {
    o.timing_spread_weight = value;
    auto bad = seed;
    rejects([&] { place(bad, o); });
  }
}
} // namespace
void abacus_weights() {
  // Null weights are unit weights, and large opt-in weights must not reach the
  // fixed-cluster sentinel M (that aborted inside the row solver).
  const std::vector<double> init_x = {0.2, 0.4, 0.6, 0.8};
  const std::vector<double> size_x = {1, 1, 1, 1}, size_y = {1, 1, 1, 1};
  const auto row = [&](const double *weights) {
    std::vector<double> x = {0, 1, 2, 3};
    std::vector<int> nodes = {0, 1, 2, 3};
    std::vector<DreamPlace::AbacusCluster<double>> clusters(4);
    check(DreamPlace::abacusPlaceRowCPU(
              init_x.data(), size_x.data(), size_y.data(), weights, x.data(),
              1.0, 1.0, 0.0, 10.0, 4, 4, nodes.data(), clusters.data(), 4),
          "abacus row failed");
    return x;
  };
  const std::vector<double> ones(4, 1), large(4, 1e6);
  const auto unit = row(nullptr);
  check(unit == row(ones.data()), "null abacus weights differ from unit");
  const auto heavy = row(large.data());
  for (int i = 0; i < 4; ++i)
    near(heavy[i], unit[i], 1e-9, "uniform abacus weights moved the optimum");
}
void after_macro_legalize_writeback() {
  auto d = design();
  Options o;
  o.bins_x = o.bins_y = 8;
  o.max_iterations = 40;
  int calls = 0;
  o.after_macro_legalize = [&](Design &snapshot) {
    snapshot.node_legalize_weights.assign(snapshot.node_x.size(), 2);
    ++calls;
  };
  check(place(d, o).legal && calls == 1, "after_macro_legalize placement");
  check(d.node_legalize_weights == std::vector<double>(d.node_x.size(), 2),
        "after_macro_legalize edits were dropped");
}
void macro_only_legalization() {
  for (bool halos : {false, true}) {
    Design d;
    d.layout = {0, 0, 120, 80, 1, 10};
    d.num_movable = 4;
    d.node_x = {40, 40, 40.25, 40.25, 90};
    d.node_y = {20, 20, 20.5, 20.5, 0};
    d.node_size_x = {20, 20, 2, 2, 10};
    d.node_size_y = {20, 20, 10, 30, 80};
    d.node_macro = {1, 1, 0, 0, 1};
    if (halos) {
      d.node_halo_x = {2, 2, 0, 0, 0};
      d.node_halo_y = {10, 10, 0, 0, 0};
    }
    d.build_maps();
    const auto original = d;
    Options o;
    o.max_iterations = 0;
    o.bins_x = o.bins_y = 4;
    o.add_fillers = o.detailed_place = false;
    o.macro_snap = halos;
    o.macro_only_legalize = true;
    bool before = false, after = false, observed = false;
    o.before_legalize = [&](Design &snapshot) {
      check(!macros_are_legal(snapshot), "macro fixture must overlap");
      before = true;
    };
    o.after_macro_legalize = [&](Design &) { after = true; };
    o.observer = [&](const Design &, PlacementStage stage, int) {
      check(stage != PlacementStage::detailed_placement,
            "macro-only mode performed cell detail");
      observed |= stage == PlacementStage::legalization;
    };
    const auto result = place(d, o);
    check(before && after && observed && macros_are_legal(d),
          "macro-only legalization must retain macro checks and callbacks");
    check(!result.legal && !is_legal(d),
          "provisional cells must not be reported fully legal");
    check(result.alignment.cell_legalization_seconds == 0 &&
              result.detailed_history.empty(),
          "macro-only mode ran cell legalization or detail");
    for (size_t i = 2; i < d.node_x.size(); ++i)
      check(d.node_x[i] == original.node_x[i] &&
                d.node_y[i] == original.node_y[i],
            "macro-only legalization moved a cell or fixed obstacle");
    o.macro_only_legalize = false;
    o.before_legalize = {};
    o.after_macro_legalize = {};
    o.observer = {};
    check(place(d, o).legal && is_legal(d),
          "full legalization must subsequently legalize the cells");
    o.macro_only_legalize = o.detailed_place = true;
    rejects([&] { place(d, o); });
    o.legalize = o.detailed_place = false;
    rejects([&] { place(d, o); });
  }
}
void virtual_macro_pulls() {
  auto d = design();
  d.node_macro.assign(d.node_x.size(), 0);
  d.node_macro[0] = 1;
  const auto before = d;
  Options o;
  o.max_iterations = 30;
  o.min_iterations = 100;
  o.legalize = o.detailed_place = o.add_fillers = false;
  o.bins_x = o.bins_y = 8;
  o.macro_pull_interval = 1;
  int calls = 0;
  o.macro_pull_update =
      [&](const Design &snapshot,
          int iteration) -> std::optional<std::vector<MacroPull>> {
    check(snapshot.pin2node == before.pin2node &&
              snapshot.net_mask == before.net_mask,
          "virtual pulls changed physical connectivity");
    ++calls;
    if (iteration == 5 || iteration == 15)
      return std::vector<MacroPull>{{0, iteration == 5 ? 45. : 15., 25., .5}};
    if (iteration == 25)
      return std::vector<MacroPull>{};
    return std::nullopt;
  };
  place(d, o);
  check(calls == 30,
        "virtual pulls require independent safe iteration boundaries");
  check(d.node_x.size() == before.node_x.size() &&
            d.pin2node == before.pin2node &&
            d.flat_net2pin == before.flat_net2pin &&
            d.net_weights == before.net_weights,
        "virtual pull leaked nodes or edges");
  for (size_t i = d.num_movable; i < d.node_x.size(); ++i) {
    near(d.node_x[i], before.node_x[i], 0, "virtual pull moved fixed node");
  }
  auto left = before, right = before;
  o.macro_pull_update =
      [](const Design &,
         int iteration) -> std::optional<std::vector<MacroPull>> {
    if (iteration == 1)
      return std::vector<MacroPull>{{0, 12, 25, 1}};
    return std::nullopt;
  };
  place(left, o);
  o.macro_pull_update =
      [](const Design &,
         int iteration) -> std::optional<std::vector<MacroPull>> {
    if (iteration == 1)
      return std::vector<MacroPull>{{0, 48, 25, 1}};
    return std::nullopt;
  };
  place(right, o);
  check(right.node_x[0] > left.node_x[0],
        "virtual pull failed to attract macro toward endpoint");
  // Wire-only best-feasible rollback must not discard a soft macro objective.
  auto guided = o;
  guided.max_iterations = 80;
  guided.min_iterations = 100;
  guided.stop_overflow = 1;
  guided.density_upper_bound = true;
  auto guided_design = before;
  // Wire-indifferent nets make the first feasible snapshot look equally good
  // even while a separate macro pull is still settling.
  guided_design.net_weights.assign(guided_design.flat_net2pin_start.size() - 1,
                                   0);
  const auto guided_result = place(guided_design, guided);
  check(guided_result.selected_iteration == guided_result.iterations,
        "wire-only rollback discarded a macro-guided placement");
  // A deliberately overfilled layout never meets the zero-overflow gate.
  auto dense = before;
  for (int i = 0; i < dense.num_movable; ++i) {
    dense.node_size_x[i] = 20;
    dense.node_size_y[i] = 5;
  }
  calls = 0;
  o.macro_pull_after_density = true;
  o.stop_overflow = 0;
  o.target_density = .01;
  o.macro_pull_update = [&](const Design &,
                            int) -> std::optional<std::vector<MacroPull>> {
    ++calls;
    return std::nullopt;
  };
  place(dense, o);
  check(calls == 0, "macro pull ignored native density gate");
  o.target_density = .7;
  o.stop_overflow = 1;
  auto ready = before;
  place(ready, o);
  check(calls > 0, "density-ready macro callback did not run");
  o.macro_pull_after_density = false;
  o.max_iterations = 60;
  o.min_iterations = 1;
  o.macro_pull_settle_iterations = 20;
  o.macro_pull_update =
      [](const Design &,
         int iteration) -> std::optional<std::vector<MacroPull>> {
    if (iteration == 20)
      return std::vector<MacroPull>{{0, 48, 25, .5}};
    return std::nullopt;
  };
  auto settling = before;
  const auto settled = place(settling, o);
  check(settled.iterations >= 40,
        "GP stopped before frozen anchor settling steps");
  o.post_density_factor = 1;
  o.max_iterations = 80;
  o.macro_pull_update =
      [](const Design &,
         int iteration) -> std::optional<std::vector<MacroPull>> {
    if (iteration == 40)
      return std::vector<MacroPull>{{0, 48, 25, .5}};
    return std::nullopt;
  };
  auto deadline = before;
  const auto continued = place(deadline, o);
  check(continued.iterations >= 60,
        "continuation deadline cut short a new macro pull");
  auto bad = before;
  o.macro_pull_update = [](const Design &, int) {
    return std::optional<std::vector<MacroPull>>{{{6, 30, 25, .5}}};
  };
  rejects([&] { place(bad, o); });
  check(bad.node_x == before.node_x,
        "invalid pull partially committed placement");
}
void size_neutral_macro_pull() {
  Design d;
  d.layout = {0, 0, 10000, 10000, 1, 1};
  d.num_movable = 3;
  d.node_x = {500, 3000, 6000};
  d.node_y = {500, 3000, 6000};
  d.node_size_x = {20, 200, 20};
  d.node_size_y = {40, 400, 40};
  d.node_macro = {1, 1, 1};
  d.pin2node = {0, 1, 1, 1, 2};
  d.pin_offset_x = d.pin_offset_y = {0, 0, 0, 0, 0};
  d.flat_net2pin = {0, 1, 2, 3, 4};
  d.flat_net2pin_start = {0, 5};
  d.build_maps();
  detail::Gradient wire, density;
  wire.x = {3, 7, 5};
  wire.y = {2, 8, 4};
  density.x = {11, 1, 9};
  density.y = {4, 8, 7};
  for (const double density_weight : {0., .01, 100.}) {
    for (const double distance : {0., 10., 1000.}) {
      std::vector<MacroPull> pulls;
      for (int i = 0; i < 3; ++i)
        pulls.push_back({i, d.node_x[i] + d.node_size_x[i] / 2 + distance,
                         d.node_y[i] + d.node_size_y[i] / 2 - distance / 2,
                         .7});
      const auto forces =
          detail::macro_pull_forces(d, pulls, wire, density, density_weight);
      check(forces.size() == 3, "missing macro pull force");
      double first_x = 0, first_y = 0;
      for (int i = 0; i < 3; ++i) {
        const double pre = std::max(1., double(i == 1 ? 3 : 1) +
                                            density_weight * d.node_size_x[i] *
                                                d.node_size_y[i]);
        const double x = forces[i].x / pre, y = forces[i].y / pre;
        if (!i) {
          first_x = x;
          first_y = y;
        }
        near(x, first_x, 1e-12,
             "macro width/area/pins changed pull center response");
        near(y, first_y, 1e-12,
             "macro height/area/pins changed pull center response");
        check(distance == 0 ? x == 0 && y == 0 : x < 0 && y > 0,
              "macro pull direction or zero-offset response incorrect");
      }
    }
  }
}
void macro_sized_virtual_fillers() {
  Design design;
  design.layout = {0, 0, 1000, 1000, 1, 1};
  design.num_movable = 3;
  design.node_x = {100, 600, 400};
  design.node_y = {100, 600, 400};
  design.node_size_x = {80, 120, 2};
  design.node_size_y = {160, 120, 1};
  design.node_macro = {1, 1, 0};
  design.node_halo_x = design.node_halo_y = {8, 8, 0};
  design.build_maps();
  Options options;
  options.max_iterations = options.min_iterations = 1;
  options.bins_x = options.bins_y = 64;
  options.legalize = options.detailed_place = false;
  options.add_fillers = true;
  options.target_density = .5;
  options.filler_max_count = 4096;
  auto cell_sized = design;
  const auto small = place(cell_sized, options);
  options.filler_max_count = 0;
  options.filler_macro_fraction = .125;
  auto macro_sized = design;
  const auto large = place(macro_sized, options);
  near(large.filler_size_x, 10, 1e-12,
       "macro filler width includes halo or wrong memory");
  near(large.filler_size_y, 20, 1e-12,
       "macro filler height clipped to cell bin resolution");
  check(large.filler_count > 0 && large.filler_count < small.filler_count,
        "large fillers did not reduce particle population");
  check(std::abs(large.filler_area - small.filler_area) <=
            large.filler_size_x * large.filler_size_y +
                small.filler_size_x * small.filler_size_y,
        "macro filler sizing changed the area budget");
  check(macro_sized.node_x.size() == design.node_x.size() &&
            macro_sized.node_size_x == design.node_size_x &&
            macro_sized.node_size_y == design.node_size_y,
        "macro-sized density particles escaped or deformed physical nodes");
  options.filler_cell_space = true;
  for (const double density : {.2, .5}) {
    options.target_density = density;
    for (const double fraction : {1.5, 1., .75, .5}) {
      options.filler_area_fraction = fraction;
      auto input = design;
      const auto budgeted = place(input, options);
      const double expected =
          (1 - density) * (1000000 - 80 * 160 - 120 * 120) * fraction;
      check(budgeted.filler_area <= expected &&
                expected - budgeted.filler_area < 200,
            "cell-space filler complement, multiplier or halo exclusion "
            "incorrect");
      check(input.node_size_x == design.node_size_x &&
                input.node_x.size() == design.node_x.size(),
            "cell-space fillers escaped into physical placement");
    }
  }
  options.filler_macro_fraction = 1.01;
  rejects([&] { place(design, options); });
}
void gp_ab_search() {
  GpSearch policy;
  std::string why;
  check(select_gp_trial({{1., -1., 10., -1.}, {.8, -100., 1000., -1000.}},
                        policy, &why) == 1 &&
            why == "rsmt_rudy",
        "RUDY must rank before timing");
  check(select_gp_trial({{1., -10., 10., -1.}, {1.005, -1., 100., -100.}},
                        policy, &why) == 1 &&
            why == "wns",
        "WNS must decide within RUDY tie band");
  check(select_gp_trial({{1., -1., 100., -1.}, {1., -1., 80., -100.}}, policy,
                        &why) == 1 &&
            why == "hpwl",
        "HPWL must rank before TNS");
  check(select_gp_trial({{1., -1., 100., -100.}, {1., -1., 100., -1.}}, policy,
                        &why) == 1 &&
            why == "tns",
        "TNS must break final tie");
  check(select_gp_trial({{1., -1., 100., -1.}, {1., -1., 100., -1.}}, policy) ==
            0,
        "ties must keep unchanged control");
  rejects([&] { select_gp_trial({{NAN, 0, 0, 0}}, policy); });
  auto base = design();
  base.node_macro.assign(base.node_x.size(), 0);
  base.node_macro[6] = 1;
  Options o;
  o.vanished_memory = true;
  o.max_iterations = o.min_iterations = 300;
  o.bins_x = o.bins_y = 8;
  o.local_density_bins = 16;
  o.local_density_weight = .67;
  o.local_density_target = .05;
  o.local_density_jitter = 0;
  o.target_density = .5;
  o.add_fillers = true;
  o.filler_max_count = 64;
  o.legalize = true;
  o.detailed_place = false;
  o.routability = o.improved_routability = true;
  o.route.model = RouteModel::rsmt;
  o.route_capacity_h = o.route_capacity_v = 1;
  o.after_macro_legalize = [](Design &) {
  }; // host finalization hook, excluded from trials
  o.gp_search.enabled = true;
  o.gp_search.rounds = 2;
  o.num_threads = 3;
  const auto main_thread = std::this_thread::get_id();
  int timing_calls = 0, normal_calls = 0;
  double expected_weight = 1;
  o.checkpoint_interval = 20;
  o.checkpoint = [&](const Design &, int it) {
    ++normal_calls;
    expected_weight = 1 + .001 * it;
    GpUpdate u;
    u.net_weights.assign(base.net_weights.size(), expected_weight);
    return u;
  };
  std::vector<Design> scored;
  o.trial_timing = [&](const Design &d) {
    check(std::this_thread::get_id() == main_thread,
          "STA callback ran in a worker");
    check(d.net_weights == std::vector<double>(base.net_weights.size(),
                                               expected_weight) &&
              d.pin2node == base.pin2node && d.node_size_x == base.node_size_x,
          "trial changed net weights, topology or dimensions");
    ++timing_calls;
    scored.push_back(d);
    return std::pair{-hpwl(d), -2 * hpwl(d)};
  };
  std::vector<GpTrialSample> trace;
  std::optional<Design> selected;
  o.trial_observer = [&](const GpTrialSample &s) {
    check(s.iteration == 220 + int(trace.size()) * 40 && s.steps == 3,
          "incorrect branch phase or horizon");
    check(s.settings.size() >= 2 && s.settings.size() <= 3 &&
              s.settings.size() == s.scores.size(),
          "missing hold/change candidate");
    check(s.selected == select_gp_trial(s.scores, o.gp_search),
          "wrong selected trial");
    for (size_t i = 0; i < s.scores.size(); ++i) {
      const auto &physical = scored[scored.size() - s.scores.size() + i];
      RudyOptions ro{o.bins_x,           o.bins_y, o.route_capacity_h,
                     o.route_capacity_v, true,     true};
      ro.model = RudyModel::rsmt;
      const auto route = rudy(physical, ro);
      const auto &q = s.scores[i];
      check(
          q.route_unroutable_nets == route.unroutable_nets &&
              q.route_blocked_pins == route.blocked_pins &&
              q.route_measured_nets == route.measured_nets &&
              q.route_peak_utilization == route.max_utilization &&
              q.route_detour_wirelength == route.detour_wirelength &&
              std::abs(q.route_congestion_overflow +
                       double(q.route_unroutable_nets + q.route_blocked_pins) /
                           (o.bins_x * o.bins_y) -
                       q.route_overflow) < 1e-12,
          "trial route diagnostics disagree with candidate geometry");
    }
    selected = scored[scored.size() - s.scores.size() + s.selected];
    trace.push_back(s);
  };
  o.observer = [&](const Design &d, PlacementStage stage, int it) {
    if (selected && stage == PlacementStage::global_placement &&
        (it == 223 || it == 263)) {
      check(d.node_x == selected->node_x && d.node_y == selected->node_y,
            "winner coordinates were not adopted exactly");
    }
    for (size_t i = base.num_movable; i < base.node_x.size(); ++i)
      check(d.node_x[i] == base.node_x[i] && d.node_y[i] == base.node_y[i],
            "trial moved a fixed node");
  };
  auto first = base;
  const auto result = place(first, o);
  check(result.legal && result.iterations == 300 && trace.size() == 2 &&
            timing_calls >= 6,
        "A/B search failed final legality or schedule");
  check(normal_calls == 15, "losing trials leaked checkpoint updates");
  const auto first_trace = trace;
  trace.clear();
  scored.clear();
  selected.reset();
  auto second = base;
  place(second, o);
  check(first.node_x == second.node_x && first.node_y == second.node_y,
        "parallel A/B trials are not repeatable");
  for (size_t i = 0; i < trace.size(); ++i)
    check(trace[i].selected == first_trace[i].selected,
          "A/B selection is nondeterministic");
  // Unlimited search must continue past the old 12-decision cap, but leave
  // at least 20 iterations after the last trial for final settling.
  trace.clear();
  scored.clear();
  selected.reset();
  o.gp_search.rounds = 0;
  o.max_iterations = o.min_iterations = 1000;
  auto unlimited = base;
  const auto unlimited_result = place(unlimited, o);
  check(unlimited_result.legal && unlimited_result.iterations == 1000 &&
            trace.size() == 19 && trace.back().iteration == 940 &&
            trace.back().iteration + trace.back().steps + 20 <= 1000,
        "unlimited A/B search stopped early or consumed final settling");
  o.gp_search.rounds = -1;
  auto invalid_rounds = base;
  rejects([&] { place(invalid_rounds, o); });
  o.gp_search.rounds = 2;
  o.max_iterations = o.min_iterations = 300;
  // Trial timing must use the same physical coordinates as normal checkpoints,
  // even when the optimizer reserves halo padding around the cells.
  auto haloed = base;
  haloed.node_halo_x.assign(base.node_x.size(), 0);
  haloed.node_halo_y.assign(base.node_x.size(), 0);
  haloed.node_halo_x[0] = 1;
  haloed.node_halo_y[0] = 1;
  auto halo_options = o;
  halo_options.num_threads = 1;
  halo_options.trial_observer = {};
  halo_options.observer = {};
  int halo_scores = 0;
  halo_options.trial_timing = [&](const Design &physical) {
    ++halo_scores;
    check(physical.node_size_x == base.node_size_x &&
              physical.node_size_y == base.node_size_y &&
              physical.pin_offset_x == base.pin_offset_x &&
              physical.pin_offset_y == base.pin_offset_y &&
              physical.node_pad_x.empty() && physical.node_pad_y.empty(),
          "trial timing received halo-padded geometry");
    return std::pair{-hpwl(physical), -2 * hpwl(physical)};
  };
  place(haloed, halo_options);
  check(halo_scores >= 6, "halo trial timing was not exercised");
  // The warm-up also gates speculative timing, not only normal checkpoints.
  auto delayed = base;
  auto delayed_options = o;
  delayed_options.checkpoint_start = 260;
  delayed_options.observer = {};
  int delayed_trials = 0, delayed_checkpoints = 0, delayed_scores = 0;
  delayed_options.checkpoint = [&](const Design &, int iteration) {
    check(iteration >= 260, "normal timing bypassed warm-up");
    ++delayed_checkpoints;
    return GpUpdate{};
  };
  delayed_options.trial_timing = [&](const Design &) {
    ++delayed_scores;
    return std::pair{0., 0.};
  };
  delayed_options.trial_observer = [&](const GpTrialSample &sample) {
    check(sample.iteration >= 260, "trial timing bypassed warm-up");
    ++delayed_trials;
  };
  place(delayed, delayed_options);
  check(delayed_trials == 1 && delayed_checkpoints > 0 && delayed_scores > 0,
        "warm-up disabled timing instead of delaying it");
  auto topology = base;
  auto topology_options = o;
  topology_options.checkpoint = [](const Design &d, int) {
    GpUpdate update;
    update.replacement = d;
    return update;
  };
  rejects([&] { place(topology, topology_options); });
  check(topology.node_x == base.node_x && topology.pin2node == base.pin2node,
        "rejected A/B topology update changed caller state");
  auto invalid_weight_options = o;
  invalid_weight_options.local_density_weight = 2;
  auto invalid_weight = base;
  rejects([&] { place(invalid_weight, invalid_weight_options); });
  auto failed = base;
  o.trial_timing = [](const Design &) -> std::pair<double, double> {
    throw std::runtime_error("host timing failed");
  };
  rejects([&] { place(failed, o); });
  check(failed.node_x == base.node_x && failed.node_y == base.node_y &&
            failed.net_weights == base.net_weights,
        "failed trial scoring partially committed caller state");
  auto bad = base;
  o.gp_search.steps = 4;
  rejects([&] { place(bad, o); });
  check(bad.node_x == base.node_x, "invalid A/B search changed caller state");
}
void route_policy_helpers() {
  for (auto model : {RouteModel::off, RouteModel::rudy, RouteModel::rsmt,
                     RouteModel::automatic, RouteModel::uroute})
    check(route_model_from_name(route_model_name(model)) == model,
          "route model name does not round-trip");
  rejects([] { route_model_from_name("automatic"); });
  RoutePolicy p;
  validate_route_policy(p);
  check(!route_latches(p) && !uroute_requested(p), "default policy latches");
  p.probe = true;
  check(route_latches(p) && !uroute_requested(p), "RUDY probe must latch");
  p.model = RouteModel::rsmt;
  check(!route_latches(p), "explicit RSMT probe must not latch");
  p.probe_model = RouteModel::uroute;
  check(route_latches(p) && uroute_requested(p), "URoute probe must latch");
  p = {};
  p.uroute_ladder = true;
  check(!uroute_requested(p), "ladder needs the automatic model");
  p.model = RouteModel::automatic;
  check(route_latches(p) && uroute_requested(p), "automatic ladder");
  p.uroute_switch_overflow = p.switch_overflow * 2;
  rejects([&] { validate_route_policy(p); });
  p = {};
  p.probe_model = RouteModel::rudy;
  rejects([&] { validate_route_policy(p); });
}
int main() {
  try {
    density_plateau_hands_off_without_claiming_convergence();
    setup_timing_does_not_veto_convergence();
    route_policy_helpers();
    size_neutral_macro_pull();
    macro_sized_virtual_fillers();
    gp_ab_search();
    virtual_macro_pulls();
    spectral();
    power2_grid();
    timing_spreading();
    far_fixed_nodes();
    threaded_determinism();
    gradients();
    hpwl_case();
    native_flow();
    transparent_macros();
    vanished_memories();
    moving_cells();
    moving_terminals();
    offgrid_macro_swap();
    same_bin_swap();
    weighted_detailed_placement();
    weighted_detailed_placement_wide_die();
    placement_masks();
    placer_instances();
    concurrent_instances();
    macros();
    zero_area_fixed_terminals<float>();
    zero_area_fixed_terminals<double>();
    fractional_io_detailed_placement();
    continued_fixed_terminals();
    mixed_sizes();
    relocated_macro_obstacles();
    halos();
    failures_and_empty();
    abacus_weights();
    after_macro_legalize_writeback();
    macro_only_legalization();
  } catch (const std::exception &e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  std::cout << "native placement checks passed\n";
}
