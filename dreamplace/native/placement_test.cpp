// BSD-3-Clause. See LICENSE. Checks deliberately remain active under -c opt.
#include "dreamplace/native/placement.hpp"

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
  // Matching only solves sets of >= 3 equal-size cells sharing a bin, and (upstream
  // quirk) needs two such sets to run at all. Cells 0-3 are the 1-wide set that must
  // swap; cells 4-6 are a 2-wide filler set. One cell per row keeps k_reorder from
  // doing the swap first.
  d.layout = {0, 0, w, 7, 1, 1};
  d.num_movable = 7;
  d.node_x = {1, 3, 5, 7, 1, 4, 7, 0, w};
  d.node_y = {0, 1, 2, 3, 4, 5, 6, 0, 0};
  d.node_size_x = {1, 1, 1, 1, 2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 1, 1, 1, 1, 0, 0};
  // net0: cell0-right (heavy), net1: cell3-left (heavy), net2: cell1-left, net3: cell2-right
  d.pin2node = {0, 8, 3, 7, 1, 7, 2, 8};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.net_weights = {1000, 1000, 1, 1};  // weight * die span exceeds INT_MAX unscaled
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
// depend on thread arrival order. Repeated identical runs must agree exactly for
// a fixed worker count, including when OpenMP is enabled (--config=openmp).
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
// Fixed nodes only need finite coordinates; one far above the core used to turn
// an out-of-range double into a negative row index inside is_legal().
void far_fixed_nodes() {
  for (const double far : {1e12, -1e12, 3e9}) {
    Design d;
    d.layout = {0, 0, 40, 20, 1, 1};
    d.num_movable = 4;
    d.node_x = {2, 5, 8, 20, far, 30};
    d.node_y = {2, 6, 2, 9, far, 5};
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
void power2_grid() {
  // ISM sizes its bin grid with ceil_power2(bins / 8); fewer than eight bins
  // gives zero, whose log2 used to be cast to int and used as a shift count.
  using DreamPlace::ceil_power2;
  const int expected[] = {1, 1, 2, 4, 4, 8, 8, 8, 8, 16};
  for (int v = 0; v < 10; ++v)
    check(ceil_power2(v) == expected[v], "ceil_power2 mismatch");
  check(ceil_power2(-3) == 1, "ceil_power2 of a negative count");
}
} // namespace
int main() {
  try {
    spectral();
    power2_grid();
    far_fixed_nodes();
    threaded_determinism();
    gradients();
    hpwl_case();
    native_flow();
    offgrid_macro_swap();
    same_bin_swap();
    weighted_detailed_placement();
    weighted_detailed_placement_wide_die();
    placer_instances();
    concurrent_instances();
    macros();
    continued_fixed_terminals();
    mixed_sizes();
    relocated_macro_obstacles();
    halos();
    failures_and_empty();
  } catch (const std::exception &e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  std::cout << "native placement checks passed\n";
}
