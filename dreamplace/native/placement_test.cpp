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

using namespace DreamPlace::native;
namespace {
void check(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
void near(double a, double b, double eps, const char* message) { check(std::abs(a - b) <= eps, message); }
template <typename F>
void rejects(F f) {
  bool rejected = false;
  try {
    f();
  } catch (const std::exception&) {
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
      rho[x * ny + y] =
          std::cos(std::numbers::pi * 2 * (x + 0.5) / nx) * std::cos(std::numbers::pi * 3 * (y + 0.5) / ny);
  s.solve(rho, phi, fx, fy);
  for (int x = 0; x < nx; ++x)
    for (int y = 0; y < ny; ++y) {
      double ax = std::numbers::pi * 2 * (x + 0.5) / nx, ay = std::numbers::pi * 3 * (y + 0.5) / ny;
      near(phi[x * ny + y], rho[x * ny + y] / (wx * wx + wy * wy), 1e-12, "Poisson eigenmode");
      near(fx[x * ny + y], wx * std::sin(ax) * std::cos(ay) / (wx * wx + wy * wy), 1e-12, "Poisson x field");
      near(fy[x * ny + y], wy * std::cos(ax) * std::sin(ay) / (wx * wx + wy * wy), 1e-12, "Poisson y field");
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
      auto& pos = axis ? d.node_y : d.node_x;
      double old = pos[i], eps = 1e-5;
      pos[i] = old + eps;
      double plus = detail::wirelength(d, 1.7, 1).value;
      pos[i] = old - eps;
      double minus = detail::wirelength(d, 1.7, 1).value;
      pos[i] = old;
      near((plus - minus) / (2 * eps), axis ? g.y[i] : g.x[i], 2e-7, "weighted WA gradient");
    }
  // Large absolute coordinates must not cause exp() overflow.
  for (auto& x : d.node_x) x += 1e6;
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
  check(r.iterations > 0 && r.iterations <= o.max_iterations, "iteration accounting");
  check(std::isfinite(r.hpwl) && std::isfinite(r.overflow), "finite placement metrics");
  check(d.node_x.size() == before.node_x.size() && d.pin2node == before.pin2node, "fillers must not escape");
  for (size_t i = d.num_movable; i < d.node_x.size(); ++i) {
    check(d.node_x[i] == before.node_x[i] && d.node_y[i] == before.node_y[i], "fixed obstacle/terminal moved");
  }
  auto again = before;
  auto rr = place(again, o);
  check(d.node_x == again.node_x && d.node_y == again.node_y && rr.hpwl == r.hpwl, "deterministic placement");
  auto observed = before;
  Options traced = o;
  traced.observer_interval = 7;  // Deliberately not aligned with 20-step convergence checks.
  std::vector<PlacementStage> stages;
  Design last;
  int previous_iteration = 0;
  traced.observer = [&](const Design& frame, PlacementStage stage, int iteration) {
    frame.validate();
    check(frame.node_x.size() == before.node_x.size() && frame.pin2node == before.pin2node,
          "observer must see physical nodes/connectivity without fillers");
    for (size_t i = frame.num_movable; i < frame.node_x.size(); ++i)
      check(frame.node_x[i] == before.node_x[i] && frame.node_y[i] == before.node_y[i], "observer fixed nodes moved");
    check(iteration >= previous_iteration, "observer iteration order");
    previous_iteration = iteration;
    if (stage == PlacementStage::initial)
      check(frame.node_x == before.node_x && frame.node_y == before.node_y && iteration == 0, "observer initial state");
    stages.push_back(stage);
    last = frame;
  };
  auto tr = place(observed, traced);
  check(
      observed.node_x == d.node_x && observed.node_y == d.node_y && tr.hpwl == r.hpwl && tr.iterations == r.iterations,
      "observation must not change placement");
  check(stages.front() == PlacementStage::initial && stages.back() == PlacementStage::finished &&
            stages[stages.size() - 2] == PlacementStage::detailed_placement &&
            stages[stages.size() - 3] == PlacementStage::legalization &&
            stages[stages.size() - 4] == PlacementStage::global_placement,
        "observer phase sequence");
  check(last.node_x == observed.node_x && last.node_y == observed.node_y, "observer final state");
  traced.observer = [](const Design&, PlacementStage stage, int) {
    if (stage == PlacementStage::finished) throw std::runtime_error("observer abort");
  };
  auto aborted = before;
  rejects([&] { place(aborted, traced); });
  check(aborted.node_x == before.node_x && aborted.node_y == before.node_y, "observer abort must not write back");
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
  check(gp.node_x != before.node_x || gp.node_y != before.node_y, "GP should move cells");
  auto lg = before;
  o.legalize = true;
  o.detailed_place = true;
  o.max_iterations = 0;
  check(place(lg, o).legal, "legalization-only path");
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
  const Placer& configured = global;
  auto d = design();
  check(configured.place(d).iterations == 17, "placer instance must own its options");
  check(configured.options().max_iterations == 17 && legalizer.options().legalize &&
            legalizer.options().max_iterations == 0,
        "separate placer instances must have independent configuration");
  check(legalizer.place(d).legal && is_legal(d), "instance must support legalizing a caller-owned placement");
  // Reconfigure and reuse one object; failed runs must preserve caller data.
  global.options().max_iterations = -1;
  auto before = d;
  rejects([&] { global.place(d); });
  check(d.node_x == before.node_x && d.node_y == before.node_y, "invalid instance options must not change design");
  global.options() = legalizer.options();
  check(global.place(d).legal, "placer instance must allow reconfiguration and reuse");
  check(legalizer.options().max_iterations == 0, "reconfiguring another instance must not change its source");
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
    auto& d = inputs[i];
    double offset = 100 * i;
    d.layout.xl += offset;
    d.layout.xh += offset;
    for (auto& x : d.node_x) x += offset;
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
        placer.options().observer = [&](const Design& frame, PlacementStage stage, int) {
          if (stage == PlacementStage::initial) start.arrive_and_wait();
          check(frame.node_x.size() == inputs[i].node_x.size(), "parallel snapshot leaked another run's nodes");
          for (size_t j = frame.num_movable; j < frame.node_x.size(); ++j)
            check(frame.node_x[j] == inputs[i].node_x[j] && frame.node_y[j] == inputs[i].node_y[j],
                  "parallel run moved another instance's fixed node");
        };
        auto d = inputs[i];
        auto result = placer.place(d);
        return Run{std::move(d), result};
      });
    for (size_t i = 0; i < count; ++i) {
      auto run = jobs[i].get();
      check(run.result.legal && is_legal(run.design) && std::isfinite(run.result.hpwl),
            "parallel instances must produce legal finite placements");
      if (threads == 1)
        check(run.design.node_x == reference[i].design.node_x && run.design.node_y == reference[i].design.node_y &&
                  run.result.hpwl == reference[i].result.hpwl &&
                  run.result.iterations == reference[i].result.iterations,
              "independent parallel calls must match their serial references");
    }
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
  // Real pins with no net, duplicate pins, nonmonotonic CSR and bad maps fail before kernels.
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
  check(d.node_x.back() == original.node_x.back() && d.node_y.back() == original.node_y.back(),
        "macro moved fixed obstacle");
  auto shifted = original;
  shifted.layout = {-50, -40, -10, -20, 1, 1};
  for (auto& x : shifted.node_x) x -= 50;
  for (auto& y : shifted.node_y) y -= 40;
  check(place(shifted, o).legal, "macro legalization with negative core origin");
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
  check(d.node_x == before.node_x && d.node_y == before.node_y, "failed legalization mutated design");
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
}  // namespace
int main() {
  try {
    spectral();
    gradients();
    hpwl_case();
    native_flow();
    placer_instances();
    concurrent_instances();
    macros();
    failures_and_empty();
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  std::cout << "native placement checks passed\n";
}
