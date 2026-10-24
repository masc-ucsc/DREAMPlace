// BSD-3-Clause. See LICENSE.
#include "local_density.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool v) {
  if (!v)
    throw std::runtime_error("local density check failed");
}
int main() {
  Design d;
  d.layout = {0, 0, 64, 64, 1, 1};
  d.num_movable = 2;
  d.node_x = {29.5, 33.5};
  d.node_y = {31.5, 31.5};
  d.node_size_x = {1, 1};
  d.node_size_y = {1, 1};
  d.node_macro = {0, 0};
  d.build_maps();
  auto g = detail::LocalDensity(d, 64, 4).gradient(d);
  check(g.x[0] > 0 &&
        g.x[1] <
            0); // descent moves the cells apart despite no physical overlap
  check(std::abs(g.x[0] + g.x[1]) < 1e-9 && std::abs(g.y[0]) < 1e-9);
  // Optional macro metadata is not a reliable physical-node count. Keep the
  // local source/receivers unchanged when global GP appends virtual particles.
  auto no_metadata = d;
  no_metadata.node_macro.clear();
  detail::LocalDensity real_field(no_metadata, 64, 4);
  const auto real_source = real_field.cell_density(no_metadata);
  const auto real_gradient = real_field.gradient(no_metadata);
  auto particles = no_metadata;
  particles.node_x.push_back(32);
  particles.node_y.push_back(32);
  particles.node_size_x.push_back(20);
  particles.node_size_y.push_back(20);
  particles.build_maps();
  check(real_field.cell_density(particles) == real_source);
  const auto particle_gradient = real_field.gradient(particles);
  check(particle_gradient.x[0] == real_gradient.x[0] &&
        particle_gradient.y[1] == real_gradient.y[1] &&
        particle_gradient.x[2] == 0 && particle_gradient.y[2] == 0);
  auto far = d;
  far.node_x = {21.5, 41.5};
  auto f = detail::LocalDensity(far, 64, 4).gradient(far);
  check(std::abs(f.x[0]) < std::abs(g.x[0]));
  auto quiet = detail::LocalDensity(d, 64, 4, 10).gradient(d);
  check(quiet.x[0] == 0 && quiet.y[1] == 0);
  auto active = detail::LocalDensity(d, 64, 4, .1).gradient(d);
  check(active.x[0] > 0 && active.x[1] < 0);
  // An independent local threshold can activate spreading below the global
  // target. Inheritance must retain the old trajectory exactly.
  Options threshold;
  threshold.add_fillers = threshold.legalize = threshold.detailed_place = false;
  threshold.max_iterations = threshold.min_iterations = 1;
  threshold.bins_x = threshold.bins_y = threshold.local_density_bins = 64;
  threshold.target_density = .5;
  threshold.local_density_weight = .75;
  threshold.local_density_jitter = 0;
  auto inherited_design = d, high_design = d, low_design = d;
  const auto inherited_result = place(inherited_design, threshold);
  threshold.local_density_target = .5;
  const auto high_result = place(high_design, threshold);
  threshold.local_density_target = .05;
  const auto low_result = place(low_design, threshold);
  check(inherited_design.node_x == high_design.node_x &&
        inherited_design.node_y == high_design.node_y);
  check(inherited_result.local_density_force == 0 &&
        high_result.local_density_force == 0 &&
        low_result.local_density_force > 0);
  // A transparent fixed macro must leave local capacity identical to a
  // zero-area terminal, while opaque capacity remains blocked.
  auto obstacle = d;
  obstacle.node_x.push_back(26);
  obstacle.node_y.push_back(26);
  obstacle.node_size_x.push_back(14);
  obstacle.node_size_y.push_back(14);
  obstacle.node_macro.push_back(1);
  obstacle.build_maps();
  auto terminals = obstacle;
  terminals.node_size_x.back() = terminals.node_size_y.back() = 0;
  detail::LocalDensity obstacle_field(obstacle, 64, 4, .1);
  const auto visible = obstacle_field.gradient(obstacle);
  const auto transparent = obstacle_field.gradient(obstacle, false, {}, true);
  const auto zero_body =
      detail::LocalDensity(terminals, 64, 4, .1).gradient(terminals);
  check(transparent.x == zero_body.x && transparent.y == zero_body.y);
  check(visible.x != transparent.x || visible.y != transparent.y);
  check(transparent.x.back() == 0 && transparent.y.back() == 0);
  // A half-sized obstacle must be the same field as an explicitly centered
  // half-sized rectangle, not a half-charge full-sized rectangle.
  auto half = obstacle;
  half.node_x.back() += 3.5;
  half.node_y.back() += 3.5;
  half.node_size_x.back() *= .5;
  half.node_size_y.back() *= .5;
  const auto grown = obstacle_field.gradient(obstacle, false, {}, false, .5);
  const auto actual_half = obstacle_field.gradient(half);
  check(grown.x == actual_half.x && grown.y == actual_half.y);
  const auto overlap_grown = detail::fixed_macro_overlap(obstacle, 8, 8, 1, .5);
  const auto overlap_half = detail::fixed_macro_overlap(half, 8, 8, 1);
  check(overlap_grown.x == overlap_half.x && overlap_grown.y == overlap_half.y);
  auto macro = d;
  macro.node_macro = {1, 1};
  auto m = detail::LocalDensity(macro, 64, 4).gradient(macro);
  check(m.x[0] == 0 && m.y[1] == 0);
  // Only standard cells emit the macro field. A large halo-sized box must
  // feel nearby cells at its edges even when its center is many bins away.
  Design pressure;
  pressure.layout = {0, 0, 128, 128, 1, 1};
  pressure.num_movable = 1;
  pressure.node_x = {40};
  pressure.node_y = {48};
  pressure.node_size_x = {32};
  pressure.node_size_y = {32};
  pressure.node_macro = {1};
  for (int k = 0; k < 8; ++k) {
    pressure.node_x.push_back(73);
    pressure.node_y.push_back(56.5 + 2 * k);
    pressure.node_size_x.push_back(1);
    pressure.node_size_y.push_back(1);
    pressure.node_macro.push_back(0);
  }
  pressure.build_maps();
  const auto pm =
      detail::LocalDensity(pressure, 128, 4, .5).gradient(pressure, true);
  check(pm.x[0] > 0 && std::abs(pm.y[0]) < 1e-8); // descent moves left
  for (int i = 1; i < int(pm.x.size()); ++i)
    check(pm.x[i] == 0 && pm.y[i] == 0); // recipient macros only
  auto left = pressure;
  for (int i = 1; i < int(left.node_x.size()); ++i)
    left.node_x[i] = 38;
  const auto pl = detail::LocalDensity(left, 128, 4, .5).gradient(left, true);
  check(pl.x[0] < 0);
  auto distant = pressure;
  for (int i = 1; i < int(distant.node_x.size()); ++i)
    distant.node_x[i] = 100;
  const auto pd =
      detail::LocalDensity(distant, 128, 4, .5).gradient(distant, true);
  check(std::abs(pd.x[0]) < std::abs(pm.x[0]));
  // A tiny amount of distant cell area must push a much larger memory.
  // Two largest-memory side lengths give a broad field, independently of
  // the short-range cell-to-cell solver. Sources accumulate without a
  // density threshold or a requirement to match the SRAM's area.
  const double long_reach = detail::macro_cell_screening_length(distant, 2);
  check(long_reach == 64);
  detail::LocalDensity broad(distant, 128, long_reach, .5);
  const auto long_force = broad.gradient(distant, true);
  check(long_force.x[0] > 2 * pd.x[0]);
  auto twice_source = broad.cell_density(distant);
  for (auto &v : twice_source)
    v *= 2;
  const auto twice_force = broad.gradient(distant, true, twice_source);
  check(std::abs(twice_force.x[0] - 2 * long_force.x[0]) < 1e-8);
  for (size_t i = 1; i < long_force.x.size(); ++i)
    check(long_force.x[i] == 0 && long_force.y[i] == 0);
  auto padded = distant;
  padded.node_pad_x.assign(padded.node_x.size(), 0);
  padded.node_pad_y.assign(padded.node_x.size(), 0);
  padded.node_pad_x[0] = 3;
  padded.node_pad_y[0] = 7;
  padded.node_size_x[0] += 6;
  padded.node_size_y[0] += 14;
  check(detail::macro_cell_screening_length(padded, 2) == long_reach);
  auto no_cells = pressure;
  no_cells.node_macro.assign(no_cells.node_x.size(), 1);
  const auto pn =
      detail::LocalDensity(no_cells, 128, 4, .5).gradient(no_cells, true);
  check(pn.x[0] == 0 && pn.y[0] == 0); // no macro-to-macro field
  auto locked = pressure;
  locked.num_movable = 0;
  const auto pf =
      detail::LocalDensity(locked, 128, 4, .5).gradient(locked, true);
  check(pf.x[0] == 0 && pf.y[0] == 0);
  // Prefix-sum pull is unbiased on uniform density, points down a density
  // slope, and never uses off-core space or macros as an empty-space source.
  detail::LocalDensity local(pressure, 128, 4, .5);
  std::vector<double> uniform(128 * 128, .25), slope(uniform.size());
  auto uniform_pull = local.density_pull(pressure, uniform, 4);
  check(std::abs(uniform_pull.x[0]) < 1e-9 &&
        std::abs(uniform_pull.y[0]) < 1e-9);
  for (int x = 0; x < 128; ++x)
    for (int y = 0; y < 128; ++y)
      slope[x * 128 + y] = (x + .5) / 128;
  const auto uphill = local.density_pull(pressure, slope, 4);
  check(uphill.x[0] > 0 && std::abs(uphill.y[0]) < 1e-9);
  for (int i = 1; i < int(uphill.x.size()); ++i)
    check(uphill.x[i] == 0 && uphill.y[i] == 0);
  for (auto &value : slope)
    value = 1 - value;
  check(local.density_pull(pressure, slope, 4).x[0] < 0);
  auto edge = pressure;
  edge.node_x[0] = 0;
  check(local.density_pull(edge, uniform, 4).x[0] == 0);
  edge.node_x[0] = .3;
  check(std::abs(local.density_pull(edge, uniform, 4).x[0]) < 1e-9);
  const auto raw = local.cell_density(pressure);
  const auto recomputed = local.gradient(pressure, true);
  const auto reused = local.gradient(pressure, true, raw);
  check(recomputed.x == reused.x && recomputed.y == reused.y);
  check(local.density_pull(no_cells, local.cell_density(no_cells), 4).x[0] ==
        0);
  check(local.density_pull(locked, uniform, 4).x[0] == 0);
  auto fixed = d;
  fixed.num_movable = 1;
  auto fix = detail::LocalDensity(fixed, 64, 4).gradient(fixed);
  check(fix.x[1] == 0);
  Options o;
  o.add_fillers = o.legalize = o.detailed_place = false;
  o.max_iterations = 40;
  auto a = d, b = d;
  auto ra = place(a, o);
  o.local_density_weight = 0;
  auto rb = place(b, o);
  check(a.node_x == b.node_x && a.node_y == b.node_y && ra.hpwl == rb.hpwl);
  o.local_density_weight = 2;
  o.local_density_bins = 64;
  o.local_density_jitter = 0;
  auto c = d;
  place(c, o);
  check(c.node_x[1] - c.node_x[0] >= d.node_x[1] - d.node_x[0]);
  check(detail::dynamic_bins(d, .5, 512) == std::pair{32, 32});
  auto wide = d;
  wide.layout.xh = 128;
  const auto bins = detail::dynamic_bins(wide, .5, 512);
  check(bins.first == 2 * bins.second);
  auto capped = d;
  capped.layout.xh = capped.layout.yh = 10000;
  check(detail::dynamic_bins(capped, .5, 256) == std::pair{256, 256});
  // Frozen coordinates must survive both new controls.
  o.dynamic_bins = true;
  o.bin_max = 64;
  auto fi = fixed;
  place(fi, o);
  check(fi.node_x[1] == fixed.node_x[1] && fi.node_y[1] == fixed.node_y[1]);
  Options popt;
  popt.add_fillers = popt.legalize = popt.detailed_place = false;
  popt.local_density_bins = 128;
  popt.max_iterations = popt.min_iterations = 20;
  auto baseline = pressure, disabled = pressure, pulling = pressure;
  const auto before = place(baseline, popt);
  popt.macro_density_pull_weight = 0;
  const auto off = place(disabled, popt);
  check(baseline.node_x == disabled.node_x &&
        baseline.node_y == disabled.node_y && before.hpwl == off.hpwl);
  auto pressure_a = pressure, pressure_b = pressure;
  auto cell_pressure = popt;
  cell_pressure.macro_cell_density_weight = .25;
  const auto pressure_result = place(pressure_a, cell_pressure);
  cell_pressure.local_density_range = 1;
  const auto independent = place(pressure_b, cell_pressure);
  check(pressure_a.node_x == pressure_b.node_x &&
        pressure_a.node_y == pressure_b.node_y);
  check(pressure_result.macro_cell_density_screening_length == 64 &&
        pressure_result.macro_cell_density_coefficient > 0 &&
        pressure_result.macro_cell_density_force > 0 &&
        independent.macro_cell_density_force ==
            pressure_result.macro_cell_density_force);
  popt.macro_density_pull_weight = .1;
  const auto on = place(pulling, popt);
  check(on.macro_density_pull_seconds > 0 && on.macro_density_pull_force > 0);
  check(pulling.node_x[1] == pressure.node_x[1] &&
        pulling.node_y[1] == pressure.node_y[1]);
  auto anchored = pressure;
  CentroidGroup anchored_group;
  anchored_group.name = "anchored";
  anchored_group.nodes = {0};
  anchored.centroid_groups = {anchored_group};
  anchored.centroid_targets = {{56, 64}};
  anchored.node_centroid_color.assign(anchored.node_x.size(), -1);
  for (size_t i = 0; i < anchored.node_x.size(); ++i)
    anchored.centroid_node_area.push_back(anchored.node_size_x[i] *
                                          anchored.node_size_y[i]);
  popt.centroid_weight = .25;
  popt.centroid_distance_attraction = true;
  const auto bounded = place(anchored, popt);
  check(bounded.centroid.attraction_coefficient > 0);
  check(bounded.macro_density_pull_anchor_ratio > 0 &&
        bounded.macro_density_pull_anchor_ratio <=
            popt.macro_density_pull_weight + 1e-12);
  // The force must actually spread a connected overlapping cloud in GP.
  Design cloud;
  cloud.layout = {0, 0, 64, 64, 1, 1};
  cloud.num_movable = 24;
  for (int i = 0; i < 25; ++i) {
    cloud.node_x.push_back(31.5);
    cloud.node_y.push_back(31.5);
    cloud.node_size_x.push_back(i == 24 ? 0 : 1);
    cloud.node_size_y.push_back(i == 24 ? 0 : 1);
    cloud.node_macro.push_back(0);
  }
  for (int i = 0; i < 24; ++i) {
    cloud.pin2node.insert(cloud.pin2node.end(), {i, 24});
    cloud.pin_offset_x.insert(cloud.pin_offset_x.end(), {.5, 0});
    cloud.pin_offset_y.insert(cloud.pin_offset_y.end(), {.5, 0});
    cloud.flat_net2pin.insert(cloud.flat_net2pin.end(), {2 * i, 2 * i + 1});
    cloud.flat_net2pin_start.push_back(2 * i + 2);
  }
  cloud.build_maps();
  Options copt;
  copt.add_fillers = copt.legalize = copt.detailed_place = false;
  copt.max_iterations = copt.min_iterations = 100;
  copt.local_density_bins = 64;
  copt.target_density = .5;
  auto packed = cloud, spread = cloud;
  place(packed, copt);
  copt.local_density_weight = 4;
  place(spread, copt);
  const auto radius = [](const Design &v) {
    double sum = 0;
    for (int i = 0; i < v.num_movable; ++i)
      sum += std::hypot(v.node_x[i] - 31.5, v.node_y[i] - 31.5);
    return sum / v.num_movable;
  };
  check(radius(spread) > radius(packed) + .5);
  // The independent target inherits exactly by default and a lower target
  // must increase local pressure on the same connected cloud.
  auto inherited = cloud;
  copt.local_density_target = copt.target_density;
  place(inherited, copt);
  check(inherited.node_x == spread.node_x && inherited.node_y == spread.node_y);
  auto looser = cloud;
  copt.local_density_target = .05;
  place(looser, copt);
  check(radius(looser) > radius(spread));
  auto combined = obstacle;
  Options trans;
  trans.bins_x = trans.bins_y = trans.local_density_bins = 64;
  trans.max_iterations = trans.min_iterations = 40;
  trans.local_density_weight = 2;
  trans.local_density_target = .05;
  trans.local_density_jitter = 0;
  trans.macro_transparency.always = true;
  trans.add_fillers = trans.legalize = trans.detailed_place = false;
  const auto combined_result = place(combined, trans);
  trans.macro_transparency.always = false;
  auto counterpart = terminals;
  const auto counterpart_result = place(counterpart, trans);
  check(combined.node_x == counterpart.node_x &&
        combined.node_y == counterpart.node_y);
  check(combined_result.overflow == counterpart_result.overflow);
  std::cout << "local density checks passed\n";
}
