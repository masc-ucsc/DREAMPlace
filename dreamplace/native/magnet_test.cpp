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
    throw std::runtime_error("magnet assertion at line " +
                             std::to_string(at.line()));
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 200, 200, 1, 1};
  d.num_movable = 7;
  d.node_x = {20, 40, 20, 40, 33, 33, 80};
  d.node_y = {20, 20, 80, 80, 20, 80, 80};
  d.node_size_x = {10, 10, 10, 10, 2, 2, 100};
  d.node_size_y = {10, 10, 10, 10, 2, 2, 100};
  d.node_macro = {1, 1, 1, 1, 0, 0, 0};
  d.node_halo_x = d.node_halo_y = {2, 2, 2, 2, 0, 0, 0};
  // Different masters/classes must still pair when the outer shapes match.
  d.node_align_class = {0, 1, 2, 3, -1, -1, -1};
  d.build_maps();
  return d;
}
MagnetSearch search(const Design &d, MagnetOptions o = {}) {
  return magnet_candidates(d, o, MagnetDensity(d, 100, 100));
}
void density_and_selection() {
  const auto d = fixture();
  const MagnetDensity grid(d, 100, 100);
  check(std::abs(grid.density({0, 0, 200, 200}) - 10008. / 40000) < 1e-12);
  check(std::abs(grid.density({32, 20, 36, 22}) - .5) < 1e-12);
  const auto a = search(d), b = search(d);
  check(a.merges.size() == 2 && b.merges.size() == 2 && a.exact_density);
  for (int i = 0; i < 2; ++i) {
    check(a.merges[i].members == b.merges[i].members);
    check(a.merges[i].poses == b.merges[i].poses);
  }
  check(a.merges[0].members == (std::vector<int>{0, 1}));
  check(a.merges[1].members == (std::vector<int>{2, 3}));
  const auto selected = std::find_if(a.candidates.begin(), a.candidates.end(),
                                     [](const auto &c) { return c.selected; });
  check(selected != a.candidates.end());
  const auto &c = *selected;
  check(c.axis == Axis::x && c.gap == 6);
  check(std::abs(c.gap_density - 4. / 84) < 1e-12);
  check(c.merge.poses[0] == std::pair{23., 20.} &&
        c.merge.poses[1] == std::pair{37., 20.});
  check(std::abs(grid.global_cell_density() - 10008. / 39600) < 1e-12);
  check(a.density_limit == .25 * grid.global_cell_density());
  auto duplicate_macro = d;
  duplicate_macro.node_x[1] = duplicate_macro.node_x[0];
  check(std::abs(MagnetDensity(duplicate_macro, 10, 10).global_cell_density() -
                 10008. / 39700) < 1e-12);
  auto clipped_macro = d;
  clipped_macro.node_x[0] = -5;
  check(std::abs(MagnetDensity(clipped_macro, 10, 10).global_cell_density() -
                 10008. / 39650) < 1e-12);
  auto fractional = d;
  for (size_t i = 0; i < fractional.node_x.size(); ++i) {
    fractional.node_x[i] += .123456789;
    fractional.node_y[i] += .987654321;
  }
  check(search(fractional).merges.size() == 2);
  MagnetOptions o;
  o.density_fraction = 1. / 16;
  check(search(d, o).merges.empty());
  o = {};
  o.max_per_round = 1;
  check(search(d, o).merges.size() == 1);
  o = {};
  o.max_gap = .4;
  check(search(d, o).merges.empty());
  auto shape = d;
  shape.node_size_x[1] = 11;
  check(search(shape).merges.size() == 1);
  auto offset = d;
  offset.node_y[1] += 7; // half the padded orthogonal extent is excluded
  check(search(offset).merges.size() == 1);
  auto obstacle = d;
  obstacle.node_x.push_back(30);
  obstacle.node_y.push_back(18);
  obstacle.node_size_x.push_back(3);
  obstacle.node_size_y.push_back(3);
  obstacle.node_macro.push_back(1);
  obstacle.node_halo_x.push_back(0);
  obstacle.node_halo_y.push_back(0);
  obstacle.node_align_class.push_back(-1);
  obstacle.build_maps();
  const auto blocked = search(obstacle);
  check(blocked.merges.size() == 1);
  check(std::any_of(blocked.candidates.begin(), blocked.candidates.end(),
                    [](const auto &v) { return v.reason == "fixed_overlap"; }));
  auto overlapping = d;
  overlapping.node_x[1] = 25;
  check(search(overlapping).merges.size() ==
        2); // broad phase includes overlapping equal macros
}
void group_shapes_and_independence() {
  auto d = fixture();
  d.num_movable = 5;
  d.node_x = {20, 40, 60, 80, 100};
  d.node_y.assign(5, 20);
  d.node_size_x.assign(5, 10);
  d.node_size_y.assign(5, 10);
  d.node_macro.assign(5, 1);
  d.node_halo_x.assign(5, 2);
  d.node_halo_y.assign(5, 2);
  d.node_align_class.assign(5, 0);
  d.build_maps();
  const auto first = search(d);
  check(first.merges.size() == 2);
  std::vector<int> used;
  for (const auto &m : first.merges)
    used.insert(used.end(), m.members.begin(), m.members.end());
  std::sort(used.begin(), used.end());
  check(std::adjacent_find(used.begin(), used.end()) == used.end());
  check(std::any_of(first.candidates.begin(), first.candidates.end(),
                    [](const auto &c) { return c.reason == "shared_member"; }));
  RigidMacroMap map;
  check(merge_macros(d, map, first.merges).applied);
  const auto second = search(d);
  check(second.merges.size() == 1);
  for (int i : second.merges[0].members)
    check(d.node_size_x[i] ==
          24); // the lone 10-wide macro cannot join a 2x1 group
  check(merge_macros(d, map, second.merges).applied);
  check(map.groups.size() == 1 && map.groups[0].members.size() == 4);
}
void offset_gap_regions() {
  auto d = fixture();
  d.node_y[1] += 4;
  // The cell is in the slanted corridor near the first face, below the shared
  // projection rectangle. The complete facing-edge corridor must count it.
  const auto found = search(d);
  const auto pair = std::find_if(
      found.candidates.begin(), found.candidates.end(),
      [](const auto &c) { return c.a == 0 && c.b == 1 && c.axis == Axis::x; });
  check(pair != found.candidates.end() && pair->selected);
  check(pair->gap_region ==
        (MagnetRegion{{{32, 18}, {38, 22}, {38, 36}, {32, 32}}}));
  check(std::abs(pair->gap_density - 4. / 84) < 1e-12);
  const MagnetDensity fine(d, 200, 200);
  check(std::abs(fine.region_density(pair->gap_region) - pair->gap_density) <
        1e-12);
  auto reverse = pair->gap_region;
  std::reverse(reverse.begin(), reverse.end());
  check(std::abs(fine.region_density(reverse) - pair->gap_density) < 1e-12);
  MagnetOptions o;
  o.density_fraction = 1. / 16;
  const auto rejected = search(d, o);
  check(std::any_of(rejected.candidates.begin(), rejected.candidates.end(),
                    [](const auto &c) {
                      return c.a == 0 && c.b == 1 && c.reason == "gap_density";
                    }));
  // A cell inside the corridor's bounding box but outside its slanted edge
  // contributes no area. Both axes use the same polygon clipping semantics.
  d.node_x[4] = 32;
  d.node_y[4] = 35;
  d.node_size_x[4] = d.node_size_y[4] = 1;
  for (int axis = 0; axis < 2; ++axis) {
    const auto empty = search(d);
    check(std::any_of(
        empty.candidates.begin(), empty.candidates.end(), [](const auto &c) {
          return c.a == 0 && c.b == 1 && c.selected && c.gap_density == 0;
        }));
    const MagnetDensity grid(d, 200, 200);
    for (const auto &c : empty.candidates)
      if (c.a == 0 && c.b == 1)
        check(grid.region_density(c.gap_region) == 0);
    std::swap(d.node_x, d.node_y);
    std::swap(d.node_size_x, d.node_size_y);
  }
  // Constant physical density is invariant under offset and partial-bin/core
  // clipping, even on a coarse, non-square grid.
  d = {};
  d.layout = {0, 0, 200, 200, 1, 1};
  d.num_movable = 1;
  d.node_x = d.node_y = {0};
  d.node_size_x = d.node_size_y = {200};
  d.build_maps();
  const MagnetDensity uniform(d, 7, 13);
  for (const auto &region :
       {pair->gap_region,
        MagnetRegion{{{-10, 0}, {40, 20}, {40, 80}, {-10, 60}}},
        MagnetRegion{{{20, 20}, {20, 20}, {20, 40}, {20, 40}}}}) {
    const bool degenerate = region[0] == region[1];
    check(std::abs(uniform.region_density(region) - (degenerate ? 0. : 1.)) <
          1e-12);
  }
}
void boundary_candidate_translation() {
  Design d;
  d.layout = {0, 0, 200, 200, 1, 1};
  d.num_movable = 2;
  d.node_x = {20, 20};
  d.node_y = {2, 15};
  d.node_size_x = d.node_size_y = {10, 10};
  d.node_halo_x = d.node_halo_y = {2, 2};
  d.node_macro = {1, 1};
  d.node_align_class = {1, 7};
  d.build_maps();
  const auto found = search(d);
  check(found.merges.size() == 1);
  const auto selected =
      std::find_if(found.candidates.begin(), found.candidates.end(),
                   [](const auto &c) { return c.selected; });
  check(selected != found.candidates.end() && selected->axis == Axis::y);
  check(selected->center_shift_x == 0 && selected->center_shift_y == .5);
  check(selected->merged_box.yl == 0 && selected->merged_box.yh == 28);
  check(selected->merge.poses ==
        (std::vector<std::pair<double, double>>{{20, 2}, {20, 16}}));
  auto grouped = d;
  RigidMacroMap map;
  check(merge_macros(grouped, map, found.merges).applied);
  check(macros_are_legal(expand_macros(grouped, d, map)));
  d.layout.yh = 26; // the padded union cannot fit: translation cannot help
  const auto rejected = search(d);
  check(rejected.merges.empty());
  check(std::any_of(rejected.candidates.begin(), rejected.candidates.end(),
                    [](const auto &c) {
                      return c.axis == Axis::y && c.reason == "outside_core";
                    }));
}
void other_macros_block_the_gap_region() {
  auto d = fixture();
  d.node_y[1] += 4;
  // This differently sized movable macro is inside the slanted corridor and
  // outside the shared-projection rectangle. Density excludes macros, so it
  // needs its own intersection check independent of shape and hierarchy.
  d.node_macro[4] = 1;
  auto found = search(d);
  check(found.merges.size() == 1);
  check(std::any_of(found.candidates.begin(), found.candidates.end(),
                    [](const auto &c) {
                      return c.a == 0 && c.b == 1 && c.reason == "macro_in_gap";
                    }));
  // Bounding-box corners outside the actual quadrilateral do not obstruct it.
  d.node_x[4] = 32;
  d.node_y[4] = 35;
  d.node_size_x[4] = d.node_size_y[4] = 1;
  found = search(d);
  check(found.merges.size() == 2);
  check(std::any_of(
      found.candidates.begin(), found.candidates.end(),
      [](const auto &c) { return c.a == 0 && c.b == 1 && c.selected; }));
}
void movement_and_events() {
  auto d = fixture();
  MagnetOptions options;
  options.settle_window = 2;
  options.settle_displacement = .02;
  MagnetStability gate(options);
  for (int i = 0; i <= 2; ++i) {
    const auto movement = gate.observe(d, d.node_x, d.node_y, i);
    check(gate.gate(i) == (i == 2));
    if (i == 2)
      check(movement.window_complete && gate.armed());
  }
  int iteration = 2;
  for (auto event :
       {MagnetEvent::timing, MagnetEvent::pull, MagnetEvent::terminal,
        MagnetEvent::snap, MagnetEvent::round}) {
    gate.event(iteration, event);
    check(!gate.settled(iteration));
    for (int j = 0; j <= 2; ++j) {
      gate.observe(d, d.node_x, d.node_y, iteration + j);
      check(gate.gate(iteration + j) == (j == 2));
    }
    iteration += 2;
  }
  // Returning to the same position cannot hide motion inside the window.
  d.node_x[0] += 1;
  auto movement = gate.observe(d, d.node_x, d.node_y, ++iteration);
  check(std::abs(movement.maximum - .1) < 1e-12);
  check(std::abs(movement.rms - .05) < 1e-12);
  d.node_x[0] -= 1;
  movement = gate.observe(d, d.node_x, d.node_y, ++iteration);
  check(std::abs(movement.maximum - .1) < 1e-12 && !gate.gate(iteration));
}
void per_macro_motion_uses_own_extent() {
  auto d = fixture();
  d.node_size_x[0] = d.node_size_y[0] = 100;
  MagnetOptions options;
  options.settle_window = 2;
  MagnetStability stability(options);
  for (int i = 0; i <= 2; ++i) {
    d.node_x[0] = 20 + i;
    stability.observe(d, d.node_x, d.node_y, i);
  }
  check(std::abs(stability.macro_movement(0) - .02) < 1e-12);
  check(stability.macro_settled(0) && stability.macro_settled(1));
  d.node_x[1] += 1;
  stability.observe(d, d.node_x, d.node_y, 3);
  check(stability.macro_settled(0) && !stability.macro_settled(1));
}
void many_macros() {
  Design d;
  d.layout = {0, 0, 1000, 3000, 1, 1};
  d.num_movable = 1040;
  for (int i = 0; i < 1040; ++i) {
    d.node_x.push_back(20 + (i % 40) * 20);
    d.node_y.push_back(20 + (i / 40) * 60);
    d.node_size_x.push_back(10);
    d.node_size_y.push_back(10);
    d.node_halo_x.push_back(2);
    d.node_halo_y.push_back(2);
    d.node_macro.push_back(1);
  }
  d.build_maps();
  const auto found = search(d);
  check(!found.exact_density && found.merges.size() == 520);
  RigidMacroMap map;
  check(merge_macros(d, map, found.merges).applied && map.groups.size() == 520);
  std::cout << "1040-macro candidate search: " << found.seconds << " s, "
            << found.merges.size() << " independent pairs\n";
}
} // namespace
int main() {
  try {
    density_and_selection();
    group_shapes_and_independence();
    offset_gap_regions();
    boundary_candidate_translation();
    other_macros_block_the_gap_region();
    movement_and_events();
    per_macro_motion_uses_own_extent();
    many_macros();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
