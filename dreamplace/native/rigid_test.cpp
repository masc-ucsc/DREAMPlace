// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rigid.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("rigid assertion at line " +
                             std::to_string(at.line()));
}
template <class F> void rejects(F f) {
  try {
    f();
  } catch (const std::invalid_argument &) {
    return;
  }
  throw std::runtime_error("expected invalid rigid operation");
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 200, 200, 1, 1};
  d.num_movable = 5;
  d.node_x = {20, 34, 20, 34, 32, 150, 0};
  d.node_y = {20, 20, 38, 38, 20, 150, 0};
  d.node_size_x = {10, 10, 10, 10, 1, 10, 0};
  d.node_size_y = {12, 12, 12, 12, 1, 12, 0};
  d.node_macro = {1, 1, 1, 1, 0, 1, 0};
  d.node_halo_x = {2, 2, 2, 2, 0, 2, 0};
  d.node_halo_y = {3, 3, 3, 3, 0, 3, 0};
  d.node_orientation = {Orientation::FN, Orientation::FS, Orientation::N,
                        Orientation::S,  Orientation::N,  Orientation::FN,
                        Orientation::N};
  d.node_align_class = {0, 1, 0, 1, -1, 0, -1};
  d.node_snap_class = {0, 1, 0, 1, -1, 0, -1};
  d.pin2node = {0, 1, 5, 2, 3, 4, 0, 1};
  d.pin_offset_x = {1, 2, 3, 4, 5, .5, 6, 7};
  d.pin_offset_y = {2, 3, 4, 5, 6, .5, 7, 8};
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 3, 6, 8};
  d.build_maps();
  d.validate();
  return d;
}
MacroMerge request(const Design &d, std::vector<int> members) {
  MacroMerge r;
  r.members = std::move(members);
  for (int i : r.members)
    r.poses.emplace_back(d.node_x[i], d.node_y[i]);
  return r;
}
void pin_positions(const Design &a, const Design &b) {
  check(a.pin2node.size() == b.pin2node.size());
  for (size_t p = 0; p < a.pin2node.size(); ++p) {
    check(std::abs(a.node_x[a.pin2node[p]] + a.pin_offset_x[p] -
                   b.node_x[b.pin2node[p]] - b.pin_offset_x[p]) < 1e-9);
    check(std::abs(a.node_y[a.pin2node[p]] + a.pin_offset_y[p] -
                   b.node_y[b.pin2node[p]] - b.pin_offset_y[p]) < 1e-9);
  }
}
void geometry_and_expansion() {
  for (bool vertical : {false, true})
    for (bool reverse : {false, true}) {
      const auto source = fixture();
      auto d = source;
      RigidMacroMap map;
      std::vector<int> members{0, vertical ? 2 : 1};
      if (reverse)
        std::reverse(members.begin(), members.end());
      const auto r = request(d, members);
      const auto result = merge_macros(d, map, std::vector{r});
      check(result.applied && result.reason.empty());
      check(d.node_x.size() == 6 && d.num_movable == 4 &&
            map.groups.size() == 1);
      check(d.node_size_x[0] == (vertical ? 10 : 24));
      check(d.node_size_y[0] == (vertical ? 30 : 12));
      check(d.node_halo_x[0] == 2 && d.node_halo_y[0] == 3);
      check(d.node_orientation_fixed[0] &&
            d.node_group_id[0] == map.groups[0].id);
      check(d.node_align_class[0] > 1 && d.node_snap_class[0] > 1);
      check(result.previous_node[0] == -1);
      pin_positions(source, d);
      const auto expanded = expand_macros(d, source, map);
      check(expanded.node_x == source.node_x &&
            expanded.node_y == source.node_y);
      check(expanded.node_size_x == source.node_size_x &&
            expanded.node_size_y == source.node_size_y);
      check(expanded.node_orientation == source.node_orientation);
      check(expanded.pin_offset_x == source.pin_offset_x &&
            expanded.pin_offset_y == source.pin_offset_y);
      check(expanded.net_mask == source.net_mask);
      check(expanded.node_group_id[members[0]] ==
            expanded.node_group_id[members[1]]);
      // An exact halo-tiling pair has the same density charge and overflow on
      // several grids, both with and without scaled macro occupancy.
      for (bool scaled : {false, true})
        for (int bins : {4, 8, 16}) {
          auto before = source, after = d;
          Options o;
          o.max_iterations = 0;
          o.bins_x = o.bins_y = bins;
          o.target_density = .1;
          o.add_fillers = o.legalize = o.detailed_place = false;
          o.scale_macro_density = scaled;
          const auto a = place(before, o), b = place(after, o);
          check(std::abs(a.overflow - b.overflow) < 1e-12);
        }
      const auto x = d.pin_offset_x, y = d.pin_offset_y;
      rejects([&] { set_orientation(d, 0, Orientation::FN); });
      check(d.node_orientation[0] == Orientation::N && d.pin_offset_x == x &&
            d.pin_offset_y == y);
      set_orientation(d, 0, Orientation::N);
      d.node_orientation[0] =
          Orientation::FS; // reject bypassing the public API, too
      rejects([&] { (void)expand_macros(d, source, map); });
    }
}
void atomic_batch_and_composition() {
  const auto source = fixture();
  auto d = source;
  RigidMacroMap map;
  const auto batch = std::vector{request(d, {0, 1}), request(d, {2, 3})};
  check(merge_macros(d, map, batch).applied);
  check(map.groups.size() == 2 &&
        d.node_align_class[0] == d.node_align_class[1]);
  check(d.node_snap_class[0] == d.node_snap_class[1]);
  const auto merge = request(d, {1, 0});
  check(merge_macros(d, map, std::vector{merge}).applied);
  check(map.groups.size() == 1 && map.groups[0].members.size() == 4);
  check(d.node_size_x[0] == 24 && d.node_size_y[0] == 30);
  d.node_x[0] += 5;
  d.node_y[0] += 7;
  auto expanded = expand_macros(d, source, map);
  for (int i = 0; i < 4; ++i) {
    check(expanded.node_x[i] == source.node_x[i] + 5 &&
          expanded.node_y[i] == source.node_y[i] + 7);
    check(expanded.node_orientation[i] == source.node_orientation[i]);
    check(expanded.node_group_id[i] == map.groups[0].id);
  }
  check(expanded.node_x[4] == source.node_x[4] &&
        expanded.node_y[4] == source.node_y[4]);
  check(expanded.node_x[5] == source.node_x[5] &&
        expanded.node_y[5] == source.node_y[5]);
  pin_positions(d, expanded);
  for (int i = 0; i < 4; ++i) {
    expanded.node_x[i] += 3;
    expanded.node_y[i] -= 2;
  }
  const auto contracted = contract_macros(expanded, d, source, map);
  check(contracted.node_x[0] == d.node_x[0] + 3);
  check(contracted.node_y[0] == d.node_y[0] - 2);
  expanded.node_x[0] += 1;
  rejects([&] { (void)contract_macros(expanded, d, source, map); });
}
void atomic_rejection() {
  for (int mode = 0; mode < 6; ++mode) {
    auto d = fixture();
    RigidMacroMap map;
    auto a = request(d, {0, 1}), b = request(d, {2, 3});
    if (mode == 0)
      b.poses[1].first += 1; // unfilled outer rectangle
    if (mode == 1)
      b.poses = {{-2, 38}, {12, 38}}; // outside core including halo
    if (mode == 2)
      b.poses = {{140, 150}, {154, 150}}; // fixed obstacle
    if (mode == 3)
      b.members[0] = 5; // fixed member
    if (mode == 4)
      b.members[0] = 0; // shared member
    if (mode == 5)
      b.poses = a.poses; // independent members, overlapping merged rectangles
    const auto source = d;
    const auto result = merge_macros(d, map, std::vector{a, b});
    check(!result.applied && !result.reason.empty() &&
          result.group_ids.empty());
    check(d.node_x == source.node_x && d.node_y == source.node_y &&
          d.pin2node == source.pin2node &&
          d.pin_offset_x == source.pin_offset_x &&
          d.node_size_x == source.node_size_x &&
          d.node_orientation == source.node_orientation &&
          d.flat_node2pin == source.flat_node2pin &&
          d.net_mask == source.net_mask);
    check(map.source_to_node.empty() && map.groups.empty() &&
          map.next_id == 1 && map.shape_classes.empty());
  }
}
void projected_cohesion_and_terminal_updates() {
  auto source = fixture();
  source.cohesion_groups = {{{0, 2, 4}, 2, {}}, {{1, 3}, .7, {}}};
  source.validate();
  auto d = source;
  RigidMacroMap map;
  check(merge_macros(d, map, std::vector{request(d, {0, 1})}).applied);
  const auto before = detail::cohesion(source, .3),
             after = detail::cohesion(d, .3);
  check(std::abs(before.value - after.value) < 1e-9);
  std::vector<double> x(d.node_x.size()), y(d.node_x.size());
  for (size_t i = 0; i < source.node_x.size(); ++i) {
    x[map.source_to_node[i]] += before.x[i];
    y[map.source_to_node[i]] += before.y[i];
  }
  for (size_t i = 0; i < x.size(); ++i) {
    check(std::abs(x[i] - after.x[i]) < 1e-9 &&
          std::abs(y[i] - after.y[i]) < 1e-9);
    auto hi = d, lo = d;
    hi.node_x[i] += .001;
    lo.node_x[i] -= .001;
    const double fd =
        (detail::cohesion(hi, .3).value - detail::cohesion(lo, .3).value) /
        .002;
    check(std::abs(fd - after.x[i]) < 1e-7);
  }
  const int terminal = map.source_to_node.back();
  d.node_x[terminal] = 200;
  d.node_y[terminal] = 50;
  const auto expanded = expand_macros(d, source, map);
  check(expanded.node_x.back() == 200 && expanded.node_y.back() == 50);
  check(expanded.cohesion_groups[0].offsets.empty());
}
void projected_hierarchy_moments_and_forces() {
  auto source = fixture();
  source.centroid_color_names = {"left", "right"};
  source.node_centroid_color = {0, 1, 0, 1, -1, -1, -1};
  source.centroid_node_area = {2, 7, 3, 11, 1, 13, 0};
  CentroidGroup root, left, right;
  root.name = "root";
  root.nodes = {4};
  root.colors = {0, 1};
  left.name = "root.left";
  left.parent = 0;
  left.colors = {0};
  right.name = "root.right";
  right.parent = 0;
  right.colors = {1};
  source.centroid_groups = {root, left, right};
  source.validate();
  auto grouped = source;
  RigidMacroMap map;
  for (int level = 0; level < 3; ++level) {
    const auto members =
        level == 1 ? std::vector{map.source_to_node[2], map.source_to_node[3]}
                   : std::vector{0, 1};
    check(merge_macros(grouped, map, std::vector{request(grouped, members)})
              .applied);
    const auto before_moments = detail::centroid_summary(source);
    const auto after_moments = detail::centroid_summary(grouped);
    for (size_t k = 0; k < before_moments.size(); ++k) {
      const auto &a = before_moments[k], &b = after_moments[k];
      check(std::abs(a.x - b.x) < 1e-12 && std::abs(a.y - b.y) < 1e-12 &&
            std::abs(a.radius - b.radius) < 1e-12 && a.area == b.area &&
            a.cells == b.cells && a.macros == b.macros);
    }
    for (int mode = 0; mode < 7; ++mode) {
      auto original = source, contracted = grouped;
      if (mode == 2)
        for (auto *d : {&original, &contracted})
          for (auto &g : d->centroid_groups)
            g.repulsion_area = 10000;
      if (mode >= 5) {
        original.centroid_targets = {{60, 60}, {75, 50}, {30, 80}};
        contracted.centroid_targets = original.centroid_targets;
      }
      const auto energy = [&](const Design &d) {
        switch (mode) {
        case 0:
        case 5:
          return detail::centroid_attraction(d, .3);
        case 1:
          return detail::centroid_circle_attraction(d, .3, .5, 1.2, false);
        case 2:
          return detail::centroid_repulsion(d, .3, .5);
        case 3:
          return detail::centroid_parent_attraction(d, .3, .5, 1.5, nullptr,
                                                    nullptr);
        case 4:
          return detail::centroid_alignment(d, .3, .2, 4, .5, nullptr, nullptr);
        default:
          return detail::centroid_circle_attraction(d, .3, .5, 1.2, true);
        }
      };
      const auto before = energy(original), after = energy(contracted);
      check(std::abs(before.value - after.value) < 1e-9);
      std::vector<double> x(grouped.node_x.size()), y(x.size());
      for (size_t i = 0; i < source.node_x.size(); ++i) {
        x[map.source_to_node[i]] += before.x[i];
        y[map.source_to_node[i]] += before.y[i];
      }
      for (size_t i = 0; i < x.size(); ++i) {
        check(std::abs(x[i] - after.x[i]) < 1e-9 &&
              std::abs(y[i] - after.y[i]) < 1e-9);
        auto hi = contracted, lo = contracted;
        hi.node_x[i] += .001;
        lo.node_x[i] -= .001;
        const double fd = (energy(hi).value - energy(lo).value) / .002;
        check(std::abs(fd - after.x[i]) < 1e-6);
      }
    }
    detail::Gradient wire, density;
    wire.x = source.centroid_node_area;
    wire.y = source.centroid_node_area;
    density.x.assign(wire.x.size(), 0);
    density.y = density.x;
    const auto before =
        detail::centroid_force_reference(source, wire, density, 0);
    detail::Gradient rigid_wire, rigid_density;
    rigid_wire.x.assign(grouped.node_x.size(), 0);
    rigid_wire.y = rigid_wire.x;
    rigid_density.x = rigid_wire.x;
    rigid_density.y = rigid_wire.x;
    for (size_t i = 0; i < wire.x.size(); ++i) {
      rigid_wire.x[map.source_to_node[i]] += wire.x[i];
      rigid_wire.y[map.source_to_node[i]] += wire.y[i];
    }
    const auto after =
        detail::centroid_force_reference(grouped, rigid_wire, rigid_density, 0);
    check(std::abs(before.first - after.first) < 1e-12 &&
          std::abs(before.second - after.second) < 1e-12);
  }
  auto bad = grouped;
  bad.centroid_groups[0].areas.pop_back();
  rejects([&] { bad.validate(); });
  bad = grouped;
  bad.centroid_groups[0].offsets[0].first = INFINITY;
  rejects([&] { bad.validate(); });
}
} // namespace
int main() {
  try {
    geometry_and_expansion();
    atomic_batch_and_composition();
    atomic_rejection();
    projected_cohesion_and_terminal_updates();
    projected_hierarchy_moments_and_forces();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cout << "rigid macro geometry, pins, density, rejection and expansion "
               "passed\n";
}
