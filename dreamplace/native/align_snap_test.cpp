// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("snap assertion");
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 80, 80, 1, 1};
    d.num_movable = 3;
    d.node_x = {10, 10.6, 10.2, 30};
    d.node_y = {2, 12, 24, 32};
    d.node_size_x.assign(4, 4);
    d.node_size_y.assign(4, 4);
    d.node_halo_x.assign(4, 1);
    d.node_halo_y.assign(4, 1);
    d.node_macro.assign(4, 1);
    d.node_align_class = {0, 0, 1, 0};
    d.build_maps();
    auto snap = detail::macro_align_snap(d, .25, 100);
    check(snap.groups.size() == 1 &&
          snap.groups[0].nodes == std::vector<int>({0, 1}));
    check(snap.constraints.size() == 1 && snap.constraints[0].a == 0 &&
          snap.constraints[0].b == 1);
    auto fixed_x = d.node_x[3], fixed_y = d.node_y[3];
    Options o;
    o.max_iterations = 0;
    o.macro_snap = true;
    o.detailed_place = false;
    auto r = place(d, o);
    check(is_legal(d) && d.node_x[0] == d.node_x[1]);
    check(r.alignment.exact_groups == 1 && r.alignment.soft_violations == 0 &&
          r.alignment.members == 2);
    check(d.node_x[3] == fixed_x && d.node_y[3] == fixed_y);
    // A one-row SRAM still belongs to macro legalization in the opt-in path.
    Design tiny;
    tiny.layout = {0, 0, 20, 10, 1, 1};
    tiny.num_movable = 2;
    tiny.node_x = {1, 1.4};
    tiny.node_y = {1, 4};
    tiny.node_size_x = {4, 4};
    tiny.node_size_y = {1, 1};
    tiny.node_macro = {1, 1};
    tiny.node_align_class = {0, 0};
    tiny.build_maps();
    auto tiny_result = place(tiny, o);
    check(is_legal(tiny) && tiny.node_x[0] == tiny.node_x[1] &&
          tiny_result.alignment.exact_groups == 1);
    // Snap is soft: a deliberately conflicting pair must remain separated.
    Design conflict;
    conflict.layout = {0, 0, 20, 8, 1, 1};
    conflict.num_movable = 2;
    conflict.node_x = {0, 6};
    conflict.node_y = {0, 0};
    conflict.node_size_x = {4, 4};
    conflict.node_size_y = {8, 8};
    conflict.node_macro = {1, 1};
    conflict.build_maps();
    AlignmentSnap forced;
    forced.groups = {{0, Axis::x, 2, {0, 1}}};
    forced.constraints = {{0, 1, 0, 100}};
    auto before = conflict;
    detail::legalize(conflict, 8, 8, &forced);
    auto report = detail::macro_alignment_report(before, conflict, forced);
    check(is_legal(conflict) && report.soft_violations == 1 &&
          report.broken_groups == 1);
    auto malformed = forced;
    malformed.constraints[0].b = 2;
    bool rejected = false;
    try {
      detail::legalize(conflict, 8, 8, &malformed);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected);
    // Regularity uses an independent hierarchy/size class, ignoring the
    // deliberately conflicting synthesis-color classes. Fixed macros stay out.
    Design regular;
    regular.layout = {0, 0, 100, 100, 1, 1};
    regular.num_movable = 4;
    regular.node_x = {10, 10.7, 10.2, 30, 50};
    regular.node_y = {5, 23, 41.2, 64, 70};
    regular.node_size_x.assign(5, 4);
    regular.node_size_y.assign(5, 4);
    regular.node_macro.assign(5, 1);
    regular.node_halo_x.assign(5, 1);
    regular.node_halo_y.assign(5, 1);
    regular.node_align_class = {0, 1, 2, 0, 0};
    regular.node_snap_class = {7, 7, 7, 8, 7};
    regular.build_maps();
    auto plain = regular, annotated = regular;
    plain.node_snap_class.clear();
    auto old_options = o;
    old_options.max_iterations = 40;
    old_options.macro_snap = false;
    auto old_result = place(plain, old_options);
    auto annotated_result = place(annotated, old_options);
    check(plain.node_x == annotated.node_x &&
          plain.node_y == annotated.node_y &&
          old_result.hpwl == annotated_result.hpwl);
    const auto scoped_alignment = detail::macro_align_snap(regular, .25, 100);
    check(!scoped_alignment.groups.empty());
    for (const auto &constraint : scoped_alignment.constraints)
      check(constraint.a < 3 && constraint.b < 3);
    o.macro_snap_regular = true;
    const auto regular_result = place(regular, o);
    check(is_legal(regular) && regular.node_x[0] == regular.node_x[1] &&
          regular.node_x[1] == regular.node_x[2]);
    check(regular.node_x[4] == 50 && regular.node_y[4] == 70);
    check(regular_result.alignment.regularity.feasible &&
          regular_result.alignment.regularity.circles.size() == 3);
    check(std::abs(regular_result.alignment.regularity.circles[0].radius -
                   std::sqrt(36. / std::acos(-1.))) < 1e-8);
    bool spacing = false;
    for (const auto &relation : regular_result.alignment.regularity.relations)
      spacing |= relation.kind == "spacing";
    check(spacing);
    for (const auto &constraint : regular_result.alignment.constraints)
      check(constraint.a < 3 && constraint.b < 3);
    // Lattice rounding is shared by soft relative offsets and anchor targets.
    for (const auto &constraint : regular_result.alignment.constraints) {
      const auto &pos = constraint.axis == 0 ? regular.node_x : regular.node_y;
      check(std::abs(pos[constraint.b] - pos[constraint.a] -
                     constraint.offset) < 1e-6);
    }
    // A nonzero soft offset that conflicts with hard separation must break.
    conflict = before;
    forced.groups.clear();
    forced.constraints = {{0, 1, 0, 100, 2}};
    detail::legalize(conflict, 8, 8, &forced);
    auto offset_report =
        detail::macro_alignment_report(before, conflict, forced);
    check(is_legal(conflict) && offset_report.soft_violations == 1);
    auto invalid = regular;
    invalid.node_size_x[1] += 1;
    bool bad_size = false;
    try {
      detail::macro_regular_snap(invalid, {}, 100);
    } catch (const std::invalid_argument &) {
      bad_size = true;
    }
    check(bad_size);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
