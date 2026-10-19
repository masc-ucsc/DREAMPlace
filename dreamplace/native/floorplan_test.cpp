// BSD-3-Clause. See LICENSE. Assertions stay enabled in optimized builds.
#include "dreamplace/native/floorplan.hpp"
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool b, const char *s) {
  if (!b)
    throw std::runtime_error(s);
}
void near(double a, double b, double tol, const char *s) {
  check(std::abs(a - b) < tol, s);
}
FloorplanCircle circle(double x, double y, double r = 8, int family = 0,
                       int parent = -1) {
  return {"c", x, y, r, r * r, family, parent};
}
void legal(const FloorplanSnap &r) {
  check(r.feasible && r.maximum_overlap < 1e-6,
        "snap returned overlapping targets");
  for (size_t i = 0; i < r.circles.size(); ++i)
    for (size_t j = 0; j < i; ++j) {
      const auto &a = r.circles[i], &b = r.circles[j];
      if (a.family == b.family)
        check(std::hypot(a.x - b.x, a.y - b.y) + 1e-5 >= a.radius + b.radius,
              "circle reach used one radius");
    }
}
int main() {
  std::vector<FloorplanCircle> mesh;
  for (int y = 0; y < 3; ++y)
    for (int x = 0; x < 3; ++x)
      mesh.push_back(
          circle(30 + x * 25 + (y - 1) * .7, 30 + y * 25 + (x - 1) * .6));
  auto m = snap_floorplan(mesh, 0, 0, 120, 120);
  legal(m);
  for (int row = 0; row < 3; ++row) {
    near(m.circles[3 * row].y, m.circles[3 * row + 2].y, 1e-4,
         "mesh row not aligned");
    near(m.circles[3 * row + 1].x * 2,
         m.circles[3 * row].x + m.circles[3 * row + 2].x, 1e-4,
         "mesh gaps differ");
  }
  for (int col = 0; col < 3; ++col)
    near(m.circles[col].x, m.circles[6 + col].x, 1e-4,
         "mesh column not aligned");
  // Three rows with 1,2,3 members. The middle row is half-pitch staggered.
  std::vector<FloorplanCircle> triangle = {
      circle(60.7, 25),   circle(47, 50.4),   circle(72.5, 49.7),
      circle(34.9, 75.5), circle(60.3, 74.6), circle(85.5, 75.2)};
  auto t = snap_floorplan(triangle, 0, 0, 120, 120);
  legal(t);
  near(t.circles[1].y, t.circles[2].y, 1e-4, "triangle two-member row lost");
  near(t.circles[0].x, (t.circles[1].x + t.circles[2].x) / 2, 1e-4,
       "triangle symmetry axis not shared");
  near(t.circles[0].x, t.circles[4].x, 1e-4, "triangle central column lost");
  near(t.circles[2].x - t.circles[1].x, t.circles[4].x - t.circles[3].x, 1e-4,
       "triangle spacing not shared");
  // A root, two children and four leaves; two distinct pitches remain useful.
  std::vector<FloorplanCircle> tree = {circle(70.5, 20),   circle(42.7, 50.5),
                                       circle(98, 49.6),   circle(33.4, 80.1),
                                       circle(58.9, 79.5), circle(83.7, 80.4),
                                       circle(109.1, 80)};
  auto tr = snap_floorplan(tree, 0, 0, 145, 110);
  legal(tr);
  near(tr.circles[0].x, (tr.circles[1].x + tr.circles[2].x) / 2, 1e-4,
       "tree children not centered on root");
  near(tr.circles[0].x, (tr.circles[3].x + tr.circles[6].x) / 2, 1e-4,
       "tree leaves not centered");
  // Projection overlap below 50% does not select a common vertical line.
  auto stagger =
      snap_floorplan({circle(30, 30), circle(39, 60)}, 0, 0, 100, 100);
  legal(stagger);
  near(stagger.circles[0].x, 30, 1e-6, "less-than-half projection snapped");
  near(stagger.circles[1].x, 39, 1e-6, "less-than-half projection snapped");
  // Preserve two independently translated families, rather than one global
  // grid.
  auto families = snap_floorplan({circle(20, 20, 5, 0), circle(40, 20.2, 5, 0),
                                  circle(67, 70, 5, 1), circle(89, 70.3, 5, 1)},
                                 0, 0, 120, 100);
  legal(families);
  near(families.circles[0].x, 20, 1e-6, "unrelated family changed line origin");
  near(families.circles[2].x, 67, 1e-6, "unrelated family changed line origin");
  // Exact regular input is an identity; repeatability and affine unit
  // invariance.
  auto perfect = snap_floorplan(m.circles, 0, 0, 120, 120);
  legal(perfect);
  for (size_t i = 0; i < mesh.size(); ++i) {
    near(perfect.circles[i].x, m.circles[i].x, 1e-4, "perfect mesh moved");
    near(perfect.circles[i].y, m.circles[i].y, 1e-4, "perfect mesh moved");
  }
  auto scaled = mesh;
  for (auto &c : scaled) {
    c.x = 1000 * c.x - 300;
    c.y = 1000 * c.y + 500;
    c.radius *= 1000;
    c.area *= 1e6;
  }
  auto scale = snap_floorplan(scaled, -300, 500, 119700, 120500);
  legal(scale);
  for (size_t i = 0; i < mesh.size(); ++i)
    near(scale.circles[i].x, m.circles[i].x * 1000 - 300, .02,
         "DBU scaling changed snap");
  // Bounded-span grouping must not collapse a chain of individually close x's.
  auto chain = snap_floorplan(
      {circle(20, 20, 10), circle(26, 50, 10), circle(32, 80, 10)}, 0, 0, 100,
      120);
  legal(chain);
  check(std::abs(chain.circles[0].x - chain.circles[2].x) > 1,
        "single-link chain collapsed");
  // An impossible family is returned unchanged with an explicit failure.
  FloorplanOptions none;
  none.max_displacement = 0;
  std::vector<FloorplanCircle> impossible = {circle(30, 30, 20),
                                             circle(35, 35, 20)};
  auto bad = snap_floorplan(impossible, 0, 0, 100, 100, none);
  check(!bad.feasible && bad.infeasible_families == 1,
        "infeasible hard circles silently accepted");
  near(bad.circles[0].x, 30, 1e-12, "failed fit changed inputs");
  // Parent target translations carry descendants before their own type is fit.
  std::vector<FloorplanCircle> nested = {
      circle(50, 40, 20, 0), circle(55, 120, 20, 0), circle(50, 55, 5, 1, 0),
      circle(55, 135, 5, 1, 1)};
  auto nested_snap = snap_floorplan(nested, 0, 0, 160, 160);
  check(nested_snap.feasible, "nested target snap failed");
  near(nested_snap.circles[2].x - nested[2].x,
       nested_snap.circles[0].x - nested[0].x, 1e-5,
       "child did not inherit parent target translation");
  near(nested_snap.circles[3].x - nested[3].x,
       nested_snap.circles[1].x - nested[1].x, 1e-5,
       "second child did not inherit parent translation");
  // Repeated child types in different candidates must not share snap lines,
  // gap classes or symmetry axes, even when their coordinates nearly match.
  std::vector<FloorplanCircle> siblings = {circle(50, 40, 20, 0),
                                           circle(50, 140, 20, 0)};
  for (int parent = 0; parent < 2; ++parent)
    for (int row = 0; row < 2; ++row)
      for (int col = 0; col < 3; ++col)
        siblings.push_back(circle(30 + parent * 2 + col * (20 + parent * 2),
                                  30 + parent * 100 + row * 20, 5, 1, parent));
  auto scoped = snap_floorplan(siblings, 0, 0, 180, 180);
  legal(scoped);
  int child_relations = 0;
  for (const auto &relation : scoped.relations) {
    const auto &first = siblings[relation.nodes.front()];
    for (int n : relation.nodes)
      check(siblings[n].family == first.family &&
                siblings[n].parent == first.parent,
            "snap relation crossed candidate parents");
    child_relations += first.parent >= 0;
  }
  check(child_relations > 0, "sibling snap was disabled entirely");
  near(scoped.circles[2].x, 30, 1e-5, "first candidate line changed");
  near(scoped.circles[8].x, 32, 1e-5, "second candidate copied another line");
  near(scoped.circles[3].x - scoped.circles[2].x, 20, 1e-5,
       "first candidate gap changed");
  near(scoped.circles[9].x - scoped.circles[8].x, 22, 1e-5,
       "gap class crossed parents");
  // Unequal radii also require the sum of both radii.
  auto unequal =
      snap_floorplan({circle(30, 30, 10), circle(55, 31, 20)}, 0, 0, 120, 120);
  legal(unequal);
  // FP boundaries reserve the entire circle, including corner cases and
  // children translated by their parent's snap. Macro center-only mode remains
  // opt-in unchanged through the default FloorplanOptions.
  FloorplanOptions bounded;
  bounded.contain_circles = true;
  const auto contained = [&](const FloorplanSnap &result) {
    check(result.feasible && result.maximum_boundary_violation == 0,
          "circle containment not feasible");
    for (const auto &c : result.circles)
      check(c.x - c.radius >= 0 && c.x + c.radius <= 100 &&
                c.y - c.radius >= 0 && c.y + c.radius <= 100,
            "centroid circle crosses core edge");
  };
  auto corners = snap_floorplan({circle(2, 3, 10), circle(98, 97, 10),
                                 circle(2, 97, 10), circle(98, 3, 10)},
                                0, 0, 100, 100, bounded);
  contained(corners);
  legal(corners);
  auto edge_child = snap_floorplan(
      {circle(2, 30, 10, 0), circle(5, 70, 10, 0), circle(90, 40, 8, 1, 0)}, 0,
      0, 100, 100, bounded);
  contained(edge_child);
  bounded.max_displacement = 0;
  contained(snap_floorplan({circle(0, 0, 10)}, 0, 0, 100, 100, bounded));
  auto oversized =
      snap_floorplan({circle(50, 50, 60)}, 0, 0, 100, 100, bounded);
  check(!oversized.feasible && oversized.maximum_boundary_violation > 0,
        "oversized circle silently shrunk to fit");
  // Frozen targets have a real translation force and a correct analytic
  // gradient.
  Design d;
  d.layout = {0, 0, 100, 100, 1, 1};
  d.num_movable = 2;
  d.node_x = {20, 40};
  d.node_y = {30, 50};
  d.node_size_x = d.node_size_y = {2, 2};
  d.node_centroid_color = {-1, -1};
  d.centroid_node_area = {4, 4};
  d.centroid_groups = {{"g", -1, 0, 0, 0, 0, {0, 1}, {}}};
  d.centroid_targets = {{60, 70}};
  d.build_maps();
  d.validate();
  const auto g = detail::centroid_attraction(d, 1);
  check(g.x[0] + g.x[1] < 0, "frozen target has no translation force");
  const double eps = 1e-4;
  for (int i = 0; i < 2; ++i) {
    d.node_x[i] += eps;
    const double plus = detail::centroid_attraction(d, 1).value;
    d.node_x[i] -= 2 * eps;
    const double minus = detail::centroid_attraction(d, 1).value;
    d.node_x[i] += eps;
    near((plus - minus) / (2 * eps), g.x[i], 1e-5, "frozen target gradient");
  }
  d.centroid_targets.clear();
  const auto moving = detail::centroid_attraction(d, 1);
  near(moving.x[0] + moving.x[1], 0, 1e-12, "moving target default changed");
  bool rejected = false;
  try {
    auto invalid = triangle;
    invalid[0].radius = std::numeric_limits<double>::quiet_NaN();
    snap_floorplan(invalid, 0, 0, 120, 120);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected, "invalid circle accepted");
  std::cout << "mesh, triangle, tree, stagger, families, chains, infeasible "
               "circles, frozen target gradient pass\n";
}
