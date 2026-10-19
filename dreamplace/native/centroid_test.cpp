// BSD-3-Clause. See LICENSE. Checks remain active under -c opt.
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool b, const char *s) {
  if (!b)
    throw std::runtime_error(s);
}
void near(double a, double b, double tol, const char *s) {
  check(std::abs(a - b) <= tol, s);
}
template <class F> void rejects(F f) {
  bool bad = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    bad = true;
  }
  check(bad, "invalid centroid metadata accepted");
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 1000, 1000, 1, 1};
  d.num_movable = 8;
  d.node_x = {100, 160, 220, 600, 660, 720, 300, 900, 950};
  d.node_y = {200, 230, 260, 500, 530, 560, 400, 100, 950};
  d.node_size_x = d.node_size_y = {30, 30, 30, 30, 30, 30, 1, 1, 1};
  d.node_macro = {1, 1, 1, 1, 1, 1, 0, 0, 0};
  d.centroid_node_area = {900, 900, 900, 900, 900, 900, 1, 1, 1};
  d.node_centroid_color = {0, 1, 2, 0, 3, 4, 0, 0, -1};
  d.centroid_color_names = {"shared", "a1", "a2", "b1", "b2"};
  d.pin2node = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  d.pin_offset_x.assign(9, 1);
  d.pin_offset_y.assign(9, 1);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  d.flat_net2pin_start = {0, 4, 9};
  d.build_maps();
  std::vector<std::string> paths = {"top.a.m0",    "top.a.m1",      "top.a.m2",
                                    "top.b.m0",    "top.b.m1",      "top.b.m2",
                                    "top.a.logic", "outside.logic", "fixed"};
  std::vector<std::string> masters = {"ram", "ram", "ram", "ram", "ram",
                                      "ram", "and", "and", "and"};
  d.centroid_groups = detail::centroid_candidates(d, paths, masters, 2000, 2);
  d.validate();
  check(d.centroid_groups.size() == 2,
        "single-master repeated candidate missing");
  near(d.centroid_groups[0].equivalent_cells, 2701, 0,
       "equivalent area ignores logic/macros");
  check(detail::centroid_candidates(d, paths, masters, 2700, 2).empty(),
        "threshold is not strict or replication omitted");
  auto mismatch = masters;
  mismatch[5] = "other";
  check(detail::centroid_candidates(d, paths, mismatch, 2000, 2).empty(),
        "nonreplicated candidate accepted");
  auto small = paths;
  small[2] = "other.m2";
  small[5] = "else.m2";
  check(detail::centroid_candidates(d, small, masters, 0, 2).empty(),
        "two-macro candidate accepted");
  auto stats = detail::centroid_summary(d);
  check(stats[0].cells == 6 && stats[0].macros == 4,
        "outside shared-color cells/macros omitted from centroid");
  return d;
}
void candidate_copy_checks() {
  auto discover = [](const std::vector<int> &copies, int minimum = 4,
                     bool wrapper = false) {
    Design d;
    d.layout = {0, 0, 1000, 1000, 1, 1};
    std::vector<std::string> paths, masters;
    for (size_t parent = 0; parent < copies.size(); ++parent)
      for (int copy = 0; copy < copies[parent]; ++copy)
        for (int macro = 0; macro < 3; ++macro) {
          paths.push_back("top.p" + std::to_string(parent) + ".u" +
                          std::to_string(copy) + (wrapper ? ".wrapper" : "") +
                          ".m" + std::to_string(macro));
          masters.push_back("ram");
          d.node_x.push_back(0);
          d.node_y.push_back(0);
          d.node_size_x.push_back(10);
          d.node_size_y.push_back(10);
          d.node_macro.push_back(1);
        }
    d.num_movable = paths.size();
    // A unit cell defines the equivalent-area unit, outside all scopes.
    d.node_x.push_back(0);
    d.node_y.push_back(0);
    d.node_size_x.push_back(1);
    d.node_size_y.push_back(1);
    d.node_macro.push_back(0);
    paths.push_back("unit");
    masters.push_back("and");
    if (minimum == 4)
      return detail::centroid_candidates(d, paths, masters, 299);
    return detail::centroid_candidates(d, paths, masters, 299, minimum);
  };
  check(discover({3}).empty(), "default admits fewer than four siblings");
  auto four = discover({4});
  check(four.size() == 4, "four siblings under one parent rejected");
  for (const auto &g : four)
    check(g.copy_parent == "top.p0" && g.copies_per_parent == 4 &&
              g.scope_macros == 3 && g.equivalent_cells == 300,
          "wrong immediate parent/copy census");
  check(discover({4}, 5).empty(), "copy threshold parameter ignored");
  auto three = discover({3}, 3);
  check(three.size() == 3, "three-copy override rejected");
  check(discover({2, 2}).empty(), "copies pooled across two parents");
  check(discover({3, 3}).empty(), "six global copies replace four per parent");
  auto local = discover({4, 3});
  check(local.size() == 4, "underfilled parent affected qualifying parent");
  for (const auto &g : local)
    check(g.copy_parent == "top.p0" && g.copies_per_parent == 4,
          "three-copy parent admitted using copies elsewhere");
  auto wrapped = discover({4}, 4, true);
  check(wrapped.size() == 4, "unary alias hides eligible outer siblings");
  for (const auto &g : wrapped)
    check(g.name.find(".wrapper") == std::string::npos &&
              g.copy_parent == "top.p0",
          "singleton wrapper chosen instead of repeated outer scope");
  auto nested = discover({4, 4, 4, 4});
  check(nested.size() == 20, "nested four-copy scopes omitted");
  int roots = 0, children = 0;
  for (const auto &g : nested) {
    check(g.copies_per_parent == 4, "nested copy census changed");
    if (g.parent < 0)
      ++roots;
    else {
      ++children;
      check(g.name.starts_with(nested[g.parent].name + "."),
            "nested selected parent index invalid");
    }
  }
  check(roots == 4 && children == 16, "nested hierarchy relations lost");
  for (int bad : {-1, 0, 1})
    rejects([&] { discover({4}, bad); });
}
void repulsion_checks(const Design &d) {
  auto g = detail::centroid_repulsion(d, .3, .01);
  check(g.value > 0, "overlapping same-family circles have no energy");
  double sx = 0, sy = 0;
  for (int n = 0; n < 9; ++n) {
    sx += g.x[n];
    sy += g.y[n];
    for (int axis = 0; axis < 2; ++axis) {
      auto lo = d, hi = d;
      (axis ? lo.node_y : lo.node_x)[n] -= 1e-3;
      (axis ? hi.node_y : hi.node_x)[n] += 1e-3;
      const double fd = (detail::centroid_repulsion(hi, .3, .01).value -
                         detail::centroid_repulsion(lo, .3, .01).value) /
                        .002;
      near(fd, axis ? g.y[n] : g.x[n], 1e-5,
           "shared-color repulsion gradient mismatch");
    }
  }
  near(sx, 0, 1e-9, "centroid pair repulsion loses x force balance");
  near(sy, 0, 1e-9, "centroid pair repulsion loses y force balance");
  auto other = d;
  other.centroid_groups[1].family = 1;
  const auto off = detail::centroid_repulsion(other, .3, .01);
  check(off.value == 0 && off.x == std::vector<double>(9) &&
            off.y == std::vector<double>(9),
        "different types repel");
  auto shifted = d;
  for (int n = 0; n < 9; ++n) {
    shifted.node_x[n] += 71;
    shifted.node_y[n] -= 22;
  }
  near(detail::centroid_repulsion(shifted, .3, .01).value, g.value, 1e-8,
       "repulsion uses an absolute anchor");
  auto padded = d;
  for (int n = 0; n < 6; ++n) {
    padded.node_x[n] -= 10;
    padded.node_y[n] -= 10;
    padded.node_size_x[n] += 20;
    padded.node_size_y[n] += 20;
  }
  near(detail::centroid_repulsion(padded, .3, .01).value, g.value, 1e-8,
       "repulsion counts halo as anchored physical area");
  near(detail::centroid_cell_area(padded), detail::centroid_cell_area(d), 0,
       "equivalent-cell area changes with halos");
  auto summary = detail::centroid_summary(d);
  near(summary[0].area, 3602, 0, "shared color outside area omitted/doubled");
  // Nested candidates: the parent contributes attraction to the entire child,
  // while only the two same-type children repel each other.
  auto nested = d;
  CentroidGroup parent;
  parent.name = "top";
  parent.family = 99;
  parent.colors = {0, 1, 2, 3, 4};
  auto child_a = d.centroid_groups[0], child_b = d.centroid_groups[1];
  child_a.parent = child_b.parent = 0;
  nested.centroid_groups = {parent, child_a, child_b};
  nested.validate();
  const auto nr = detail::centroid_repulsion(nested, .3, .01);
  near(nr.value, g.value, 1e-8, "parent repels its nested child");
  for (int n = 0; n < 9; ++n) {
    near(nr.x[n], g.x[n], 1e-9, "nested parent adds x repulsion");
    near(nr.y[n], g.y[n], 1e-9, "nested parent adds y repulsion");
  }
  const auto ns = detail::centroid_summary(nested);
  auto parent_only = nested;
  parent_only.centroid_groups.resize(1);
  const auto pg = detail::centroid_attraction(parent_only, .3);
  double pull_x = 0, pull_y = 0;
  for (int n : {0, 1, 2, 3, 6, 7}) {
    pull_x += pg.x[n];
    pull_y += pg.y[n];
  }
  near(pull_x, .6 * ns[1].area / ns[0].area * (ns[1].x - ns[0].x), 1e-9,
       "child centroid does not inherit parent x attraction");
  near(pull_y, .6 * ns[1].area / ns[0].area * (ns[1].y - ns[0].y), 1e-9,
       "child centroid does not inherit parent y attraction");

  Design two;
  two.layout = {0, 0, 1000, 1000, 1, 1};
  two.num_movable = 2;
  two.node_x = {10, 20};
  two.node_y = {10, 5};
  two.node_size_x = two.node_size_y = {10, 20};
  two.node_macro = {0, 1};
  two.centroid_node_area = {100, 400};
  two.node_centroid_color = {-1, -1};
  CentroidGroup a, b;
  a.nodes = {0};
  b.nodes = {1};
  two.centroid_groups = {a, b};
  two.build_maps();
  two.validate();
  // One standard cell + a four-equivalent-cell macro, using the same area
  // unit as candidate eligibility. Circle area is A / requested density.
  near(detail::centroid_cell_area(two), 100, 0, "incorrect cell-area unit");
  const double reach = 30 / std::sqrt(std::acos(-1.));
  const auto pair = detail::centroid_repulsion(two, 2, 1);
  near(pair.value, 2 * std::pow(reach - 15, 2), 1e-10,
       "circle radius is not sqrt(area / pi / density)");
  check(pair.x[0] > 0 && pair.x[1] == -pair.x[0],
        "gradient does not push centroids apart");
  near(detail::centroid_repulsion(two, 2, .25).value,
       2 * std::pow(2 * reach - 15, 2), 1e-9,
       "requested density does not scale circle area");
  two.node_x[1] = reach + 6; // distance reach+1
  const auto outside = detail::centroid_repulsion(two, 2, 1);
  check(outside.value == 0 && outside.x[0] == 0 && outside.x[1] == 0,
        "repulsion persists outside the sum of radii");
  two.node_x[1] = 5; // exactly coincident centers
  const auto tie = detail::centroid_repulsion(two, 2, 1);
  const auto repeat = detail::centroid_repulsion(two, 2, 1);
  check((tie.x[0] != 0 || tie.y[0] != 0) && tie.x[0] == -tie.x[1] &&
            tie.y[0] == -tie.y[1] && tie.x == repeat.x && tie.y == repeat.y,
        "coincident centroid tie is not deterministically separated");
  // Equal radii: force must persist when neither center is inside the
  // other's circle. It stops at TWO radii, not one.
  two.node_size_x = two.node_size_y = {10, 10};
  two.centroid_node_area = {100, 100};
  const double radius = std::sqrt(100 / std::acos(-1.));
  two.node_y = {10, 10};
  two.node_x = {10, 10 + 1.5 * radius};
  const auto between = detail::centroid_repulsion(two, 1, 1);
  check(between.value > 0 && between.x[0] > 0,
        "repulsion incorrectly stops at one radius");
  near(detail::centroid_circle_overlaps(two, 1).maximum, .25, 1e-12,
       "circle penetration uses one radius instead of radius sum");
  two.node_x[1] = 10 + 2.01 * radius;
  check(detail::centroid_repulsion(two, 1, 1).value == 0 &&
            detail::centroid_circle_overlaps(two, 1).overlaps == 0,
        "equal circles repel beyond two radii");

  // A density-legal pair cannot terminate while its centroid circles overlap.
  two.node_x = {10, 12};
  Options constrained;
  constrained.centroid_weight = .25;
  constrained.centroid_repulsion_weight = 5;
  constrained.centroid_adaptive_repulsion = true;
  constrained.centroid_ramp = 20;
  constrained.max_iterations = 400;
  constrained.bins_x = constrained.bins_y = 8;
  constrained.target_density = 1;
  constrained.add_fillers = constrained.legalize = constrained.detailed_place =
      false;
  // A net supplies a real reference force for the native force normalization.
  two.pin2node = {0, 1};
  two.pin_offset_x = two.pin_offset_y = {0, 0};
  two.flat_net2pin = {0, 1};
  two.flat_net2pin_start = {0, 2};
  two.build_maps();
  auto capped = two;
  constrained.max_iterations = 40;
  const auto cap = place(capped, constrained);
  check(!cap.converged || cap.centroid.circles_converged,
        "density-only convergence hides overlapping centroid circles");
  constrained.max_iterations = 400;
  auto spread = two;
  const auto separated = place(spread, constrained);
  check(separated.centroid.circles_converged && separated.converged,
        "adaptive penalty does not separate a feasible equal-circle pair");
  check(separated.centroid.effective_repulsion_weight >= 5 &&
            !separated.history.empty(),
        "missing circle penalty diagnostics");
  rejects([&] { detail::centroid_repulsion(two, 1, 0); });
  rejects([&] { detail::centroid_repulsion(two, -1, .5); });
}
void circle_attraction_checks(Design d) {
  // Exercise shared colors, nested anchors, explicit members and fixed nodes.
  CentroidGroup extra;
  extra.nodes = {2, 7, 8};
  extra.repulsion_area = 100;
  d.centroid_groups.push_back(extra);
  for (bool frozen : {false, true}) {
    if (frozen)
      d.centroid_targets = {{300, 400}, {500, 450}, {600, 650}};
    const auto g = detail::centroid_circle_attraction(d, .3, .5, 4);
    double sx = 0, sy = 0;
    for (size_t n = 0; n < d.node_x.size(); ++n) {
      sx += g.x[n];
      sy += g.y[n];
      for (int axis = 0; axis < 2; ++axis) {
        auto hi = d, lo = d;
        (axis ? hi.node_y : hi.node_x)[n] += .001;
        (axis ? lo.node_y : lo.node_x)[n] -= .001;
        const double fd =
            (detail::centroid_circle_attraction(hi, .3, .5, 4).value -
             detail::centroid_circle_attraction(lo, .3, .5, 4).value) /
            .002;
        near(fd, axis ? g.y[n] : g.x[n], 1e-4,
             "circle attraction finite-difference mismatch");
      }
    }
    if (!frozen) {
      near(sx, 0, 1e-9, "moving circle anchors x translation");
      near(sy, 0, 1e-9, "moving circle anchors y translation");
      auto shifted = d;
      for (size_t n = 0; n < d.node_x.size(); ++n) {
        shifted.node_x[n] += 99;
        shifted.node_y[n] -= 57;
      }
      near(detail::centroid_circle_attraction(shifted, .3, .5, 4).value,
           g.value, 1e-8, "moving circle depends on world origin");
    } else {
      check(std::hypot(sx, sy) > 1, "frozen anchor loses restoring force");
      const auto macros =
          detail::centroid_circle_attraction(d, .3, .5, 4, true);
      for (size_t n = 0; n < d.node_x.size(); ++n) {
        const bool active = n < size_t(d.num_movable) && d.node_macro[n];
        near(macros.x[n], active ? g.x[n] : 0, 0,
             "macro-only anchor changes macro strength or pulls cells");
        near(macros.y[n], active ? g.y[n] : 0, 0,
             "macro-only anchor changes macro strength or pulls cells");
        for (int axis = 0; axis < 2; ++axis) {
          auto hi = d, lo = d;
          (axis ? hi.node_y : hi.node_x)[n] += .001;
          (axis ? lo.node_y : lo.node_x)[n] -= .001;
          const double fd =
              (detail::centroid_circle_attraction(hi, .3, .5, 4, true).value -
               detail::centroid_circle_attraction(lo, .3, .5, 4, true).value) /
              .002;
          near(fd, axis ? macros.y[n] : macros.x[n], 1e-4,
               "macro-only anchor finite-difference mismatch");
        }
      }
    }
    auto padded = d;
    for (int n = 0; n < 6; ++n) {
      padded.node_x[n] -= 10;
      padded.node_y[n] -= 10;
      padded.node_size_x[n] += 20;
      padded.node_size_y[n] += 20;
    }
    near(detail::centroid_circle_attraction(padded, .3, .5, 4).value, g.value,
         1e-8, "halo changes circle attraction");
    check(detail::centroid_attraction_reference(padded, .5) ==
              detail::centroid_attraction_reference(d, .5),
          "halo changes spring calibration");
    const auto quadratic = detail::centroid_attraction(d, .3);
    const auto zero_outside = detail::centroid_circle_attraction(d, .3, .5, 0);
    near(quadratic.value, zero_outside.value, 1e-8,
         "zero outside stiffness is not quadratic");
    for (size_t n = 0; n < d.node_x.size(); ++n) {
      near(quadratic.x[n], zero_outside.x[n], 1e-9, "quadratic x mismatch");
      near(quadratic.y[n], zero_outside.y[n], 1e-9, "quadratic y mismatch");
    }
  }
  // Radius ten, one member: force is 2d inside and 2(d+4(d-R)) outside.
  d.centroid_groups = {extra};
  d.centroid_groups[0].nodes = {0};
  d.centroid_groups[0].repulsion_area = 50 * std::acos(-1.);
  d.centroid_targets = {{115, 215}};
  for (const auto &[distance, expected] :
       {std::pair{5., 10.}, {10., 20.}, {20., 120.}}) {
    d.node_x[0] = 100 + distance;
    const auto g = detail::centroid_circle_attraction(d, 1, .5, 4);
    near(g.x[0], expected, 1e-9, "restoring force does not grow with distance");
    near(g.y[0], 0, 1e-9, "radial spring has tangential force");
  }
  d.node_x[0] = 110 - 1e-7;
  const auto left = detail::centroid_circle_attraction(d, 1, .5, 4);
  d.node_x[0] = 110 + 1e-7;
  const auto right = detail::centroid_circle_attraction(d, 1, .5, 4);
  near(left.x[0], right.x[0], 2e-6, "force discontinuity at circle boundary");
  // A single member is its own moving centroid, even when moment arithmetic
  // leaves a tiny residual. There must be no self-pull to amplify in GP.
  d.centroid_targets.clear();
  d.node_x[0] = 1398357.29110317;
  d.node_y[0] = 741953.81801293;
  const auto singleton = detail::centroid_circle_attraction(d, 1, .5, 4);
  near(singleton.x[0], 0, 0, "singleton has moving-centroid x self-pull");
  near(singleton.y[0], 0, 0, "singleton has moving-centroid y self-pull");
  rejects([&] { detail::centroid_circle_attraction(d, 1, .5, -1); });
  rejects([&] { detail::centroid_attraction_reference(d, 0); });
}
void fixed_spring_gp() {
  Design d;
  d.layout = {0, 0, 1000, 1000, 1, 1};
  d.num_movable = 1;
  d.node_x = {95, 900};
  d.node_y = {495, 500};
  d.node_size_x = d.node_size_y = {10, 1};
  d.node_macro = {1, 0};
  d.centroid_node_area = {100, 1};
  d.node_centroid_color = {-1, -1};
  CentroidGroup group;
  group.nodes = {0};
  group.repulsion_area = 5000 * std::acos(-1.); // R=100 at density .5
  d.centroid_groups = {group};
  d.centroid_targets = {{100, 500}};
  d.pin2node = {0, 1};
  d.pin_offset_x = d.pin_offset_y = {5, 0};
  d.flat_net2pin = {0, 1};
  d.flat_net2pin_start = {0, 2};
  d.build_maps();
  Options o;
  o.max_iterations = 200;
  o.bins_x = o.bins_y = 8;
  o.add_fillers = o.legalize = o.detailed_place = false;
  o.density_upper_bound = true;
  o.target_density = .5;
  o.centroid_weight = .25;
  o.centroid_ramp = 10000; // frozen spring must not restart this ramp
  auto old = d;
  place(old, o);
  auto disabled = d;
  o.centroid_attraction_outside = 999;
  place(disabled, o);
  check(old.node_x == disabled.node_x && old.node_y == disabled.node_y,
        "new controls change disabled legacy path");
  o.centroid_distance_attraction = true;
  o.centroid_attraction_outside = 4;
  const auto result = place(d, o);
  const double drift = std::hypot(d.node_x[0] - 95, d.node_y[0] - 495);
  const double old_drift = std::hypot(old.node_x[0] - 95, old.node_y[0] - 495);
  std::cout << "frozen spring GP drift " << drift << " legacy " << old_drift
            << " calibration " << result.centroid.attraction_coefficient
            << '\n';
  check(drift < 200 && old_drift > drift * 2,
        "fixed spring does not bound GP drift");
  check(result.centroid.attraction_coefficient > 0 &&
            result.centroid.attraction_calibration_iteration == 0,
        "spring calibration missing or restarted");
  check(d.centroid_targets == old.centroid_targets && d.node_x[1] == 900,
        "target/fixed node moved");
  o.centroid_attraction_outside = -1;
  rejects([&] { place(d, o); });
}
void parent_attraction_checks() {
  auto d = fixture();
  d.centroid_groups[0].colors = {0, 1, 2, 3, 4};
  d.centroid_groups[0].repulsion_area = std::acos(-1.) * .5 * 40 * 40;
  d.centroid_groups[1].parent = 0;
  d.centroid_groups[1].family += 1;
  d.validate();
  std::uint64_t edges = 0;
  double distance = 0;
  const auto g =
      detail::centroid_parent_attraction(d, 2, .5, 4, &edges, &distance);
  check(edges == 1 && distance > 1 && g.value > 0, "parent edge missing");
  double sx = 0, sy = 0;
  for (int n = 0; n < 9; ++n) {
    sx += g.x[n];
    sy += g.y[n];
    for (int axis = 0; axis < 2; ++axis) {
      auto lo = d, hi = d;
      (axis ? lo.node_y : lo.node_x)[n] -= 1e-3;
      (axis ? hi.node_y : hi.node_x)[n] += 1e-3;
      const double fd =
          (detail::centroid_parent_attraction(hi, 2, .5, 4).value -
           detail::centroid_parent_attraction(lo, 2, .5, 4).value) /
          .002;
      near(fd, axis ? g.y[n] : g.x[n], 1e-5 * std::max(1., std::abs(fd)),
           "parent edge chain rule mismatch with nested/shared members");
    }
  }
  near(sx, 0, 1e-7, "parent edge x force not balanced");
  near(sy, 0, 1e-7, "parent edge y force not balanced");
  auto shifted = d;
  for (int n = 0; n < 9; ++n) {
    shifted.node_x[n] += 30;
    shifted.node_y[n] -= 20;
  }
  near(detail::centroid_parent_attraction(shifted, 2, .5, 4).value, g.value,
       1e-6, "parent edge pins absolute position");
  auto inside = d;
  inside.centroid_groups[0].repulsion_area = std::acos(-1.) * .5 * 1000 * 1000;
  const auto weak = detail::centroid_parent_attraction(inside, 2, .5, 4);
  const auto quadratic = detail::centroid_parent_attraction(inside, 2, .5, 0);
  check(weak.x == quadratic.x && weak.y == quadratic.y,
        "outside penalty active inside parent circle");
  const auto outside_linear = detail::centroid_parent_attraction(d, 2, .5, 0);
  check(std::abs(g.x[4]) > 10 * std::abs(outside_linear.x[4]),
        "far child has no progressive restoring force");
  auto unrelated = d;
  unrelated.centroid_groups[1].parent = -1;
  check(detail::centroid_parent_attraction(unrelated, 2, .5, 4).value == 0,
        "unrelated centroids acquired parent attraction");
  Options o;
  o.max_iterations = 80;
  o.bins_x = o.bins_y = 8;
  o.add_fillers = o.legalize = o.detailed_place = false;
  o.centroid_weight = .25;
  o.centroid_ramp = 20;
  auto baseline = d, disabled = d;
  place(baseline, o);
  o.centroid_parent_outside = 8;
  place(disabled, o);
  check(baseline.node_x == disabled.node_x &&
            baseline.node_y == disabled.node_y,
        "zero parent weight changes GP");
  o.centroid_parent_weight = 1;
  auto enabled = d;
  const auto result = place(enabled, o);
  bool active = false;
  for (const auto &sample : result.history)
    active |=
        sample.centroid_parent_edges == 1 && sample.centroid_parent_force > 0;
  check(active && result.centroid.parent_seconds > 0,
        "parent force not applied during GP");
  check(enabled.node_x[8] == d.node_x[8] && enabled.node_y[8] == d.node_y[8],
        "parent force moved fixed node");
  o.centroid_parent_weight = -1;
  rejects([&] { place(enabled, o); });
  d.centroid_targets = {{1, 2}, {3, 4}};
  rejects([&] { detail::centroid_parent_attraction(d, 1, .5, 4); });
}
void alignment_checks() {
  auto d = fixture();
  for (int i = 3; i < 6; ++i)
    d.node_y[i] -= 300;
  d.node_y[4] += 12;
  std::uint64_t pairs = 0;
  double residual = 0;
  for (double sigma : {0., .25}) {
    const auto g =
        detail::centroid_alignment(d, 2, 3, sigma, .5, &pairs, &residual);
    check(pairs == 1 && residual > 0 && g.value > 0,
          "centroid row pair missing");
    double sx = 0, sy = 0;
    for (int n = 0; n < 9; ++n) {
      sx += g.x[n];
      sy += g.y[n];
      for (int axis = 0; axis < 2; ++axis) {
        auto lo = d, hi = d;
        (axis ? lo.node_y : lo.node_x)[n] -= 1e-3;
        (axis ? hi.node_y : hi.node_x)[n] += 1e-3;
        const double fd =
            (detail::centroid_alignment(hi, 2, 3, sigma, .5).value -
             detail::centroid_alignment(lo, 2, 3, sigma, .5).value) /
            .002;
        near(fd, axis ? g.y[n] : g.x[n], 1e-7,
             "shared-member alignment chain rule mismatch");
      }
      near(g.x[n], 0, 0, "row alignment applies an x force");
    }
    near(sx, 0, 1e-10, "alignment x force imbalance");
    near(sy, 0, 1e-10, "alignment y force imbalance");
    auto shifted = d;
    for (int n = 0; n < 9; ++n) {
      shifted.node_x[n] += 50;
      shifted.node_y[n] -= 70;
    }
    near(detail::centroid_alignment(shifted, 2, 3, sigma, .5).value, g.value,
         1e-9, "centroid alignment pins absolute locations");
    auto padded = d;
    for (int n = 0; n < 6; ++n) {
      padded.node_x[n] -= 5;
      padded.node_y[n] -= 5;
      padded.node_size_x[n] += 10;
      padded.node_size_y[n] += 10;
    }
    near(detail::centroid_alignment(padded, 2, 3, sigma, .5).value, g.value,
         1e-9, "halo changes centroid alignment circle");
  }
  auto column = d;
  std::swap(column.node_x, column.node_y);
  detail::centroid_alignment(column, 2, 3, .25, .5, &pairs);
  check(pairs == 1, "column alignment pair missing");
  auto other = d;
  other.centroid_groups[1].family += 1;
  check(detail::centroid_alignment(other, 2, 3, .25, .5).value == 0,
        "alignment crosses candidate types");
  other = d;
  other.centroid_groups[1].parent = 0;
  check(detail::centroid_alignment(other, 2, 3, .25, .5).value == 0,
        "alignment crosses parent scopes");
  detail::centroid_alignment(d, 2, 3, .25, .001, &pairs);
  check(pairs == 0, "alignment competes with overlapping circles");
  auto far = d;
  for (int n = 3; n < 6; ++n)
    far.node_y[n] += 1000;
  detail::centroid_alignment(far, 2, 3, .25, .5, &pairs);
  check(pairs == 0, "non-overlapping projections attract");
  detail::Gradient wire, density;
  wire.x.assign(9, 1);
  wire.y.assign(9, 2);
  density.x.assign(9, 0);
  density.y.assign(9, 0);
  const auto reference = detail::centroid_force_reference(d, wire, density, 1);
  near(reference.first, 6, 0, "centroid reference is not total member force");
  near(reference.second, 12, 0, "centroid force reference y mismatch");
  Options o;
  o.max_iterations = 80;
  o.bins_x = o.bins_y = 8;
  o.add_fillers = o.legalize = o.detailed_place = false;
  o.centroid_weight = .25;
  o.centroid_repulsion_weight = 5;
  o.centroid_distance_attraction = true;
  o.centroid_ramp = 20;
  auto baseline = d, disabled = d;
  place(baseline, o);
  o.centroid_align_sigma = .5;
  o.centroid_align_start = 0;
  o.centroid_align_full = 60;
  place(disabled, o);
  check(disabled.node_x == baseline.node_x &&
            disabled.node_y == baseline.node_y,
        "zero centroid alignment changes existing FP");
  o.centroid_align_weight = 25;
  auto enabled = d;
  const auto result = place(enabled, o);
  check(result.iterations > 60 && result.centroid.alignment_seconds > 0,
        "centroid alignment ramp or accounting omitted");
  bool force = false;
  for (const auto &sample : result.history)
    force |= sample.centroid_align_force > 0;
  check(force, "centroid alignment has no applied force");
  check(enabled.node_x[8] == d.node_x[8] && enabled.node_y[8] == d.node_y[8],
        "centroid alignment moved fixed node");
  o.centroid_align_weight = -1;
  rejects([&] { place(enabled, o); });
  o.centroid_align_weight = 25;
  o.centroid_align_full = o.centroid_align_start;
  rejects([&] { place(enabled, o); });
  other = d;
  other.centroid_targets = {{1, 2}, {3, 4}};
  rejects([&] { detail::centroid_alignment(other, 2, 3, .25, .5); });
}

int main() {
  try {
    parent_attraction_checks();
    alignment_checks();
    candidate_copy_checks();
    auto d = fixture();
    repulsion_checks(d);
    circle_attraction_checks(d);
    fixed_spring_gp();
    auto grad = detail::centroid_attraction(d, .3);
    auto a = d, b = d;
    a.centroid_groups.resize(1);
    b.centroid_groups.erase(b.centroid_groups.begin());
    auto ga = detail::centroid_attraction(a, .3),
         gb = detail::centroid_attraction(b, .3);
    for (int n = 0; n < 9; ++n) {
      near(grad.x[n], ga.x[n] + gb.x[n], 1e-10,
           "shared color loses one of its anchors");
      near(grad.y[n], ga.y[n] + gb.y[n], 1e-10,
           "shared color loses one of its anchors");
    }
    check(ga.x[7] != 0 && gb.x[7] != 0,
          "outside standard cell does not get both candidate anchors");
    // A third overlapping candidate contributes a third anchor, with no
    // ownership.
    auto triple = d;
    triple.centroid_groups.push_back(d.centroid_groups[0]);
    auto gt = detail::centroid_attraction(triple, .3);
    near(gt.x[7], grad.x[7] + ga.x[7], 1e-10,
         "n anchors reduced to at most two");
    double sx = 0, sy = 0;
    for (int n = 0; n < 9; ++n) {
      sx += grad.x[n];
      sy += grad.y[n];
    }
    near(sx, 0, 1e-9, "absolute x anchor introduced");
    near(sy, 0, 1e-9, "absolute y anchor introduced");
    auto shifted = d;
    for (int n = 0; n < 9; ++n) {
      shifted.node_x[n] += 70;
      shifted.node_y[n] -= 20;
    }
    near(detail::centroid_attraction(shifted, .3).value, grad.value, 1e-8,
         "centroid is tied to an absolute coordinate");
    for (int axis = 0; axis < 2; ++axis)
      for (int n = 0; n < 8; ++n) {
        auto hi = d, lo = d;
        double eps = 1e-3;
        (axis ? hi.node_y : hi.node_x)[n] += eps;
        (axis ? lo.node_y : lo.node_x)[n] -= eps;
        double fd = (detail::centroid_attraction(hi, .3).value -
                     detail::centroid_attraction(lo, .3).value) /
                    (2 * eps);
        near(fd, axis ? grad.y[n] : grad.x[n], 1e-5,
             "finite-difference gradient mismatch");
      }
    auto padded = d;
    for (int n = 0; n < 6; ++n) {
      padded.node_x[n] -= 10;
      padded.node_y[n] -= 10;
      padded.node_size_x[n] += 20;
      padded.node_size_y[n] += 20;
    }
    near(detail::centroid_attraction(padded, .3).value, grad.value, 1e-8,
         "halo changes physical-area weights");
    auto bad = d;
    bad.centroid_groups[0].colors.push_back(0);
    rejects([&] { bad.validate(); });
    bad = d;
    bad.centroid_groups[0].nodes.push_back(0);
    rejects([&] { bad.validate(); });
    Options o;
    o.max_iterations = 40;
    o.bins_x = o.bins_y = 8;
    o.add_fillers = false;
    o.legalize = o.detailed_place = false;
    o.centroid_ramp = 20;
    auto plain = d;
    plain.centroid_groups.clear();
    place(plain, o);
    auto off = d;
    place(off, o);
    check(off.node_x == plain.node_x && off.node_y == plain.node_y,
          "disabled centroid changes GP");
    auto on = d;
    o.centroid_weight = .25;
    bool started = false;
    o.observer = [&](const Design &snapshot, PlacementStage stage, int iter) {
      if (stage == PlacementStage::centroid_start) {
        started = true;
        check(iter == 0, "R8H2 uses a preliminary GP");
        check(snapshot.node_x == d.node_x && snapshot.node_y == d.node_y,
              "R8H2 does not start from initial positions");
      }
    };
    auto result = place(on, o);
    check(started && result.centroid.completed &&
              result.centroid.iterations == result.iterations,
          "centroid force not active throughout first GP");
    check(on.node_x[8] == d.node_x[8] && on.node_y[8] == d.node_y[8],
          "fixed node moved");
    auto separated = d;
    o.centroid_repulsion_weight = 5;
    place(separated, o);
    check(separated.node_x[8] == d.node_x[8] &&
              separated.node_y[8] == d.node_y[8],
          "repulsion moved fixed cell");
    o.centroid_repulsion_weight = -1;
    rejects([&] { place(separated, o); });
    o.centroid_repulsion_weight = 0;
    GpState state;
    rejects([&] { place(on, o, &state); });
    o.centroid_weight = -1;
    rejects([&] { place(on, o); });
    std::cout
        << "candidate thresholds, all-color centroids, additive n anchors, "
           "outside cells, gradient, halo, same-type density circles and "
           "from-scratch GP pass\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
