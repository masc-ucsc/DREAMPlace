// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <set>
#include <source_location>
#include <stdexcept>
#include <tuple>
using namespace DreamPlace::native;
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("alignment assertion at line " +
                             std::to_string(at.line()));
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 80, 80, 1, 1};
  d.num_movable = 5;
  d.node_x = {10, 11, 20, 11, 11, 11};
  d.node_y = {10, 20, 11, 11, 30, 40};
  d.node_size_x.assign(6, 4);
  d.node_size_y.assign(6, 4);
  d.node_macro = {1, 1, 1, 1, 1, 1};
  d.node_align_class = {0, 0, 0, 0, 1, 0};
  d.build_maps();
  return d;
}
// Compare a strip-index search with an independent exhaustive geometry oracle.
void strip_pairs() {
  Design d;
  d.layout = {0, 0, 1000, 1000, 1, 1};
  d.num_movable = 260;
  for (int i = 0; i < 270; ++i) {
    d.node_x.push_back((i * 37) % 473 + .25 * (i % 3));
    d.node_y.push_back((i * 83) % 479 + .25 * (i % 5));
    d.node_size_x.push_back(7 + (i % 4) * 6);
    d.node_size_y.push_back(9 + (i % 4) * 4);
    d.node_macro.push_back(i % 17 != 0);
    d.node_align_class.push_back(i % 13 == 0 ? -1 : i % 4);
  }
  d.build_maps();
  using Pair = std::tuple<int, int, Axis>;
  std::set<Pair> expected;
  for (int a = 0; a < d.num_movable; ++a)
    for (int b = a + 1; b < d.num_movable; ++b) {
      if (!d.node_macro[a] || !d.node_macro[b] || d.node_align_class[a] < 0 ||
          d.node_align_class[a] != d.node_align_class[b])
        continue;
      const bool x = d.node_x[a] < d.node_x[b] + d.node_size_x[b] &&
                     d.node_x[b] < d.node_x[a] + d.node_size_x[a];
      const bool y = d.node_y[a] < d.node_y[b] + d.node_size_y[b] &&
                     d.node_y[b] < d.node_y[a] + d.node_size_y[a];
      if (x != y)
        expected.emplace(a, b, x ? Axis::x : Axis::y);
    }
  const auto collect = [](const auto &pairs) {
    std::set<Pair> result;
    for (const auto &p : pairs)
      check(result.emplace(p.a, p.b, p.axis).second);
    return result;
  };
  check(!expected.empty());
  check(collect(detail::macro_grid_align_pairs(d)) == expected);
  // Padding does not create attraction across a non-overlapping physical edge.
  auto padded = d;
  padded.node_pad_x.assign(270, 2);
  padded.node_pad_y.assign(270, 3);
  for (int i = 0; i < 270; ++i) {
    padded.node_x[i] -= 2;
    padded.node_y[i] -= 3;
    padded.node_size_x[i] += 4;
    padded.node_size_y[i] += 6;
  }
  check(collect(detail::macro_grid_align_pairs(padded)) == expected);
  const auto pairs = detail::macro_grid_align_pairs(d);
  const auto a = detail::macro_alignment(d, pairs, 1, 2, .25, nullptr, true);
  const auto b =
      detail::macro_alignment(padded, pairs, 1, 2, .25, nullptr, true);
  check(a.x == b.x && a.y == b.y && a.value == b.value);
  // Adjacent strips: centers straddle a boundary, and must still attract.
  d = fixture();
  d.num_movable = 2;
  d.node_x[0] = 1.9;
  d.node_x[1] = 2.1;
  auto near = detail::macro_grid_align_pairs(d);
  check(near.size() == 1 && near[0].axis == Axis::x);
  auto force = detail::macro_alignment(d, near, 1, 1, .25, nullptr, true);
  check(force.x[0] < 0 && force.x[1] > 0 && force.y[0] == 0 && force.y[1] == 0);
  d.node_x[1] = d.node_x[0] + 4;
  check(detail::macro_grid_align_pairs(d)
            .empty()); // just touching is not overlap
  d.node_x[1] = d.node_x[0];
  d.node_y[1] = d.node_y[0];
  check(detail::macro_grid_align_pairs(d).empty()); // overlap on both axes
  d.node_y[1] += 10;
  d.node_align_class[1] = 1;
  check(detail::macro_grid_align_pairs(d).empty()); // other master
}
int main() {
  try {
    strip_pairs();
    auto d = fixture();
    auto pairs = detail::macro_align_pairs(d);
    auto has = [&](int a, int b, Axis axis) {
      for (auto p : pairs)
        if (p.a == a && p.b == b && p.axis == axis)
          return true;
      return false;
    };
    check(has(0, 1, Axis::x));
    check(has(0, 2, Axis::y));
    for (auto p : pairs)
      check(!(p.a == 0 && p.b == 3) && p.a != 4 && p.b != 4 && p.a != 5 &&
            p.b != 5);
    auto near = detail::macro_align_pairs(d, true);
    check(!near.empty() && near.size() <= pairs.size());
    // Fixed pair list: gradients do not differentiate the discrete refresh.
    for (double sigma : {0., .25}) {
      auto g = detail::macro_alignment(d, pairs, 2, 3, sigma);
      double sx = 0, sy = 0;
      for (int i = 0; i < 5; ++i) {
        sx += g.x[i];
        sy += g.y[i];
        for (bool y : {false, true}) {
          auto lo = d, hi = d;
          (y ? lo.node_y : lo.node_x)[i] -= 1e-5;
          (y ? hi.node_y : hi.node_x)[i] += 1e-5;
          double fd = (detail::macro_alignment(hi, pairs, 2, 3, sigma).value -
                       detail::macro_alignment(lo, pairs, 2, 3, sigma).value) /
                      2e-5;
          check(std::abs(fd - (y ? g.y[i] : g.x[i])) < 1e-7);
        }
      }
      check(std::abs(sx) < 1e-10 && std::abs(sy) < 1e-10 && g.x[5] == 0 &&
            g.y[5] == 0);
    }
    Options o;
    o.macro_align_start = 100;
    o.macro_align_full = 500;
    check(detail::macro_align_ramp(o, 99) == 0);
    check(detail::macro_align_ramp(o, 100) == 0);
    check(detail::macro_align_ramp(o, 300) == .25);
    check(detail::macro_align_ramp(o, 500) == 1);
    check(detail::macro_align_ramp(o, 700) == 1);
    // Nearest-only drops the distant partner when a nearer one exists each
    // side.
    d.num_movable = 3;
    d.node_x = {10, 11, 12};
    d.node_y = {0, 10, 20};
    d.node_size_x.assign(3, 4);
    d.node_size_y.assign(3, 4);
    d.node_macro.assign(3, 1);
    d.node_align_class.assign(3, 0);
    d.build_maps();
    check(detail::macro_align_pairs(d).size() == 3);
    check(detail::macro_align_pairs(d, true).size() == 2);
    check(detail::macro_grid_align_pairs(d, true).size() == 2);
    for (auto invalid : {-2, 0}) {
      auto bad = d;
      if (invalid == -2)
        bad.node_align_class[0] = -2;
      else
        bad.node_align_class.pop_back();
      bool rejected = false;
      try {
        bad.validate();
      } catch (const std::invalid_argument &) {
        rejected = true;
      }
      check(rejected);
    }
    for (int which = 0; which < 3; ++which) {
      Options bad;
      bad.macro_align_weight = 1;
      if (which == 0)
        bad.macro_align_sigma = -1;
      if (which == 1)
        bad.macro_align_full = bad.macro_align_start;
      if (which == 2)
        bad.macro_align_overlap_boost = -1;
      bool rejected = false;
      try {
        place(d, bad);
      } catch (const std::invalid_argument &) {
        rejected = true;
      }
      check(rejected);
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
