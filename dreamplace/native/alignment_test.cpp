// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("alignment assertion");
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
int main() {
  try {
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
