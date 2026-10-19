// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("orientation assertion failed");
}
template <typename F> void rejects(F f) {
  bool bad = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    bad = true;
  }
  check(bad);
}
Design design() {
  Design d;
  d.layout = {0, 0, 40, 40, 1, 1};
  d.num_movable = 1;
  d.node_x = {10, 30, 0};
  d.node_y = {10, 30, 20};
  d.node_size_x = {12, 4, 0};
  d.node_size_y = {8, 6, 0};
  d.pin2node = {0, 0, 1, 2};
  d.pin_offset_x = {1, 9, 2, 0};
  d.pin_offset_y = {2, 6, 3, 0};
  d.flat_net2pin = {0, 2, 1, 3};
  d.flat_net2pin_start = {0, 2, 4};
  d.build_maps();
  return d;
}
} // namespace
int main() {
  try {
    const std::array all{Orientation::N, Orientation::FN, Orientation::FS,
                         Orientation::S};
    for (auto a : all) {
      check(compose(a, a) == Orientation::N);
      for (auto b : all)
        for (auto c : all)
          check(compose(a, compose(b, c)) == compose(compose(a, b), c));
      for (int x = 0; x <= 12; ++x)
        for (int y = 0; y <= 8; ++y) {
          auto p = orient_pin(x, y, 12, 8, a);
          check(p.first >= 0 && p.first <= 12 && p.second >= 0 &&
                p.second <= 8);
          auto back = orient_pin(p.first, p.second, 12, 8, a);
          check(back.first == x && back.second == y);
          for (auto b : all) {
            auto q = orient_pin(p.first, p.second, 12, 8, b);
            check(q == orient_pin(x, y, 12, 8, compose(a, b)));
          }
        }
    }
    auto d = design();
    const auto original = d;
    for (auto a : all) {
      set_orientation(d, 0, a);
      const auto offset = orient_pin(1, 2, 12, 8, a);
      check(d.pin_offset_x[0] == offset.first &&
            d.pin_offset_y[0] == offset.second);
      const auto once = d;
      set_orientation(d, 0, a);
      check(d.pin_offset_x == once.pin_offset_x &&
            d.pin_offset_y == once.pin_offset_y);
      check(d.node_x == original.node_x && d.node_y == original.node_y);
      check(d.node_size_x == original.node_size_x &&
            d.node_size_y == original.node_size_y);
      check(d.pin_offset_x[2] == original.pin_offset_x[2]);
    }
    set_orientation(d, 0, Orientation::N);
    check(d.pin_offset_x == original.pin_offset_x &&
          d.pin_offset_y == original.pin_offset_y);
    set_orientation(d, 0, Orientation::S);
    set_orientation(d, 1, Orientation::FN);
    const auto oriented = d;
    Options o;
    o.max_iterations = 0;
    o.detailed_place = false;
    check(place(d, o).legal);
    check(d.node_orientation == oriented.node_orientation &&
          d.pin_offset_x == oriented.pin_offset_x &&
          d.pin_offset_y == oriented.pin_offset_y);
    check(d.node_x[1] == oriented.node_x[1] &&
          d.node_y[1] == oriented.node_y[1]);
    rejects([&] { set_orientation(d, 2, Orientation::S); });
    rejects([&] { set_orientation(d, 0, static_cast<Orientation>(4)); });
    rejects([&] { set_orientation(d, -1, Orientation::N); });
    check(d.pin_offset_x == oriented.pin_offset_x &&
          d.pin_offset_y == oriented.pin_offset_y);
    // The unchecked variant performs the same transform and rejections.
    {
      auto checked = design(), fast = checked;
      for (auto a : all)
        for (auto b : all) {
          set_orientation(checked, 0, a);
          set_orientation_unchecked(fast, 0, a);
          set_orientation(checked, 0, b);
          set_orientation_unchecked(fast, 0, b);
          check(checked.pin_offset_x == fast.pin_offset_x &&
                checked.pin_offset_y == fast.pin_offset_y &&
                checked.node_orientation == fast.node_orientation);
        }
      const auto before = fast;
      rejects([&] { set_orientation_unchecked(fast, 2, Orientation::S); });
      rejects([&] {
        set_orientation_unchecked(fast, 0, static_cast<Orientation>(4));
      });
      rejects([&] { set_orientation_unchecked(fast, 3, Orientation::N); });
      check(fast.pin_offset_x == before.pin_offset_x &&
            fast.pin_offset_y == before.pin_offset_y &&
            fast.node_orientation == before.node_orientation);
    }
    auto invalid = d;
    invalid.node_orientation.pop_back();
    rejects([&] { invalid.validate(); });
    invalid = d;
    invalid.node_orientation[0] = static_cast<Orientation>(255);
    rejects([&] { invalid.validate(); });
    // Default seed remains the historical 1000. Explicit same seed repeats;
    // another filler seed changes density forces without modifying fixed nodes.
    auto a = design(), b = a, c = a;
    o.max_iterations = 1;
    o.legalize = false;
    o.bins_x = o.bins_y = 8;
    place(a, o);
    o.seed = 1000;
    place(b, o);
    o.seed = 7;
    place(c, o);
    check(a.node_x == b.node_x && a.node_y == b.node_y);
    check(a.node_x != c.node_x || a.node_y != c.node_y);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
