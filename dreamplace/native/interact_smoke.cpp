// BSD-3-Clause. See LICENSE.
// Options::interact: neutral until it pins, pins hold cells and macros, releases free them.
#include <cmath>
#include <iostream>
#include <source_location>
#include <stdexcept>

#include "dreamplace/native/placement.hpp"
using namespace DreamPlace::native;
void check(bool ok, std::source_location where = std::source_location::current()) {
  if (!ok) {
    throw std::runtime_error("interact assertion failed at line " + std::to_string(where.line()));
  }
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 80, 40, 1, 1};
  d.num_movable = 4;  // three cells and one movable macro, then two terminals
  d.node_x = {10, 20, 30, 40, 0, 80};
  d.node_y = {10, 20, 10, 20, 20, 20};
  d.node_size_x = {2, 2, 2, 8, 0, 0};
  d.node_size_y = {1, 1, 1, 6, 0, 0};
  d.node_macro = {0, 0, 0, 1, 0, 0};
  d.pin2node = {4, 0, 0, 1, 1, 2, 2, 3, 3, 5};
  d.pin_offset_x.assign(10, 0);
  d.pin_offset_y.assign(10, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  d.flat_net2pin_start = {0, 2, 4, 6, 8, 10};
  d.build_maps();
  return d;
}
Options base() {
  Options o;
  o.max_iterations = 40;
  o.schedule_iterations = 40;
  o.legalize = o.detailed_place = false;
  o.bins_x = o.bins_y = 8;
  o.stop_overflow = 0;
  return o;
}
int main() {
  try {
    auto plain = fixture();
    place(plain, base());
    // Installed but idle: identical to no hook at all.
    auto idle = fixture();
    auto o = base();
    int calls = 0, last = 0, finishes = 0;
    o.interact_interval = 5;
    o.interact = [&](const InteractView &v) {
      if (v.first) {
        check(v.iteration == 1);
      } else if (v.finishing) {
        check(v.iteration == 40);
        ++finishes;
      } else {
        ++calls;
        last = v.iteration;
      }
      check(v.count == 6 && v.design && v.x && v.y);
      return InteractUpdate{};
    };
    place(idle, o);
    check(calls == 8 && last == 40 && finishes == 1);
    for (int i = 0; i < 6; ++i) {
      check(idle.node_x[i] == plain.node_x[i] && idle.node_y[i] == plain.node_y[i]);
    }
    // Pin a cell and the movable macro at iteration 10; both stay put until released.
    auto held = fixture();
    o = base();
    o.interact_interval = 10;
    o.interact = [&](const InteractView &v) {
      InteractUpdate u;
      if (v.iteration == 10) {
        u.pins = {{0, 12, 14, 0}, {3, 50, 12, 0}};
      } else if (v.iteration == 20 || v.iteration == 30) {
        check(v.x[0] == 12 && v.y[0] == 14 && v.x[3] == 50 && v.y[3] == 12);
        if (v.iteration == 30) {
          u.releases = {0};
        }
      }
      return u;
    };
    place(held, o);
    check(held.node_x[3] == 50 && held.node_y[3] == 12);  // macro stayed pinned to the end
    check(held.node_x[0] != 12 || held.node_y[0] != 14);  // the released cell moved again
    check(std::abs(held.node_x[1] - plain.node_x[1]) + std::abs(held.node_y[1] - plain.node_y[1]) > 0);  // others reacted
    // The host extends a solve that is about to stop; the stop rules wait for the extra steps.
    auto more = fixture();
    o = base();
    o.interact_interval = 10;
    int finishing = 0;
    o.interact = [&](const InteractView &v) {
      InteractUpdate u;
      if (v.finishing) {
        ++finishing;
        check(v.iteration == 40 || v.iteration == 55);
        u.extend = finishing == 1 ? 15 : 0;  // once
      }
      return u;
    };
    const auto result = place(more, o);
    check(finishing == 2 && result.iterations == 55);
    std::cout << "interact smoke ok\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
