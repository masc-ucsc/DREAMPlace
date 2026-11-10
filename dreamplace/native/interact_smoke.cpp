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
    // Centroid targets: dragging a group's target pulls its members; release restores the moving centroids.
    {
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
      const std::vector<std::string> paths = {"top.a.m0", "top.a.m1", "top.a.m2", "top.b.m0", "top.b.m1",
                                              "top.b.m2", "top.a.logic", "outside.logic", "fixed"};
      const std::vector<std::string> masters = {"ram", "ram", "ram", "ram", "ram", "ram", "and", "and", "and"};
      d.centroid_groups = detail::centroid_candidates(d, paths, masters, 2000, 2);
      d.validate();
      Options oc;
      oc.max_iterations = 120;
      oc.bins_x = oc.bins_y = 8;
      oc.add_fillers = oc.legalize = oc.detailed_place = false;
      oc.centroid_weight = .25;
      oc.centroid_ramp = 20;
      oc.target_density = 1;
      auto plain = d, dragged = d;
      place(plain, oc);
      const auto before = detail::centroid_summary(plain);
      check(before.size() == 2);
      oc.interact_interval = 10;
      oc.interact_centroids = true;
      int seen = 0;
      oc.interact = [&](const InteractView &v) {
        InteractUpdate u;
        if (v.centroids) {
          check(v.centroids->size() == 2);
          ++seen;
        }
        if (v.iteration == 20) {
          u.centroid_targets = {{0, 100, 900}};  // drag group 0 to the top-left corner
        }
        return u;
      };
      place(dragged, oc);
      const auto after = detail::centroid_summary(dragged);
      check(seen > 5);
      check(dragged.centroid_targets.empty());  // the caller's design is not rewritten
      const double d_plain = std::hypot(before[0].x - 100, before[0].y - 900);
      const double d_drag = std::hypot(after[0].x - 100, after[0].y - 900);
      std::cout << "centroid 0 distance to target: plain " << d_plain << " dragged " << d_drag << "\n";
      check(d_drag < d_plain - 10);  // the group's members follow its dragged target
    }
    std::cout << "interact smoke ok\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
