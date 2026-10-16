// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("GP hook assertion failed");
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 80, 40, 1, 1};
  d.num_movable = 3;
  d.node_x = {10, 20, 30, 0, 80};
  d.node_y = {10, 20, 10, 20, 20};
  d.node_size_x = {2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 0, 0};
  d.pin2node = {3, 0, 0, 1, 1, 2, 2, 4};
  d.pin_offset_x.assign(8, 0);
  d.pin_offset_y.assign(8, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7};
  d.flat_net2pin_start = {0, 2, 4, 6, 8};
  d.build_maps();
  return d;
}
int main() {
  try {
    auto full = fixture(), split = full;
    Options o;
    o.max_iterations = 40;
    o.schedule_iterations = 40;
    o.legalize = o.detailed_place = false;
    o.bins_x = o.bins_y = 8;
    o.stop_overflow = 0;
    place(full, o);
    GpState state;
    o.max_iterations = 20;
    place(split, o, &state);
    check(state.iteration == 20);
    place(split, o, &state);
    check(state.iteration == 40);
    check(full.node_x == split.node_x && full.node_y == split.node_y);
    const auto old = split;
    const auto history = state;
    split.node_size_x[0] = 3;
    bool bad = false;
    try {
      place(split, o, &state);
    } catch (const std::invalid_argument &) {
      bad = true;
    }
    check(bad && state.x == history.x && state.iteration == history.iteration);
    split = old;
    auto boxed = fixture();
    boxed.node_boxes.assign(5, {0, 0, 78, 39});
    boxed.node_boxes[0] = {4, 6, 4, 6};
    int snapshots = 0;
    o.observer_interval = 1;
    o.observer = [&](const Design &d, PlacementStage, int) {
      ++snapshots;
      check(d.node_x[0] == 4 && d.node_y[0] == 6);
    };
    place(boxed, o);
    check(snapshots > 2);
    auto order = fixture();
    order.node_x[0] = 50;
    order.node_x[1] = 10;
    order.order_edges.push_back({0, 1, Axis::x, 0, {false, 1}, "left"});
    o.observer = {};
    o.order_weight = 1000;
    o.max_iterations = 100;
    o.add_fillers = false;
    place(order, o);
    check(order.node_x[0] + 2 <= order.node_x[1] + .1);
    auto route = fixture();
    route.node_inflatable = {1, 1, 1, 0, 0};
    auto sizes = route.node_size_x;
    o.order_weight = 0;
    o.routability = true;
    o.route_capacity_h = o.route_capacity_v = .01;
    place(route, o);
    check(route.node_size_x == sizes);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
