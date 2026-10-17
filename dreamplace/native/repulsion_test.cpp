// BSD-3-Clause. See LICENSE.
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "dreamplace/native/placement.hpp"
using namespace DreamPlace::native;
void check(bool v) {
  if (!v) {
    throw std::runtime_error("repulsion assertion");
  }
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 1413, 1413, 1, 1};
    d.num_movable = 2;
    d.node_x = {10, 12};
    d.node_y = {10, 13};
    d.node_size_x = {7, 7};
    d.node_size_y = {11, 11};
    d.node_macro = {1, 1};
    d.build_maps();
    d.validate();
    check(detail::automatic_bins(d) == std::pair{64, 64});
    auto small = d;
    small.layout.xh = small.layout.yh = 32;
    check(detail::automatic_bins(small) == std::pair{32, 32});
    const auto g = detail::macro_overlap(d, 3);
    check(detail::macro_overlap_counts(d).first == 1);
    check(std::abs(g.x[0] + g.x[1]) < 1e-10);
    for (int i = 0; i < 2; ++i) {
      for (bool y : {false, true}) {
        auto lo = d, hi = d;
        auto &a = y ? lo.node_y : lo.node_x;
        auto &b = y ? hi.node_y : hi.node_x;
        a[i] -= 1e-5;
        b[i] += 1e-5;
        const double fd = (detail::macro_overlap(hi, 3).value -
                           detail::macro_overlap(lo, 3).value) /
                          2e-5;
        check(std::abs(fd - (y ? g.y[i] : g.x[i])) < 1e-7);
      }
    }
    d.node_x[1] = d.node_x[0];
    d.node_y[1] = d.node_y[0];
    check(detail::macro_overlap_counts(d).second == 1);
    const auto tie = detail::macro_overlap(d, 1);
    check(tie.x[0] != 0 && tie.x[0] == -tie.x[1]);
    Options o;
    o.legalize = false;
    o.detailed_place = false;
    o.add_fillers = false;
    o.max_iterations = 60;
    o.automatic_bins = true;
    o.exact_macro_density = true;
    o.macro_jitter = .01;
    o.macro_overlap_weight = 50;
    o.adaptive_density = true;
    // A small layout keeps the direct DCT fixture under the test time budget.
    d.layout.xh = d.layout.yh = 32;
    auto a = d, b = d;
    const auto r = place(a, o);
    place(b, o);
    check(a.node_x == b.node_x && a.node_y == b.node_y);
    check(a.node_x != d.node_x && !r.history.empty());
    d.num_movable = 1;
    auto c = d;
    place(c, o);
    check(c.node_x[1] == d.node_x[1] && c.node_y[1] == d.node_y[1]);
    Design tile;
    tile.layout = {0, 0, 20, 10, 1, 1};
    tile.num_movable = 2;
    tile.node_x = {0, 10};
    tile.node_y = {0, 0};
    tile.node_size_x = {10, 10};
    tile.node_size_y = {10, 10};
    tile.node_macro = {1, 1};
    tile.build_maps();
    Options measure;
    measure.max_iterations = 0;
    measure.legalize = false;
    measure.detailed_place = false;
    measure.target_density = .5;
    auto raw = tile;
    check(std::abs(place(raw, measure).overflow - .5) < 1e-12);
    measure.scale_macro_density = true;
    check(place(tile, measure).overflow == 0);
    Design congested;
    congested.layout = {0, 0, 10, 10, .01, .01};
    congested.num_movable = 2;
    congested.node_x = {1, 9};
    congested.node_y = {1, 9};
    congested.node_size_x = {.1, .1};
    congested.node_size_y = {.1, .1};
    congested.node_inflatable = {1, 1};
    congested.pin2node = {0, 1};
    congested.pin_offset_x = {0, 0};
    congested.pin_offset_y = {0, 0};
    congested.flat_net2pin = {0, 1};
    congested.flat_net2pin_start = {0, 2};
    congested.build_maps();
    Options route;
    route.legalize = false;
    route.detailed_place = false;
    route.add_fillers = false;
    route.max_iterations = 40;
    route.routability = route.improved_routability = true;
    route.route_capacity_h = route.route_capacity_v = 1e-6;
    const auto routing = place(congested, route);
    check(routing.history.size() == 2);
    check(routing.history[0].inflation_area_ratio ==
          1); // wait for measured overflow
    check(routing.history[1].inflation_area_ratio > 1 &&
          routing.history[1].inflation_area_ratio <= 1.05 + 1e-12);
    auto with_fillers = congested;
    with_fillers.node_x.push_back(5);
    with_fillers.node_y.push_back(5);
    with_fillers.node_size_x.push_back(0);
    with_fillers.node_size_y.push_back(0);
    with_fillers.node_inflatable.push_back(0);
    with_fillers.build_maps(); // zero-area fixed IO must retain zero charge
    route.add_fillers = true;
    route.filler_width = .5;
    const auto compensated = place(with_fillers, route);
    check(compensated.history.size() == 2);
    check(compensated.history[1].model_density < route.target_density);
    check(compensated.history[1].filler_charge_scale <
          compensated.history[0].filler_charge_scale);
    auto legacy_fillers = congested;
    route.routability_charge_compensation = false;
    const auto legacy_charge = place(legacy_fillers, route);
    check(legacy_charge.history.size() == 2 &&
          legacy_charge.history.back().filler_charge_scale == 1);
    // Disabled controls preserve the old solver exactly, including RNG
    // sequence.
    auto plain = d, masked = d;
    plain.node_macro.clear();
    Options old;
    old.max_iterations = 10;
    old.legalize = false;
    old.detailed_place = false;
    place(plain, old);
    place(masked, old);
    check(plain.node_x == masked.node_x && plain.node_y == masked.node_y);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
