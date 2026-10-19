// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok) {
    throw std::runtime_error("RUDY assertion failed at line " +
                             std::to_string(at.line()));
  }
}
void near(double a, double b,
          std::source_location at = std::source_location::current()) {
  check(std::abs(a - b) < 1e-12, at);
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
  d.layout = {10, 20, 50, 44, 1, 1};
  d.node_x = {8, 28, 16, 20, 12, 48, 15, 35, 10, 50};
  d.node_y = {18, 40, 25, 33, 22, 42, 30, 30, 20, 44};
  d.node_size_x.resize(10);
  d.node_size_y.resize(10);
  d.pin2node = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  d.pin_offset_x.resize(10);
  d.pin_offset_y.resize(10);
  d.flat_net2pin = d.pin2node;
  d.flat_net2pin_start = {0, 4, 6, 8, 10};
  d.net_weights = {1.5, 2, 1, 100};
  d.net_mask = {1, 1, 1, 0};
  d.build_maps();
  return d;
}
} // namespace
int main() {
  try {
    auto d = design();
    const auto before = d;
    for (int nx = 1; nx <= 6; ++nx) {
      for (int ny = 1; ny <= 5; ++ny) {
        RudyOptions o;
        o.bins_x = nx;
        o.bins_y = ny;
        o.horizontal_capacity = .025;
        o.vertical_capacity = .05;
        const auto r = rudy(d, o);
        check(r.measured_nets == 2 && r.degenerate_nets == 1);
        double expected_overflow = 0, maximum = 0;
        const double dx = 40. / nx, dy = 24. / ny;
        // Brute force every bin, using independently known net bboxes/weights.
        for (int x = 0; x < nx; ++x) {
          for (int y = 0; y < ny; ++y) {
            double h = 0, v = 0;
            for (int net = 0; net < 2; ++net) {
              const double xl = net ? 12 : 8, yl = net ? 22 : 18;
              const double xh = net ? 48 : 28, yh = net ? 42 : 40;
              const double wt = net ? 2 : 1.5 * 1.0828;
              const double overlap =
                  std::max(0., std::min(xh, 10 + (x + 1) * dx) -
                                   std::max(xl, 10 + x * dx)) *
                  std::max(0., std::min(yh, 20 + (y + 1) * dy) -
                                   std::max(yl, 20 + y * dy));
              h += wt * overlap / (yh - yl) / (dx * dy * o.horizontal_capacity);
              v += wt * overlap / (xh - xl) / (dx * dy * o.vertical_capacity);
            }
            near(h, r.horizontal[x * ny + y]);
            near(v, r.vertical[x * ny + y]);
            const double u = std::max(h, v);
            maximum = std::max(maximum, u);
            expected_overflow += std::max(0., u - 1) / (nx * ny);
          }
        }
        near(maximum, r.max_utilization);
        near(expected_overflow, r.overflow);
      }
    }
    check(d.node_x == before.node_x && d.node_y == before.node_y &&
          d.net_weights == before.net_weights);
    RudyOptions o;
    o.horizontal_capacity = 0;
    rejects([&] { rudy(d, o); });
    o = {};
    o.bins_x = 0;
    rejects([&] { rudy(d, o); });
    o = {};
    o.vertical_capacity = std::numeric_limits<double>::infinity();
    rejects([&] { rudy(d, o); });
    auto huge = d;
    huge.node_x[0] = std::numeric_limits<double>::max();
    huge.node_x[1] = -std::numeric_limits<double>::max();
    rejects([&] { rudy(huge); });
    auto blocked = design();
    blocked.node_macro.assign(blocked.node_x.size(), 0);
    blocked.node_macro[0] = 1;
    blocked.node_size_x[0] = 20;
    blocked.node_size_y[0] = 20;
    blocked.node_blocked_capacity_h.assign(blocked.node_x.size(), 0);
    blocked.node_blocked_capacity_v.assign(blocked.node_x.size(), 0);
    blocked.node_blocked_capacity_h[0] = .01;
    blocked.node_blocked_capacity_v[0] = .02;
    RudyOptions access{2, 2, .025, .05};
    const auto unobstructed = rudy(blocked, access);
    access.macro_blockages = true;
    const auto obstructed = rudy(blocked, access);
    check(obstructed.blocked_fraction[0] > 0);
    check(obstructed.unroutable_nets > 0 ||
          obstructed.max_utilization > unobstructed.max_utilization);
    // This macro's signal pin lies outside the core; move its access point in.
    blocked.pin_offset_x[0] = 4;
    blocked.pin_offset_y[0] = 4;
    access.pin_access = true;
    const auto pins = rudy(blocked, access);
    check(pins.pin_utilization[0] > 0);
    blocked.node_blocked_capacity_h[0] = -1;
    rejects([&] { rudy(blocked, access); });
    // A horizontal net must consume tracks and go around an opaque macro,
    // even when that macro has no LEF obstruction-layer metadata.
    Design corridor;
    corridor.layout = {0, 0, 100, 100, 1, 1};
    corridor.node_x = {5, 95, 40};
    corridor.node_y = {50, 50, 40};
    corridor.node_size_x = {0, 0, 20};
    corridor.node_size_y = {0, 0, 20};
    corridor.node_macro = {0, 0, 1};
    corridor.pin2node = {0, 1};
    corridor.pin_offset_x = corridor.pin_offset_y = {0, 0};
    corridor.flat_net2pin = {0, 1};
    corridor.flat_net2pin_start = {0, 2};
    corridor.net_weights = {1};
    corridor.build_maps();
    RudyOptions hard{10, 10, .1, .1, true, false};
    const auto wide = rudy(corridor, hard);
    check(wide.measured_nets == 1 && wide.degenerate_nets == 0);
    check(wide.unroutable_nets == 0 && wide.detour_wirelength > 0);
    near(wide.blocked_fraction[4 * 10 + 4], 1);
    near(wide.horizontal[4 * 10 + 4], 0);
    near(wide.vertical[4 * 10 + 4], 0);
    double hwire = 0, vwire = 0;
    for (int k = 0; k < 100; ++k) {
      hwire += wide.horizontal[k] * 100 * .1 * (1 - wide.blocked_fraction[k]);
      vwire += wide.vertical[k] * 100 * .1 * (1 - wide.blocked_fraction[k]);
    }
    near(hwire, 90);
    near(vwire, wide.detour_wirelength);
    // Halos reserve placement whitespace, not opaque routing area. Both the
    // public physical design and the internal padded GP view must agree.
    auto halo = corridor;
    halo.node_halo_x = halo.node_halo_y = {0, 0, 5};
    auto padded = halo;
    padded.node_pad_x = padded.node_pad_y = {0, 0, 5};
    padded.node_x[2] -= 5;
    padded.node_y[2] -= 5;
    padded.node_size_x[2] += 10;
    padded.node_size_y[2] += 10;
    for (const auto &view : {halo, padded}) {
      const auto measured = rudy(view, hard);
      check(measured.blocked_fraction == wide.blocked_fraction);
      check(measured.horizontal == wide.horizontal);
      check(measured.vertical == wide.vertical);
    }
    corridor.node_y[2] = 10;
    corridor.node_size_y[2] = 80;
    const auto narrow = rudy(corridor, hard);
    check(narrow.unroutable_nets == 0);
    check(narrow.detour_wirelength > wide.detour_wirelength);
    check(narrow.max_utilization > wide.max_utilization);
    corridor.node_y[2] = 0;
    corridor.node_size_y[2] = 100;
    const auto wall = rudy(corridor, hard);
    check(wall.unroutable_nets == 1 && wall.overflow > 0);
    // Sub-bin disconnected whitespace is not capacity for a net on this side.
    auto pocket = corridor;
    pocket.node_x[1] = 35;
    auto one_bin = hard;
    one_bin.bins_x = one_bin.bins_y = 1;
    const auto local = rudy(pocket, one_bin);
    check(local.unroutable_nets == 0);
    near(local.horizontal[0], 30 / (40 * 100 * .1));
    // A wall assembled from two macros must also disconnect the free space.
    corridor.node_size_y[2] = 50;
    corridor.node_x.push_back(40);
    corridor.node_y.push_back(50);
    corridor.node_size_x.push_back(20);
    corridor.node_size_y.push_back(50);
    corridor.node_macro.push_back(1);
    corridor.build_maps();
    check(rudy(corridor, hard).unroutable_nets == 1);
    // Overlapping footprints remove their union, not twice their shared area.
    corridor.node_y[3] = 25;
    corridor.node_size_y[3] = 50;
    const auto overlap = rudy(corridor, hard);
    double blocked_area = 0;
    for (double f : overlap.blocked_fraction)
      blocked_area += f * 100;
    near(blocked_area, 20 * 75);
    check(overlap.unroutable_nets == 0);
    Design empty;
    empty.layout = {0, 0, 10, 10, 1, 1};
    empty.build_maps();
    near(rudy(empty).overflow, 0);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
