// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok) {
  if (!ok) {
    throw std::runtime_error("RUDY assertion failed");
  }
}
void near(double a, double b) { check(std::abs(a - b) < 1e-12); }
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
    check(obstructed.horizontal[0] > unobstructed.horizontal[0]);
    check(obstructed.vertical[0] > unobstructed.vertical[0]);
    // This macro's signal pin lies outside the core; move its access point in.
    blocked.pin_offset_x[0] = 4;
    blocked.pin_offset_y[0] = 4;
    access.pin_access = true;
    const auto pins = rudy(blocked, access);
    check(pins.pin_utilization[0] > 0);
    blocked.node_blocked_capacity_h[0] = -1;
    rejects([&] { rudy(blocked, access); });
    {
      // A pin one ulp inside the right/top edge can round to bin index `bins`
      // when divided by the bin pitch; it must still land in the last bin.
      Design edge;
      edge.layout = {0, 0, 1000, 1000, 1, 1};
      edge.num_movable = 0;
      const double inside = std::nextafter(1000., 0.);
      edge.node_x = {inside, 10};
      edge.node_y = {inside, 10};
      edge.node_size_x = {0, 0};
      edge.node_size_y = {0, 0};
      edge.node_macro = {1, 0};
      edge.pin2node = {0, 1};
      edge.pin_offset_x = {0, 0};
      edge.pin_offset_y = {0, 0};
      edge.flat_net2pin = {0, 1};
      edge.flat_net2pin_start = {0, 2};
      edge.net_weights = {1};
      edge.net_mask = {1};
      edge.build_maps();
      RudyOptions o3{3, 3, 1, 1};
      o3.pin_access = true;
      const auto r = rudy(edge, o3);
      check(r.pin_utilization.size() == 9);
      check(r.pin_utilization[2 * 3 + 2] > 0);
      double elsewhere = 0;
      for (size_t i = 0; i < r.pin_utilization.size(); ++i)
        if (i != 8)
          elsewhere += r.pin_utilization[i];
      near(elsewhere, 0);
    }
    Design empty;
    empty.layout = {0, 0, 10, 10, 1, 1};
    empty.build_maps();
    near(rudy(empty).overflow, 0);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
