// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "rudy/src/parameters.h"

namespace DreamPlace::native {
namespace {
template <typename T> inline DEFINE_NET_WIRING_DISTRIBUTION_MAP_WEIGHT;
struct Rect {
  double xl, yl, xh, yh;
  double area() const { return std::max(0., xh - xl) * std::max(0., yh - yl); }
};
Rect intersect(Rect a, Rect b) {
  return {std::max(a.xl, b.xl), std::max(a.yl, b.yl), std::min(a.xh, b.xh),
          std::min(a.yh, b.yh)};
}
// Non-overlapping free rectangles: subtract the UNION of macro footprints,
// including overlapping macros during GP, without subtracting capacity twice.
void subtract(std::vector<Rect> &free, Rect obstacle) {
  std::vector<Rect> next;
  for (auto a : free) {
    auto b = intersect(a, obstacle);
    if (b.area() <= 0) {
      next.push_back(a);
      continue;
    }
    for (auto piece :
         {Rect{a.xl, a.yl, b.xl, a.yh}, Rect{b.xh, a.yl, a.xh, a.yh},
          Rect{b.xl, a.yl, b.xh, b.yl}, Rect{b.xl, b.yh, b.xh, a.yh}})
      if (piece.area() > 0)
        next.push_back(piece);
  }
  free = std::move(next);
}
} // namespace

// Obstacle-aware RUDY remains a demand estimator, not a detailed router. It
// spreads wire demand only over the exact free area of each bin. A macro that
// cuts the whole projected corridor forces that corridor around a macro edge;
// the extra Manhattan excursion contributes demand too. Full-core barriers and
// boxes with no free area are reported, never given artificial routing
// capacity.
RudyResult obstacle_rudy(const Design &d, const RudyOptions &o) {
  const auto &l = d.layout;
  const Rect core{l.xl, l.yl, l.xh, l.yh};
  const double dx = (l.xh - l.xl) / o.bins_x, dy = (l.yh - l.yl) / o.bins_y;
  const double gx = std::min(dx, l.site_width), gy = std::min(dy, l.row_height);
  const int count = o.bins_x * o.bins_y;
  auto bx = [&](double x) {
    return int(
        std::clamp(std::floor((x - l.xl) / dx), 0., double(o.bins_x - 1)));
  };
  auto by = [&](double y) {
    return int(
        std::clamp(std::floor((y - l.yl) / dy), 0., double(o.bins_y - 1)));
  };
  auto bins = [&](Rect a, const auto &visit) {
    for (int x = bx(a.xl); x <= bx(a.xh); ++x)
      for (int y = by(a.yl); y <= by(a.yh); ++y)
        visit(x * o.bins_y + y);
  };
  std::vector<Rect> obstacles;
  std::vector<int> macro(d.node_x.size(), -1);
  std::vector<std::vector<int>> blockers(count);
  std::vector<std::vector<Rect>> free(count);
  for (int x = 0; x < o.bins_x; ++x)
    for (int y = 0; y < o.bins_y; ++y)
      free[x * o.bins_y + y].push_back({l.xl + x * dx, l.yl + y * dy,
                                        l.xl + (x + 1) * dx,
                                        l.yl + (y + 1) * dy});
  for (int i = 0; i < int(d.node_macro.size()); ++i) {
    if (!d.node_macro[i])
      continue;
    // GP occupancy envelopes include halos; those halos are routing whitespace,
    // not SRAM body. Recover the physical footprint in the native padded solve.
    const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
    const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
    Rect b = intersect(core, {d.node_x[i] + px, d.node_y[i] + py,
                              d.node_x[i] + d.node_size_x[i] - px,
                              d.node_y[i] + d.node_size_y[i] - py});
    if (b.area() <= 0)
      continue;
    macro[i] = int(obstacles.size());
    obstacles.push_back(b);
    bins(b, [&](int k) {
      blockers[k].push_back(macro[i]);
      subtract(free[k], b);
    });
  }
  RudyResult r;
  r.horizontal.resize(count);
  r.vertical.resize(count);
  r.pin_utilization.resize(count);
  r.blocked_fraction.resize(count);
  std::vector<double> available(count);
  for (int k = 0; k < count; ++k) {
    for (auto a : free[k])
      available[k] += a.area();
    r.blocked_fraction[k] = std::clamp(1 - available[k] / (dx * dy), 0., 1.);
  }
  // Exact free-space connectivity, including sub-bin channels. Bins alone can
  // miss a thin full-height wall or incorrectly close a narrow legal channel.
  std::vector<int> starts(count + 1);
  for (int k = 0; k < count; ++k)
    starts[k + 1] = starts[k] + int(free[k].size());
  std::vector<int> parent(starts.back());
  std::iota(parent.begin(), parent.end(), 0);
  auto root = [&](int i) {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };
  const double epsilon = 1e-12 * std::max(dx, dy);
  auto touches = [&](Rect a, Rect b) {
    return ((std::abs(a.xh - b.xl) <= epsilon ||
             std::abs(b.xh - a.xl) <= epsilon) &&
            std::min(a.yh, b.yh) - std::max(a.yl, b.yl) > epsilon) ||
           ((std::abs(a.yh - b.yl) <= epsilon ||
             std::abs(b.yh - a.yl) <= epsilon) &&
            std::min(a.xh, b.xh) - std::max(a.xl, b.xl) > epsilon);
  };
  auto connect = [&](int a, int b) {
    for (int i = 0; i < int(free[a].size()); ++i)
      for (int j = (a == b ? i + 1 : 0); j < int(free[b].size()); ++j)
        if (touches(free[a][i], free[b][j]))
          parent[root(starts[a] + i)] = root(starts[b] + j);
  };
  for (int x = 0; x < o.bins_x; ++x)
    for (int y = 0; y < o.bins_y; ++y) {
      const int k = x * o.bins_y + y;
      connect(k, k);
      if (x + 1 < o.bins_x)
        connect(k, k + o.bins_y);
      if (y + 1 < o.bins_y)
        connect(k, k + 1);
    }
  std::vector<std::array<double, 2>> pins(d.pin2node.size());
  std::vector<int> pin_component(d.pin2node.size(), -1);
  for (int p = 0; p < int(d.pin2node.size()); ++p) {
    const int node = d.pin2node[p];
    auto &q = pins[p];
    const double x = d.node_x[node] + d.pin_offset_x[p];
    const double y = d.node_y[node] + d.pin_offset_y[p];
    if (!std::isfinite(x) || !std::isfinite(y))
      throw std::invalid_argument("RUDY pin position overflow");
    q = {std::clamp(x, l.xl, l.xh), std::clamp(y, l.yl, l.yh)};
    if (macro[node] >= 0) {
      const auto a = obstacles[macro[node]];
      // A terminal belongs at an access edge, not inside an opaque SRAM/halo.
      if (q[0] >= a.xl && q[0] <= a.xh && q[1] >= a.yl && q[1] <= a.yh) {
        std::array<double, 4> dist{q[0] - a.xl, a.xh - q[0], q[1] - a.yl,
                                   a.yh - q[1]};
        if (a.xl <= l.xl)
          dist[0] = std::numeric_limits<double>::infinity();
        if (a.xh >= l.xh)
          dist[1] = std::numeric_limits<double>::infinity();
        if (a.yl <= l.yl)
          dist[2] = std::numeric_limits<double>::infinity();
        if (a.yh >= l.yh)
          dist[3] = std::numeric_limits<double>::infinity();
        const int side =
            int(std::min_element(dist.begin(), dist.end()) - dist.begin());
        if (std::isfinite(dist[side])) {
          if (side == 0)
            q[0] = std::nextafter(a.xl, l.xl);
          if (side == 1)
            q[0] = std::nextafter(a.xh, l.xh);
          if (side == 2)
            q[1] = std::nextafter(a.yl, l.yl);
          if (side == 3)
            q[1] = std::nextafter(a.yh, l.yh);
        }
      }
    }
    const int k = bx(q[0]) * o.bins_y + by(q[1]);
    for (int i = 0; i < int(free[k].size()); ++i) {
      const auto a = free[k][i];
      if (q[0] >= a.xl && q[0] <= a.xh && q[1] >= a.yl && q[1] <= a.yh) {
        pin_component[p] = root(starts[k] + i);
        break;
      }
    }
    const bool active =
        d.net_mask[d.pin2net[p]] && d.net_weights[d.pin2net[p]] != 0;
    if (active && pin_component[p] < 0) {
      ++r.blocked_pins;
      // This is an invalid access point, not capacity inside the macro. Keep
      // the violation local so GP can inflate/expel overlapping standard cells.
      r.pin_utilization[k] = std::max(r.pin_utilization[k], 2.);
    }
    if (macro[node] >= 0 && o.pin_access && active) {
      if (available[k] > 0)
        r.pin_utilization[k] +=
            1 / ((dx + dy) *
                 std::sqrt(o.horizontal_capacity * o.vertical_capacity) *
                 available[k] / (dx * dy));
    }
  }
  std::vector<std::uint64_t> seen(obstacles.size());
  std::uint64_t epoch = 0;
  std::vector<double> overlap(count), component_area(count);
  for (int n = 0; n < int(d.net_weights.size()); ++n) {
    if (!d.net_mask[n] || d.net_weights[n] == 0)
      continue;
    Rect box{l.xh, l.yh, l.xl, l.yl};
    int component = -1;
    bool disconnected = false;
    for (int k = d.flat_net2pin_start[n]; k < d.flat_net2pin_start[n + 1];
         ++k) {
      const int pin = d.flat_net2pin[k];
      const auto p = pins[pin];
      if (pin_component[pin] < 0 ||
          (component >= 0 && component != pin_component[pin]))
        disconnected = true;
      component = pin_component[pin];
      box.xl = std::min(box.xl, p[0]);
      box.xh = std::max(box.xh, p[0]);
      box.yl = std::min(box.yl, p[1]);
      box.yh = std::max(box.yh, p[1]);
    }
    const double width = box.xh - box.xl, height = box.yh - box.yl;
    if (width == 0 && height == 0) {
      ++r.degenerate_nets;
      continue;
    }
    ++r.measured_nets;
    if (disconnected) {
      ++r.unroutable_nets;
      continue;
    }
    // A straight net still consumes tracks. Give its zero-width projection a
    // site/row-wide corridor instead of silently dropping all its demand.
    if (width == 0) {
      box.xl = std::max(l.xl, box.xl - gx / 2);
      box.xh = std::min(l.xh, box.xh + gx / 2);
    }
    if (height == 0) {
      box.yl = std::max(l.yl, box.yl - gy / 2);
      box.yh = std::min(l.yh, box.yh + gy / 2);
    }
    const Rect initial = box;
    auto free_area = [&](Rect region) {
      double area = 0;
      bins(region, [&](int k) {
        overlap[k] = component_area[k] = 0;
        for (int i = 0; i < int(free[k].size()); ++i)
          if (root(starts[k] + i) == component) {
            overlap[k] += intersect(free[k][i], region).area();
            component_area[k] += free[k][i].area();
          }
        area += overlap[k];
      });
      return area;
    };

    struct Channel {
      Rect region;
      double horizontal, vertical;
    };
    std::vector<Channel> channels;
    bool impossible = false, changed = true;
    // Every expansion reaches a new macro edge. Iterations are bounded by the
    // obstacle count; the fallback is explicit, not a route through a macro.
    for (size_t step = 0; changed && step <= obstacles.size(); ++step) {
      changed = false;
      ++epoch;
      std::vector<int> candidates;
      bins(box, [&](int k) {
        for (int id : blockers[k])
          if (seen[id] != epoch) {
            seen[id] = epoch;
            candidates.push_back(id);
          }
      });
      for (int id : candidates) {
        const auto a = obstacles[id];
        if (intersect(a, box).area() <= 0)
          continue;
        if (a.yl <= box.yl && a.yh >= box.yh && a.xl > box.xl &&
            a.xh < box.xh) {
          const double up =
              a.yh < l.yh && free_area({a.xl, a.yh, a.xh, l.yh}) > 0
                  ? std::min(l.yh, a.yh + gy) - box.yh
                  : std::numeric_limits<double>::infinity();
          const double down =
              a.yl > l.yl && free_area({a.xl, l.yl, a.xh, a.yl}) > 0
                  ? box.yl - std::max(l.yl, a.yl - gy)
                  : std::numeric_limits<double>::infinity();
          if (!std::isfinite(std::min(up, down))) {
            impossible = true;
            break;
          }
          if (up < down) {
            channels.push_back({{a.xl, a.yh, a.xh, l.yh}, a.xh - a.xl, 0});
            box.yh += up;
          } else {
            channels.push_back({{a.xl, l.yl, a.xh, a.yl}, a.xh - a.xl, 0});
            box.yl -= down;
          }
          changed = true;
        }
        if (a.xl <= box.xl && a.xh >= box.xh && a.yl > box.yl &&
            a.yh < box.yh) {
          const double right =
              a.xh < l.xh && free_area({a.xh, a.yl, l.xh, a.yh}) > 0
                  ? std::min(l.xh, a.xh + gx) - box.xh
                  : std::numeric_limits<double>::infinity();
          const double left =
              a.xl > l.xl && free_area({l.xl, a.yl, a.xl, a.yh}) > 0
                  ? box.xl - std::max(l.xl, a.xl - gx)
                  : std::numeric_limits<double>::infinity();
          if (!std::isfinite(std::min(left, right))) {
            impossible = true;
            break;
          }
          if (right < left) {
            channels.push_back({{a.xh, a.yl, l.xh, a.yh}, 0, a.yh - a.yl});
            box.xh += right;
          } else {
            channels.push_back({{l.xl, a.yl, a.xl, a.yh}, 0, a.yh - a.yl});
            box.xl -= left;
          }
          changed = true;
        }
      }
      if (impossible)
        break;
    }
    if (impossible || free_area(box) <= 0) {
      ++r.unroutable_nets;
      continue;
    }
    const double extra_x =
        std::max(0., box.xh - box.xl - (initial.xh - initial.xl));
    const double extra_y =
        std::max(0., box.yh - box.yl - (initial.yh - initial.yl));
    const double weight =
        d.net_weights[n] *
        netWiringDistributionMapWeight<double>(d.flat_net2pin_start[n + 1] -
                                               d.flat_net2pin_start[n]);
    r.detour_wirelength += weight * 2 * (extra_x + extra_y);
    auto deposit = [&](Rect region, double horizontal, double vertical) {
      const double area = free_area(region);
      if (area <= 0)
        return false;
      bins(region, [&](int k) {
        if (overlap[k] <= 0 || component_area[k] <= 0)
          return;
        // A disconnected pocket in the same bin cannot supply this net with
        // capacity. Use only the free component reachable by its terminals.
        r.horizontal[k] += weight * horizontal * overlap[k] / area /
                           (component_area[k] * o.horizontal_capacity);
        r.vertical[k] += weight * vertical * overlap[k] / area /
                         (component_area[k] * o.vertical_capacity);
      });
      return true;
    };
    double horizontal = width + 2 * extra_x, vertical = height + 2 * extra_y;
    // The crossing portion must use a channel beside the blocking macro. Do
    // not dilute it across the whole expanded box: that would make a narrow
    // passage look cheaper simply because the detour made the box larger.
    for (auto channel : channels) {
      const double h = std::min(horizontal, channel.horizontal);
      const double v = std::min(vertical, channel.vertical);
      if (deposit(channel.region, h, v)) {
        horizontal -= h;
        vertical -= v;
      } else
        impossible = true;
    }
    deposit(box, horizontal, vertical);
    if (impossible)
      ++r.unroutable_nets;
  }
  for (int k = 0; k < count; ++k) {
    const double u =
        std::max({r.horizontal[k], r.vertical[k], r.pin_utilization[k]});
    r.max_utilization = std::max(r.max_utilization, u);
    r.overflow += std::max(0., u - 1) / count;
  }
  // Geometry that cannot escape an obstacle must prevent false convergence.
  r.overflow += double(r.unroutable_nets + r.blocked_pins) / count;
  if (r.unroutable_nets || r.blocked_pins)
    r.max_utilization = std::max(r.max_utilization, 2.);
  if (!std::isfinite(r.overflow) || !std::isfinite(r.max_utilization))
    throw std::invalid_argument("obstacle RUDY overflow");
  return r;
}
} // namespace DreamPlace::native
