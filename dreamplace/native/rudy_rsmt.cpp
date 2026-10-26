// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy_tree.hpp"
#include "rudy/src/parameters.h"
#include "thirdparty/flute/flute.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
namespace DreamPlace::native {
namespace {
template <typename T> inline DEFINE_NET_WIRING_DISTRIBUTION_MAP_WEIGHT;
}
namespace detail {
std::vector<RouteEdge>
route_tree(const std::vector<std::pair<double, double>> &pins,
           const RudyOptions &o) {
  if (pins.size() > size_t(o.degree))
    return {};
  std::vector<std::pair<int, int>> points;
  for (auto [x, y] : pins) {
    if (!std::isfinite(x) || !std::isfinite(y) ||
        std::abs(x) > std::numeric_limits<int>::max() / 4. ||
        std::abs(y) > std::numeric_limits<int>::max() / 4.)
      throw std::invalid_argument("coordinates exceed FLUTE integer range");
    points.emplace_back(std::llround(x), std::llround(y));
  }
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  std::vector<RouteEdge> edges;
  auto edge = [&](double x, double y, double u, double v) {
    if (x != u || y != v)
      edges.push_back(
          {std::min(x, u), std::min(y, v), std::max(x, u), std::max(y, v)});
  };
  // Preserve fractional coordinates and the existing box model for two pins.
  if (pins.size() == 2) {
    edge(pins[0].first, pins[0].second, pins[1].first, pins[1].second);
    return edges;
  }
  const int degree = int(points.size());
  if (degree < 2)
    return edges;
  if (degree == 2) {
    edge(points[0].first, points[0].second, points[1].first, points[1].second);
    return edges;
  }
  std::vector<int> xs, ys;
  for (auto [x, y] : points) {
    xs.push_back(x);
    ys.push_back(y);
  }
  if (degree == 3) {
    auto a = xs, b = ys;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    for (int i = 0; i < 3; ++i)
      edge(xs[i], ys[i], a[1], b[1]);
    return edges;
  }
  const double span = double(*std::max_element(xs.begin(), xs.end())) -
                      *std::min_element(xs.begin(), xs.end()) +
                      double(*std::max_element(ys.begin(), ys.end())) -
                      *std::min_element(ys.begin(), ys.end());
  if (span * degree > std::numeric_limits<int>::max())
    throw std::invalid_argument("net wirelength exceeds FLUTE integer range");
  if (!o.initialize)
    throw std::runtime_error(
        "RSMT-RUDY requires a FLUTE lookup-table initializer");
  o.initialize();
  auto tree = flute::flute(degree, xs.data(), ys.data(), o.accuracy);
  std::unique_ptr<flute::Branch, decltype(&std::free)> branches(tree.branch,
                                                                &std::free);
  for (int i = 0; i < 2 * degree - 2; ++i) {
    const auto &a = tree.branch[i], &b = tree.branch[a.n];
    edge(a.x, a.y, b.x, b.y);
  }
  return edges;
}
} // namespace detail
RudyResult rsmt_rudy(const Design &d, const RudyOptions &o) {
  RudyResult r;
  const int count = o.bins_x * o.bins_y;
  r.horizontal.resize(count);
  r.vertical.resize(count);
  r.pin_utilization.resize(count);
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / o.bins_x, dy = (l.yh - l.yl) / o.bins_y;
  const double gx = std::min(dx, l.site_width), gy = std::min(dy, l.row_height);
  auto bx = [&](double x) {
    return detail::clamped_bin(x, l.xl, dx, o.bins_x);
  };
  auto by = [&](double y) {
    return detail::clamped_bin(y, l.yl, dy, o.bins_y);
  };
  if (o.pin_access)
    for (int p = 0; p < int(d.pin2node.size()); ++p) {
      int n = d.pin2node[p];
      if (d.node_macro.empty() || !d.node_macro[n] || !d.net_mask[d.pin2net[p]])
        continue;
      const double x = std::clamp(d.node_x[n] + d.pin_offset_x[p], l.xl, l.xh),
                   y = std::clamp(d.node_y[n] + d.pin_offset_y[p], l.yl, l.yh);
      r.pin_utilization[bx(x) * o.bins_y + by(y)] +=
          1 /
          ((dx + dy) * std::sqrt(o.horizontal_capacity * o.vertical_capacity));
    }
  for (int n = 0; n < int(d.net_weights.size()); ++n) {
    if (!d.net_mask[n] || d.net_weights[n] == 0)
      continue;
    std::vector<std::pair<double, double>> pins;
    detail::RouteEdge box{l.xh, l.yh, l.xl, l.yl};
    for (int k = d.flat_net2pin_start[n]; k < d.flat_net2pin_start[n + 1];
         ++k) {
      int p = d.flat_net2pin[k], node = d.pin2node[p];
      const double px = d.node_x[node] + d.pin_offset_x[p],
                   py = d.node_y[node] + d.pin_offset_y[p];
      if (!std::isfinite(px) || !std::isfinite(py))
        throw std::invalid_argument("RUDY pin position overflow");
      double x = std::clamp(px, l.xl, l.xh), y = std::clamp(py, l.yl, l.yh);
      pins.emplace_back(x, y);
      box.xl = std::min(box.xl, x);
      box.xh = std::max(box.xh, x);
      box.yl = std::min(box.yl, y);
      box.yh = std::max(box.yh, y);
    }
    auto edges = detail::route_tree(pins, o);
    double weight = d.net_weights[n];
    if (pins.size() > size_t(o.degree)) {
      edges.push_back(box);
      weight *= netWiringDistributionMapWeight<double>(pins.size());
    }
    if (edges.empty()) {
      ++r.degenerate_nets;
      continue;
    }
    ++r.measured_nets;
    for (auto e : edges) {
      double w = e.xh - e.xl, h = e.yh - e.yl;
      if (w == 0) {
        e.xl = std::max(l.xl, e.xl - gx / 2);
        e.xh = std::min(l.xh, e.xh + gx / 2);
      }
      if (h == 0) {
        e.yl = std::max(l.yl, e.yl - gy / 2);
        e.yh = std::min(l.yh, e.yh + gy / 2);
      }
      double area = (e.xh - e.xl) * (e.yh - e.yl);
      if (area <= 0)
        continue;
      for (int x = bx(e.xl); x <= bx(e.xh); ++x)
        for (int y = by(e.yl); y <= by(e.yh); ++y) {
          const double overlap =
              std::max(0., std::min(e.xh, l.xl + (x + 1) * dx) -
                               std::max(e.xl, l.xl + x * dx)) *
              std::max(0., std::min(e.yh, l.yl + (y + 1) * dy) -
                               std::max(e.yl, l.yl + y * dy));
          int k = x * o.bins_y + y;
          r.horizontal[k] +=
              weight * w * overlap / area / (dx * dy * o.horizontal_capacity);
          r.vertical[k] +=
              weight * h * overlap / area / (dx * dy * o.vertical_capacity);
        }
    }
  }
  detail::accumulate_overflow(r);
  if (!std::isfinite(r.overflow))
    throw std::invalid_argument("RSMT-RUDY overflow");
  return r;
}
} // namespace DreamPlace::native
