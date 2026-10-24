// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "rudy/src/parameters.h"

namespace DreamPlace::native {
RudyResult obstacle_rudy(const Design &, const RudyOptions &);
RudyResult rsmt_rudy(const Design &, const RudyOptions &);
namespace {
template <typename T> inline DEFINE_NET_WIRING_DISTRIBUTION_MAP_WEIGHT;
}
void detail::validate_rudy(const Design &d, const RudyOptions &o) {
  d.validate();
  if (o.bins_x < 1 || o.bins_x > 1024 || o.bins_y < 1 || o.bins_y > 1024 ||
      !std::isfinite(o.horizontal_capacity) || o.horizontal_capacity <= 0 ||
      !std::isfinite(o.vertical_capacity) || o.vertical_capacity <= 0) {
    throw std::invalid_argument("invalid RUDY grid or capacity");
  }
  if (o.degree < 2 || o.degree > 1000 || o.accuracy < 1 || o.accuracy > 10)
    throw std::invalid_argument(
        "invalid FLUTE degree (2..1000) or accuracy (1..10)");
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / o.bins_x, dy = (l.yh - l.yl) / o.bins_y;
  if (!std::isfinite(dx * dy * o.horizontal_capacity) ||
      !std::isfinite(dx * dy * o.vertical_capacity) ||
      dx * dy * o.horizontal_capacity <= 0 ||
      dx * dy * o.vertical_capacity <= 0) {
    throw std::invalid_argument("RUDY bin capacity overflow or underflow");
  }
}
RudyResult rudy(const Design &d, const RudyOptions &o) {
  detail::validate_rudy(d, o);
  RudyResult r;
  r.horizontal.resize(o.bins_x * o.bins_y);
  r.vertical.resize(o.bins_x * o.bins_y);
  r.pin_utilization.resize(o.bins_x * o.bins_y);
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / o.bins_x, dy = (l.yh - l.yl) / o.bins_y;
  if (o.macro_blockages) {
    return obstacle_rudy(d, o);
  }
  if (o.model == RudyModel::rsmt)
    return rsmt_rudy(d, o);
  auto bin = [](double coordinate, double lo, double pitch, int count) {
    // Clip before integer conversion, including far-outside terminal positions.
    return static_cast<int>(
        std::clamp(std::floor((coordinate - lo) / pitch), 0., double(count)));
  };
  std::vector<double> available_h(o.bins_x * o.bins_y, o.horizontal_capacity);
  std::vector<double> available_v(o.bins_x * o.bins_y, o.vertical_capacity);
  if (o.pin_access) {
    for (int p = 0; p < static_cast<int>(d.pin2node.size()); ++p) {
      const int node = d.pin2node[p];
      if (d.node_macro.empty() || !d.node_macro[node] ||
          !d.net_mask[d.pin2net[p]]) {
        continue;
      }
      const double x = d.node_x[node] + d.pin_offset_x[p],
                   y = d.node_y[node] + d.pin_offset_y[p];
      if (x < l.xl || x >= l.xh || y < l.yl || y >= l.yh) {
        continue;
      }
      // bin() may return `count` for a coordinate a ulp inside the edge, where
      // the division rounds up to the grid size; keep the pin in the last bin.
      const int k =
          std::min(bin(x, l.xl, dx, o.bins_x), o.bins_x - 1) * o.bins_y +
          std::min(bin(y, l.yl, dy, o.bins_y), o.bins_y - 1);
      // Access sites scale with tracks crossing bin edges, not bin area.
      r.pin_utilization[k] +=
          1 / ((dx + dy) * std::sqrt(available_h[k] * available_v[k]));
    }
  }
  for (int n = 0; n < static_cast<int>(d.net_weights.size()); ++n) {
    if (!d.net_mask[n] || d.net_weights[n] == 0) {
      continue;
    }
    double xl = std::numeric_limits<double>::infinity(), yl = xl, xh = -xl,
           yh = -xl;
    for (int k = d.flat_net2pin_start[n]; k < d.flat_net2pin_start[n + 1];
         ++k) {
      const int p = d.flat_net2pin[k], node = d.pin2node[p];
      const double x = d.node_x[node] + d.pin_offset_x[p],
                   y = d.node_y[node] + d.pin_offset_y[p];
      if (!std::isfinite(x) || !std::isfinite(y)) {
        throw std::invalid_argument("RUDY pin position overflow");
      }
      xl = std::min(xl, x);
      xh = std::max(xh, x);
      yl = std::min(yl, y);
      yh = std::max(yh, y);
    }
    const double w = xh - xl, h = yh - yl;
    if (!std::isfinite(w) || !std::isfinite(h)) {
      throw std::invalid_argument("RUDY net bbox overflow");
    }
    if (w == 0 || h == 0) {
      ++r.degenerate_nets;
      continue;
    }
    ++r.measured_nets;
    const double weight =
        d.net_weights[n] *
        netWiringDistributionMapWeight<double>(d.flat_net2pin_start[n + 1] -
                                               d.flat_net2pin_start[n]);
    const int x0 = bin(xl, l.xl, dx, o.bins_x),
              x1 = bin(xh, l.xl, dx, o.bins_x);
    const int y0 = bin(yl, l.yl, dy, o.bins_y),
              y1 = bin(yh, l.yl, dy, o.bins_y);
    for (int x = x0; x < std::min(x1 + 1, o.bins_x); ++x) {
      for (int y = y0; y < std::min(y1 + 1, o.bins_y); ++y) {
        const double a = std::max(0., std::min(xh, l.xl + (x + 1) * dx) -
                                          std::max(xl, l.xl + x * dx));
        const double b = std::max(0., std::min(yh, l.yl + (y + 1) * dy) -
                                          std::max(yl, l.yl + y * dy));
        const int i = x * o.bins_y + y;
        r.horizontal[i] +=
            weight * (a / dx) * (b / h) / dy / o.horizontal_capacity;
        r.vertical[i] += weight * (a / w) * (b / dy) / dx / o.vertical_capacity;
        if (!std::isfinite(r.horizontal[i]) || !std::isfinite(r.vertical[i])) {
          throw std::invalid_argument("RUDY demand overflow");
        }
      }
    }
  }
  for (size_t i = 0; i < r.horizontal.size(); ++i) {
    r.horizontal[i] *= o.horizontal_capacity / available_h[i];
    r.vertical[i] *= o.vertical_capacity / available_v[i];
    const double u =
        std::max({r.horizontal[i], r.vertical[i], r.pin_utilization[i]});
    r.max_utilization = std::max(r.max_utilization, u);
    r.overflow += std::max(0., u - 1) / r.horizontal.size();
  }
  if (!std::isfinite(r.overflow)) {
    throw std::invalid_argument("RUDY overflow metric overflow");
  }
  return r;
}
} // namespace DreamPlace::native
