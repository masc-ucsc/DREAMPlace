// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <vector>

namespace DreamPlace::native::detail {
MacroRoutingPressure macro_routing_pressure(const Design &d, int bx, int by,
                                            const std::vector<double> &u,
                                            double range_bins) {
  if (bx < 1 || by < 1 || u.size() != size_t(bx) * by ||
      !std::isfinite(range_bins) || range_bins <= 0 || range_bins > 16)
    throw std::invalid_argument("invalid macro routing pressure map");
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  const double rx = dx * range_bins, ry = dy * range_bins;
  struct Box {
    int node;
    double xl, yl, xh, yh;
  };
  std::vector<Box> macros;
  for (int i = 0; i < int(d.node_macro.size()); ++i) {
    if (!d.node_macro[i])
      continue;
    const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
    const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
    Box b{i, d.node_x[i] + px, d.node_y[i] + py,
          d.node_x[i] + d.node_size_x[i] - px,
          d.node_y[i] + d.node_size_y[i] - py};
    if (b.xh > b.xl && b.yh > b.yl)
      macros.push_back(b);
  }
  MacroRoutingPressure r;
  r.x.resize(d.node_x.size());
  r.y.resize(d.node_y.size());
  auto ix = [&](double x) {
    return int(std::clamp(std::floor((x - l.xl) / dx), 0., double(bx - 1)));
  };
  auto iy = [&](double y) {
    return int(std::clamp(std::floor((y - l.yl) / dy), 0., double(by - 1)));
  };
  auto bins = [&](double xl, double yl, double xh, double yh,
                  const auto &visit) {
    for (int x = ix(xl); x <= ix(xh); ++x)
      for (int y = iy(yl); y <= iy(yh); ++y)
        visit(x, y);
  };
  std::vector<std::vector<int>> grid(bx * by);
  for (int m = 0; m < int(macros.size()); ++m) {
    const auto a = macros[m];
    bins(a.xl, a.yl, a.xh, a.yh,
         [&](int x, int y) { grid[x * by + y].push_back(m); });
  }
  auto candidates = [&](double xl, double yl, double xh, double yh) {
    std::vector<int> ids;
    bins(xl, yl, xh, yh, [&](int x, int y) {
      ids.insert(ids.end(), grid[x * by + y].begin(), grid[x * by + y].end());
    });
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
  };
  // Every visible facing-edge segment is a physical channel. Occlusion is
  // computed below bin resolution: a macro behind another one is not a
  // neighbor.
  auto channel = [&](int first, int second, bool horizontal, double lo,
                     double hi, double begin, double end, double reference) {
    const double reach = horizontal ? rx : ry, gap = hi - lo;
    if (gap < 0 || gap >= reach || end <= begin || reference <= 0)
      return;
    if ((first < 0 || macros[first].node >= d.num_movable) &&
        (second < 0 || macros[second].node >= d.num_movable))
      return;
    std::vector<std::pair<double, double>> segments{{begin, end}};
    const auto blockers = horizontal ? candidates(lo, begin, hi, end)
                                     : candidates(begin, lo, end, hi);
    for (int k : blockers) {
      if (k == first || k == second)
        continue;
      const auto b = macros[k];
      const double low = horizontal ? b.xl : b.yl,
                   high = horizontal ? b.xh : b.yh;
      if (low >= hi || high <= lo)
        continue;
      const double a = horizontal ? b.yl : b.xl, z = horizontal ? b.yh : b.xh;
      std::vector<std::pair<double, double>> next;
      for (auto [s, t] : segments) {
        if (z <= s || a >= t)
          next.emplace_back(s, t);
        else {
          if (a > s)
            next.emplace_back(s, a);
          if (z < t)
            next.emplace_back(z, t);
        }
      }
      segments = std::move(next);
    }
    double force = 0;
    for (auto [s, t] : segments) {
      // A touching edge has no channel capacity. Use a narrow sampling strip;
      // overlap legality itself remains the responsibility of macro repulsion.
      const double pitch = horizontal ? dx : dy;
      const double a = gap > 0 ? lo : lo - pitch * 1e-6,
                   b = gap > 0 ? hi : hi + pitch * 1e-6;
      const double xl = horizontal ? a : s, xh = horizontal ? b : t;
      const double yl = horizontal ? s : a, yh = horizontal ? t : b;
      double pressure = 0, area = 0;
      bins(xl, yl, xh, yh, [&](int x, int y) {
        const double overlap = std::max(0., std::min(xh, l.xl + (x + 1) * dx) -
                                                std::max(xl, l.xl + x * dx)) *
                               std::max(0., std::min(yh, l.yl + (y + 1) * dy) -
                                                std::max(yl, l.yl + y * dy));
        const double utilization = std::clamp(u[x * by + y], 0., 3.);
        pressure += overlap * utilization * utilization;
        area += overlap;
      });
      if (area > 0)
        force += (t - s) / reference * pressure / area;
    }
    force *= std::pow(1 - gap / reach, 2);
    if (force == 0)
      return;
    ++r.channels;
    auto &gradient = horizontal ? r.x : r.y;
    if (first >= 0 && macros[first].node < d.num_movable)
      gradient[macros[first].node] += force;
    if (second >= 0 && macros[second].node < d.num_movable)
      gradient[macros[second].node] -= force;
  };
  for (int i = 0; i < int(macros.size()); ++i) {
    const auto a = macros[i];
    for (int j : candidates(a.xl - rx, a.yl - ry, a.xh + rx, a.yh + ry)) {
      if (j <= i)
        continue;
      const auto b = macros[j];
      const double ylo = std::max(a.yl, b.yl), yhi = std::min(a.yh, b.yh);
      const double xlo = std::max(a.xl, b.xl), xhi = std::min(a.xh, b.xh);
      if (a.xh <= b.xl)
        channel(i, j, true, a.xh, b.xl, ylo, yhi,
                std::min(a.yh - a.yl, b.yh - b.yl));
      if (b.xh <= a.xl)
        channel(j, i, true, b.xh, a.xl, ylo, yhi,
                std::min(a.yh - a.yl, b.yh - b.yl));
      if (a.yh <= b.yl)
        channel(i, j, false, a.yh, b.yl, xlo, xhi,
                std::min(a.xh - a.xl, b.xh - b.xl));
      if (b.yh <= a.yl)
        channel(j, i, false, b.yh, a.yl, xlo, xhi,
                std::min(a.xh - a.xl, b.xh - b.xl));
    }
    channel(-1, i, true, l.xl, a.xl, a.yl, a.yh, a.yh - a.yl);
    channel(i, -1, true, a.xh, l.xh, a.yl, a.yh, a.yh - a.yl);
    channel(-1, i, false, l.yl, a.yl, a.xl, a.xh, a.xh - a.xl);
    channel(i, -1, false, a.yh, l.yh, a.xl, a.xh, a.xh - a.xl);
  }
  return r;
}
} // namespace DreamPlace::native::detail
