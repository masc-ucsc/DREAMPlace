// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include "dreamplace/native/placement.hpp"

namespace DreamPlace::native::detail {
namespace {
// Deterministic spatial hash. Insert macro occupancy boxes only; fillers and
// standard cells never enter the pair enumeration. Each pair is returned once.
auto pairs(const Design &d) {
  std::vector<int> macros;
  std::vector<double> sides;
  for (int i = 0; i < static_cast<int>(d.node_macro.size()); ++i) {
    if (d.node_macro[i] && d.node_size_x[i] > 0 && d.node_size_y[i] > 0) {
      macros.push_back(i);
      sides.push_back(std::max(d.node_size_x[i], d.node_size_y[i]));
    }
  }
  std::set<std::pair<int, int>> result;
  if (macros.size() < 2) {
    return result;
  }
  std::sort(sides.begin(), sides.end());
  const double pitch = sides[sides.size() / 2];
  std::map<std::pair<int, int>, std::vector<int>> bins;
  for (int i : macros) {
    const int x0 = std::floor((d.node_x[i] - d.layout.xl) / pitch);
    const int y0 = std::floor((d.node_y[i] - d.layout.yl) / pitch);
    const int x1 =
        std::ceil((d.node_x[i] + d.node_size_x[i] - d.layout.xl) / pitch);
    const int y1 =
        std::ceil((d.node_y[i] + d.node_size_y[i] - d.layout.yl) / pitch);
    for (int x = x0; x < x1; ++x) {
      for (int y = y0; y < y1; ++y) {
        auto &bucket = bins[{x, y}];
        for (int j : bucket) {
          if (i < d.num_movable || j < d.num_movable) {
            result.emplace(j, i);
          }
        }
        bucket.push_back(i);
      }
    }
  }
  return result;
}
struct Intersection {
  double width, derivative;
};
Intersection intersection(double a, double wa, double b, double wb, int i,
                          int j) {
  const double right = std::min(a + wa, b + wb), left = std::max(a, b);
  if (right <= left) {
    return {0, 0};
  }
  const double center = (a + wa / 2) - (b + wb / 2);
  // At equal edges select a deterministic separating subgradient. Opposite
  // pair gradients preserve total force. Initialization jitter handles stacks
  // of identical connectivity without relying on floating point accidents.
  const double sign =
      center == 0 ? (i < j ? -1. : 1.) : std::copysign(1., center);
  const double dr =
      a + wa == b + wb ? (sign > 0 ? 0. : 1.) : double(a + wa < b + wb);
  const double dl = a == b ? (sign > 0 ? 1. : 0.) : double(a > b);
  return {right - left, dr - dl};
}
} // namespace
std::pair<int, int> automatic_bins(const Design &d) {
  double pitch = std::numeric_limits<double>::infinity(), area = 0;
  int cells = 0;
  for (int i = 0; i < static_cast<int>(d.node_x.size()); ++i) {
    if (d.node_size_x[i] <= 0 || d.node_size_y[i] <= 0) {
      continue;
    }
    if (!d.node_macro.empty() && d.node_macro[i]) {
      pitch = std::min({pitch, d.node_size_x[i], d.node_size_y[i]});
    } else {
      area += d.node_size_x[i] * d.node_size_y[i];
      ++cells;
    }
  }
  if (!std::isfinite(pitch)) {
    pitch =
        cells ? 4 * std::sqrt(area / cells)
              : std::max(d.layout.xh - d.layout.xl, d.layout.yh - d.layout.yl) /
                    32;
  }
  auto bins = [&](double length) {
    int n = 32;
    // Keep the native 32-bin baseline; automatic refinement is at most 2x.
    while (n < 64 && n < length / pitch) {
      n *= 2;
    }
    return n;
  };
  return {bins(d.layout.xh - d.layout.xl), bins(d.layout.yh - d.layout.yl)};
}
std::pair<int, int> dynamic_bins(const Design &d, double target, int maximum) {
  double area = 0;
  int count = 0;
  for (int i = 0; i < d.num_movable; ++i) {
    const double a = d.node_size_x[i] * d.node_size_y[i];
    if (a > 0) {
      area += a;
      ++count;
    }
  }
  if (!count)
    return {32, 32};
  const double width = d.layout.xh - d.layout.xl,
               height = d.layout.yh - d.layout.yl;
  const double ideal = width * height * target / (area / count);
  int ratio = 1;
  while (2 * ratio <= std::max(width, height) / std::min(width, height) &&
         2 * ratio <= maximum / 2)
    ratio *= 2;
  int n = 2;
  while (n < maximum / ratio && 4. * n * n * ratio <= ideal)
    n *= 2;
  return width > height ? std::pair{n * ratio, n} : std::pair{n, n * ratio};
}
Gradient macro_overlap(const Design &d, double weight) {
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  for (auto [i, j] : pairs(d)) {
    const auto x = intersection(d.node_x[i], d.node_size_x[i], d.node_x[j],
                                d.node_size_x[j], i, j);
    const auto y = intersection(d.node_y[i], d.node_size_y[i], d.node_y[j],
                                d.node_size_y[j], i, j);
    const double a = x.width * y.width;
    const double reference = std::min(d.node_size_x[i] * d.node_size_y[i],
                                      d.node_size_x[j] * d.node_size_y[j]);
    const double w = weight / (reference * std::sqrt(reference));
    g.value += w * a * a;
    const double gx = 2 * w * a * y.width * x.derivative,
                 gy = 2 * w * a * x.width * y.derivative;
    g.x[i] += gx;
    g.x[j] -= gx;
    g.y[i] += gy;
    g.y[j] -= gy;
  }
  return g;
}
// A bounded force toward the closest feasible escape from an exact fixed
// macro footprint. The Huber loss has distance units, matching HPWL. Its force
// is capped so a cell deep in a macro cannot destabilize the global step.
Gradient fixed_macro_overlap(const Design &d, int bx, int by, double weight,
                             double exclusion_scale) {
  Gradient g;
  const int n = d.node_x.size();
  g.x.resize(n);
  g.y.resize(n);
  if (weight == 0 || exclusion_scale == 0 || d.node_macro.empty())
    return g;
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  double max_w = 0, max_h = 0;
  for (int i = 0; i < d.num_movable; ++i) {
    if (!d.node_macro[i]) {
      max_w = std::max(max_w, d.node_size_x[i]);
      max_h = std::max(max_h, d.node_size_y[i]);
    }
  }
  const auto bin = detail::clamped_bin;
  std::vector<std::vector<int>> grid(bx * by);
  // Fillers appended past node_macro are never obstacles.
  const int obstacles = std::min<int>(n, d.node_macro.size());
  for (int j = d.num_movable; j < obstacles; ++j) {
    if (!d.node_macro[j] || d.node_size_x[j] <= 0 || d.node_size_y[j] <= 0)
      continue;
    // Query each cell by its lower-left point. Expansion by the maximum cell
    // extent makes every potentially intersecting obstacle reachable once.
    const double w = d.node_size_x[j] * exclusion_scale,
                 h = d.node_size_y[j] * exclusion_scale;
    const double mx = d.node_x[j] + (d.node_size_x[j] - w) / 2;
    const double my = d.node_y[j] + (d.node_size_y[j] - h) / 2;
    const int x0 = bin(mx - max_w, l.xl, dx, bx);
    const int x1 = bin(mx + w, l.xl, dx, bx);
    const int y0 = bin(my - max_h, l.yl, dy, by);
    const int y1 = bin(my + h, l.yl, dy, by);
    for (int x = x0; x <= x1; ++x)
      for (int y = y0; y <= y1; ++y)
        grid[x * by + y].push_back(j);
  }
  for (int i = 0; i < d.num_movable; ++i) {
    const double x = d.node_x[i], y = d.node_y[i];
    const double w = d.node_size_x[i], h = d.node_size_y[i];
    if (d.node_macro[i] || w <= 0 || h <= 0)
      continue;
    const auto &candidates =
        grid[bin(x, l.xl, dx, bx) * by + bin(y, l.yl, dy, by)];
    for (int j : candidates) {
      const double mw = d.node_size_x[j] * exclusion_scale,
                   mh = d.node_size_y[j] * exclusion_scale;
      const double mx = d.node_x[j] + (d.node_size_x[j] - mw) / 2;
      const double my = d.node_y[j] + (d.node_size_y[j] - mh) / 2;
      if (x + w <= mx || x >= mx + mw || y + h <= my || y >= my + mh)
        continue;
      const double distances[4] = {x + w - mx, mx + mw - x, y + h - my,
                                   my + mh - y};
      const bool feasible[4] = {mx - w >= l.xl, mx + mw <= l.xh - w,
                                my - h >= l.yl, my + mh <= l.yh - h};
      int direction = -1;
      // Stable index-derived tie-breaking separates coincident cells without
      // changing another cell's RNG sequence.
      for (int k = 0; k < 4; ++k) {
        const int side = (k + unsigned(i) * 3U + unsigned(j)) % 4;
        if (feasible[side] &&
            (direction < 0 || distances[side] < distances[direction]))
          direction = side;
      }
      if (direction < 0)
        continue; // An infeasible core is diagnosed by hard legalization.
      const double depth = distances[direction], extent = std::min(w, h);
      g.value += weight * (depth < extent ? depth * depth / (2 * extent)
                                          : depth - extent / 2);
      const double force = weight * std::min(depth / extent, 1.);
      if (direction == 0)
        g.x[i] += force;
      if (direction == 1)
        g.x[i] -= force;
      if (direction == 2)
        g.y[i] += force;
      if (direction == 3)
        g.y[i] -= force;
    }
  }
  return g;
}
std::pair<std::uint64_t, std::uint64_t> macro_overlap_counts(const Design &d) {
  std::uint64_t overlaps = 0, stacks = 0;
  for (auto [i, j] : pairs(d)) {
    if (std::min(d.node_x[i] + d.node_size_x[i],
                 d.node_x[j] + d.node_size_x[j]) <=
            std::max(d.node_x[i], d.node_x[j]) ||
        std::min(d.node_y[i] + d.node_size_y[i],
                 d.node_y[j] + d.node_size_y[j]) <=
            std::max(d.node_y[i], d.node_y[j])) {
      continue;
    }
    ++overlaps;
    if (std::abs((d.node_x[i] + d.node_size_x[i] / 2) -
                 (d.node_x[j] + d.node_size_x[j] / 2)) < d.layout.site_width &&
        std::abs((d.node_y[i] + d.node_size_y[i] / 2) -
                 (d.node_y[j] + d.node_size_y[j] / 2)) < d.layout.row_height) {
      ++stacks;
    }
  }
  return {overlaps, stacks};
}
} // namespace DreamPlace::native::detail
