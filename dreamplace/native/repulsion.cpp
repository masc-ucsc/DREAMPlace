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
