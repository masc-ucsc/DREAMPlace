// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <tuple>

namespace DreamPlace::native::detail {
namespace {
std::vector<AlignmentPair>
select_pairs(const Design &d, std::vector<AlignmentPair> pairs, bool nearest) {
  if (!nearest)
    return pairs;
  using Key = std::tuple<int, int, bool>;
  std::map<Key, std::pair<double, size_t>> best;
  for (size_t k = 0; k < pairs.size(); ++k) {
    const auto &p = pairs[k];
    // Neighbors above/below a column, or left/right along a row.
    const bool column = p.axis == Axis::x;
    const auto &pos = column ? d.node_y : d.node_x;
    const auto &size = column ? d.node_size_y : d.node_size_x;
    double delta = pos[p.b] + size[p.b] / 2 - pos[p.a] - size[p.a] / 2;
    for (auto [a, sign] :
         {std::pair{p.a, delta > 0}, std::pair{p.b, delta < 0}}) {
      Key key{a, static_cast<int>(p.axis), sign};
      auto candidate = std::pair{std::abs(delta), k};
      if (!best.contains(key) || candidate < best[key])
        best[key] = candidate;
    }
  }
  std::set<size_t> keep;
  for (const auto &[key, v] : best) {
    (void)key;
    keep.insert(v.second);
  }
  std::vector<AlignmentPair> selected;
  for (auto k : keep)
    selected.push_back(pairs[k]);
  return selected;
}
} // namespace
std::vector<AlignmentPair>
macro_align_pairs(const Design &d, bool nearest,
                  const std::vector<int> *class_override) {
  const auto &node_class =
      class_override ? *class_override : d.node_align_class;
  std::map<int, std::vector<int>> classes;
  for (int i = 0; i < d.num_movable; ++i)
    if (!node_class.empty() && node_class[i] >= 0 && !d.node_macro.empty() &&
        d.node_macro[i])
      classes[node_class[i]].push_back(i);
  std::vector<AlignmentPair> pairs;
  for (const auto &[id, nodes] : classes) {
    (void)id;
    for (size_t j = 0; j < nodes.size(); ++j)
      for (size_t k = j + 1; k < nodes.size(); ++k) {
        int a = nodes[j], b = nodes[k];
        bool x = std::min(d.node_x[a] + d.node_size_x[a],
                          d.node_x[b] + d.node_size_x[b]) >
                 std::max(d.node_x[a], d.node_x[b]);
        bool y = std::min(d.node_y[a] + d.node_size_y[a],
                          d.node_y[b] + d.node_size_y[b]) >
                 std::max(d.node_y[a], d.node_y[b]);
        if (x != y)
          pairs.push_back({a, b, x ? Axis::x : Axis::y});
      }
  }
  return select_pairs(d, std::move(pairs), nearest);
}

// Strip bins are only a broad-phase index, not snap coordinates. Search all
// intersecting center bins, including neighboring bins at strip boundaries;
// exact physical projection tests determine whether a force is applicable.
std::vector<AlignmentPair>
macro_grid_align_pairs(const Design &d, bool nearest,
                       const std::vector<int> *class_override,
                       bool padded_extents, bool overlapping_pairs) {
  struct Box {
    int node;
    double x, y, w, h;
  };
  std::map<int, std::vector<Box>> classes;
  const auto &node_class =
      class_override ? *class_override : d.node_align_class;
  double pitch_x = std::numeric_limits<double>::infinity(), pitch_y = pitch_x;
  for (int i = 0; i < d.num_movable; ++i) {
    if (node_class.empty() || node_class[i] < 0 || d.node_macro.empty() ||
        !d.node_macro[i])
      continue;
    const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
    const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
    const auto halo = [&](const std::vector<double> &values, double pitch) {
      return values.empty() ? 0 : std::ceil(values[i] / pitch - 1e-9) * pitch;
    };
    const double hx = padded_extents && d.node_pad_x.empty()
                          ? halo(d.node_halo_x, d.layout.site_width)
                          : 0;
    const double hy = padded_extents && d.node_pad_y.empty()
                          ? halo(d.node_halo_y, d.layout.row_height)
                          : 0;
    const double w = d.node_size_x[i] - (padded_extents ? 0 : 2 * px) + 2 * hx;
    const double h = d.node_size_y[i] - (padded_extents ? 0 : 2 * py) + 2 * hy;
    if (w <= 0 || h <= 0)
      continue;
    classes[node_class[i]].push_back({i, d.node_x[i] + d.node_size_x[i] / 2,
                                      d.node_y[i] + d.node_size_y[i] / 2, w,
                                      h});
    pitch_x = std::min(pitch_x, w);
    pitch_y = std::min(pitch_y, h);
  }
  std::vector<AlignmentPair> pairs;
  for (const auto &[id, nodes] : classes) {
    (void)id;
    for (bool x : {true, false}) {
      const double pitch = x ? pitch_x : pitch_y;
      const double origin = x ? d.layout.xl : d.layout.yl;
      // Floating integer keys avoid narrowing overflow for extreme legal
      // ratios of core span to macro width. No dense floorplan-sized array.
      const auto bin = [&](double p) {
        return std::floor((p - origin) / pitch);
      };
      std::map<double, std::vector<int>> bins;
      double largest = 0;
      for (int k = 0; k < int(nodes.size()); ++k) {
        const auto &a = nodes[k];
        bins[bin(x ? a.x : a.y)].push_back(k);
        largest = std::max(largest, x ? a.w : a.h);
      }
      for (const auto &a : nodes) {
        const double center = x ? a.x : a.y;
        const double reach = ((x ? a.w : a.h) + largest) / 2;
        const double last = bin(center + reach);
        for (auto it = bins.lower_bound(bin(center - reach));
             it != bins.end() && it->first <= last; ++it) {
          for (int k : it->second) {
            const auto &b = nodes[k];
            if (b.node <= a.node)
              continue;
            const bool ox = std::abs(a.x - b.x) < (a.w + b.w) / 2;
            const bool oy = std::abs(a.y - b.y) < (a.h + b.h) / 2;
            if ((overlapping_pairs || ox != oy) && (x ? ox : oy))
              pairs.push_back({a.node, b.node, x ? Axis::x : Axis::y});
          }
        }
      }
    }
  }
  std::sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b) {
    return std::tie(a.a, a.b, a.axis) < std::tie(b.a, b.b, b.axis);
  });
  return select_pairs(d, std::move(pairs), nearest);
}

double macro_align_ramp(const Options &o, int iteration) {
  if (iteration <= o.macro_align_start)
    return 0;
  if (iteration >= o.macro_align_full)
    return 1;
  double s = double(iteration - o.macro_align_start) /
             (o.macro_align_full - o.macro_align_start);
  return s * s;
}

Gradient macro_alignment(const Design &d,
                         const std::vector<AlignmentPair> &pairs,
                         double lambda_x, double lambda_y, double sigma,
                         double *residual, bool physical_extents) {
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  std::vector<int> degree(d.node_x.size());
  for (const auto &p : pairs) {
    ++degree[p.a];
    ++degree[p.b];
  }
  double total = 0;
  for (const auto &p : pairs) {
    bool x = p.axis == Axis::x;
    const auto &pos = x ? d.node_x : d.node_y;
    const auto &size = x ? d.node_size_x : d.node_size_y;
    const auto &pad = x ? d.node_pad_x : d.node_pad_y;
    const auto extent = [&](int i) {
      return size[i] - (physical_extents && !pad.empty() ? 2 * pad[i] : 0);
    };
    double w = std::min(extent(p.a), extent(p.b));
    double delta = pos[p.a] + size[p.a] / 2 - pos[p.b] - size[p.b] / 2;
    double lambda = (x ? lambda_x : lambda_y) /
                    (w * std::sqrt(double(degree[p.a]) * degree[p.b]));
    double grad;
    if (sigma > 0) {
      double sw = sigma * w, z = delta / sw, decay = std::exp(-z * z);
      g.value += lambda * sw * sw * -std::expm1(-z * z);
      grad = 2 * lambda * delta * decay;
    } else {
      g.value += lambda * delta * delta;
      grad = 2 * lambda * delta;
    }
    auto &axis = x ? g.x : g.y;
    axis[p.a] += grad;
    axis[p.b] -= grad;
    total += std::abs(delta) / w;
  }
  if (residual)
    *residual = pairs.empty() ? 0 : total / pairs.size();
  return g;
}
} // namespace DreamPlace::native::detail
