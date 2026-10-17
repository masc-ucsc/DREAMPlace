// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>

namespace DreamPlace::native::detail {
std::vector<AlignmentPair> macro_align_pairs(const Design &d, bool nearest) {
  std::map<int, std::vector<int>> classes;
  for (int i = 0; i < d.num_movable; ++i)
    if (!d.node_align_class.empty() && d.node_align_class[i] >= 0 &&
        !d.node_macro.empty() && d.node_macro[i])
      classes[d.node_align_class[i]].push_back(i);
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
                         double *residual) {
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
    double w = std::min(size[p.a], size[p.b]);
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
