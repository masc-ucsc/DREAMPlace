// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>

namespace DreamPlace::native::detail {
AlignmentSnap macro_align_snap(const Design &d, double tau,
                               std::int64_t weight) {
  if (!std::isfinite(tau) || tau < 0 || weight <= 0 || weight > 1000000)
    throw std::invalid_argument("invalid macro snap tolerance/weight");
  AlignmentSnap snap;
  const auto &classes =
      d.node_snap_class.empty() ? d.node_align_class : d.node_snap_class;
  auto pairs = macro_align_pairs(d, false, &classes);
  for (Axis axis : {Axis::x, Axis::y}) {
    bool x = axis == Axis::x;
    const auto &pos = x ? d.node_x : d.node_y;
    const auto &size = x ? d.node_size_x : d.node_size_y;
    std::vector<int> parent(d.num_movable);
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](int a) {
      while (parent[a] != a) {
        parent[a] = parent[parent[a]];
        a = parent[a];
      }
      return a;
    };
    std::vector<bool> active(d.num_movable);
    for (const auto &p : pairs) {
      if (p.axis != axis || d.node_size_x[p.a] != d.node_size_x[p.b] ||
          d.node_size_y[p.a] != d.node_size_y[p.b])
        continue;
      double delta = pos[p.a] + size[p.a] / 2 - pos[p.b] - size[p.b] / 2;
      if (std::abs(delta) <= tau * size[p.a]) {
        int a = root(p.a), b = root(p.b);
        parent[std::max(a, b)] = std::min(a, b);
        active[p.a] = active[p.b] = true;
      }
    }
    std::map<int, std::vector<int>> components;
    for (int i = 0; i < d.num_movable; ++i)
      if (active[i])
        components[root(i)].push_back(i);
    double pitch = x ? d.layout.site_width : d.layout.row_height;
    double lo = x ? d.layout.xl : d.layout.yl;
    for (auto &[id, members] : components) {
      (void)id;
      if (members.size() < 2)
        continue;
      std::vector<double> anchors;
      for (int i : members)
        anchors.push_back(pos[i]);
      std::sort(anchors.begin(), anchors.end());
      // Same-master members have equal extents. Round the median origin to
      // the placement lattice; rounding the center can produce half-site
      // origins.
      double target =
          lo +
          std::round((anchors[(anchors.size() - 1) / 2] - lo) / pitch) * pitch;
      snap.groups.push_back({classes[members.front()], axis, target, members});
      for (size_t k = 1; k < members.size(); ++k)
        snap.constraints.push_back(
            {members[k - 1], members[k], x ? 0 : 1, weight});
    }
  }
  return snap;
}

// The same fitter used for FP centroids supplies line, spacing and midpoint
// targets. Equal-area circles describe halo-padded macros; they are NOT hard
// physical fences. The ordinary rectangular macro LP retains final legality.
AlignmentSnap macro_regular_snap(const Design &d, const FloorplanOptions &o,
                                 std::int64_t weight) {
  if (weight <= 0 || weight > 1000000)
    throw std::invalid_argument("invalid macro regularity weight");
  d.validate();
  AlignmentSnap snap;
  if (d.node_snap_class.empty() || d.node_macro.empty())
    return snap;
  std::map<int, std::vector<int>> classes;
  for (int n = 0; n < d.num_movable; ++n)
    if (d.node_macro[n] && d.node_snap_class[n] >= 0)
      classes[d.node_snap_class[n]].push_back(n);
  std::vector<FloorplanCircle> circles;
  std::vector<int> nodes;
  for (const auto &[id, members] : classes) {
    if (members.size() < 2)
      continue;
    const int first = members.front();
    for (int n : members) {
      const double w = d.node_size_x[n], h = d.node_size_y[n];
      if (w != d.node_size_x[first] || h != d.node_size_y[first])
        throw std::invalid_argument("macro snap class must have one size");
      circles.push_back({std::to_string(n), d.node_x[n] + w / 2,
                         d.node_y[n] + h / 2, std::sqrt(w * h / std::acos(-1.)),
                         w * h, id, -1});
      nodes.push_back(n);
    }
  }
  snap.regularity = snap_floorplan(circles, d.layout.xl, d.layout.yl,
                                   d.layout.xh, d.layout.yh, o);
  std::vector<double> tx(nodes.size()), ty(nodes.size());
  auto lattice = [](double center, double size, double lo, double hi,
                    double pitch) {
    const double last = std::floor((hi - lo - size) / pitch);
    return lo + std::clamp(std::round((center - size / 2 - lo) / pitch), 0.,
                           std::max(0., last)) *
                    pitch;
  };
  for (size_t k = 0; k < nodes.size(); ++k) {
    int n = nodes[k];
    const auto &c = snap.regularity.circles[k];
    tx[k] = lattice(c.x, d.node_size_x[n], d.layout.xl, d.layout.xh,
                    d.layout.site_width);
    ty[k] = lattice(c.y, d.node_size_y[n], d.layout.yl, d.layout.yh,
                    d.layout.row_height);
  }
  std::set<std::tuple<int, int, int>> emitted;
  std::set<int> active;
  for (const auto &relation : snap.regularity.relations) {
    auto members = relation.nodes;
    std::sort(members.begin(), members.end());
    members.erase(std::unique(members.begin(), members.end()), members.end());
    const auto &target = relation.axis == 0 ? tx : ty;
    std::stable_sort(members.begin(), members.end(), [&](int a, int b) {
      return std::pair(target[a], a) < std::pair(target[b], b);
    });
    std::vector<int> group;
    for (int k : members) {
      active.insert(k);
      group.push_back(nodes[k]);
    }
    if (relation.kind == "alignment")
      snap.groups.push_back({circles[members.front()].family,
                             relation.axis == 0 ? Axis::x : Axis::y,
                             target[members.front()], group});
    // Keep the fitted relative positions soft and translation-invariant. Use
    // chains rather than all-pairs, including nonzero spacing/midpoint offsets.
    for (size_t k = 1; k < members.size(); ++k) {
      int a = members[k - 1], b = members[k];
      if (emitted.emplace(relation.axis, nodes[a], nodes[b]).second)
        snap.constraints.push_back(
            {nodes[a], nodes[b], relation.axis, weight, target[b] - target[a]});
    }
  }
  for (int k : active)
    snap.targets.push_back({nodes[k], tx[k], ty[k]});
  return snap;
}

AlignmentReport macro_alignment_report(const Design &before,
                                       const Design &after,
                                       const AlignmentSnap &snap) {
  AlignmentReport r;
  r.groups = snap.groups;
  r.regularity = snap.regularity;
  r.constraints = snap.constraints;
  std::vector<unsigned char> axes(before.num_movable);
  for (const auto &g : snap.groups) {
    bool x = g.axis == Axis::x;
    ++(x ? r.x_groups : r.y_groups);
    const auto &pos = x ? after.node_x : after.node_y;
    bool exact = true;
    for (int i : g.nodes) {
      axes[i] |= x ? 1 : 2;
      if (std::abs(pos[i] - pos[g.nodes.front()]) > 1e-6)
        exact = false;
    }
    ++(exact ? r.exact_groups : r.broken_groups);
  }
  for (auto v : axes) {
    r.members += v != 0;
    r.both_axes += v == 3;
  }
  for (const auto &c : snap.constraints) {
    const auto &pos = c.axis == 0 ? after.node_x : after.node_y;
    r.soft_violations += std::abs(pos[c.b] - pos[c.a] - c.offset) > 1e-6;
  }
  for (int i = 0; i < before.num_movable; ++i)
    if (!before.node_macro.empty() && before.node_macro[i])
      r.macro_displacements.push_back(
          std::hypot(after.node_x[i] - before.node_x[i],
                     after.node_y[i] - before.node_y[i]));
  return r;
}
} // namespace DreamPlace::native::detail
