// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace DreamPlace::native::detail {
AlignmentSnap macro_align_snap(const Design &d, double tau,
                               std::int64_t weight) {
  if (!std::isfinite(tau) || tau < 0 || weight <= 0 || weight > 1000000)
    throw std::invalid_argument("invalid macro snap tolerance/weight");
  AlignmentSnap snap;
  auto pairs = macro_align_pairs(d, false);
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
      if (p.axis != axis || size[p.a] != size[p.b])
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
      snap.groups.push_back(
          {d.node_align_class[members.front()], axis, target, members});
      for (size_t k = 1; k < members.size(); ++k)
        snap.constraints.push_back(
            {members[k - 1], members[k], x ? 0 : 1, weight});
    }
  }
  return snap;
}

AlignmentReport macro_alignment_report(const Design &before,
                                       const Design &after,
                                       const AlignmentSnap &snap) {
  AlignmentReport r;
  r.groups = snap.groups;
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
    r.soft_violations += std::abs(pos[c.a] - pos[c.b]) > 1e-6;
  }
  for (int i = 0; i < before.num_movable; ++i)
    if (!before.node_macro.empty() && before.node_macro[i])
      r.macro_displacements.push_back(
          std::hypot(after.node_x[i] - before.node_x[i],
                     after.node_y[i] - before.node_y[i]));
  return r;
}
} // namespace DreamPlace::native::detail
