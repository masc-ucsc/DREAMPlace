// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <cmath>

namespace DreamPlace::native::detail {
Gradient cohesion(const Design &d, double lambda) {
  Gradient out;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  if (!lambda)
    return out;
  for (const auto &group : d.cohesion_groups) {
    if (group.nodes.size() < 2 || !group.weight)
      continue;
    double area = 0, cx = 0, cy = 0;
    for (size_t k = 0; k < group.nodes.size(); ++k) {
      const int i = group.nodes[k];
      const double a = group.areas.empty() ? d.node_size_x[i] * d.node_size_y[i]
                                           : group.areas[k];
      area += a;
      double x = d.node_x[i] + d.node_size_x[i] / 2 - d.layout.xl;
      double y = d.node_y[i] + d.node_size_y[i] / 2 - d.layout.yl;
      if (!group.offsets.empty()) {
        x += group.offsets[k].first;
        y += group.offsets[k].second;
      }
      cx += a * x;
      cy += a * y;
    }
    cx /= area;
    cy /= area;
    for (size_t k = 0; k < group.nodes.size(); ++k) {
      const int i = group.nodes[k];
      const double a = group.areas.empty() ? d.node_size_x[i] * d.node_size_y[i]
                                           : group.areas[k];
      double dx = d.node_x[i] + d.node_size_x[i] / 2 - d.layout.xl - cx;
      double dy = d.node_y[i] + d.node_size_y[i] / 2 - d.layout.yl - cy;
      if (!group.offsets.empty()) {
        dx += group.offsets[k].first;
        dy += group.offsets[k].second;
      }
      const double weight = lambda * group.weight * (a / area);
      out.value += weight * (dx * dx + dy * dy);
      out.x[i] += 2 * weight * dx;
      out.y[i] += 2 * weight * dy;
    }
  }
  return out;
}
} // namespace DreamPlace::native::detail
