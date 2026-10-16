// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/check.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_set>

namespace DreamPlace::native::macro {
CheckResult check(const Design &d, std::span<const int> nodes,
                  double tolerance) {
  d.validate();
  if (!std::isfinite(tolerance) || tolerance < 0)
    throw std::invalid_argument("invalid macro-check tolerance");
  struct Box {
    double xl, yl, xh, yh;
    int id;
  };
  std::vector<Box> boxes;
  std::unordered_set<int> seen;
  CheckResult result;
  const auto &l = d.layout;
  auto aligned = [&](double coordinate, double origin, double pitch) {
    return std::abs(coordinate - origin -
                    std::round((coordinate - origin) / pitch) * pitch) <=
           tolerance;
  };
  for (int id : nodes) {
    if (id < 0 || id >= static_cast<int>(d.node_x.size()) ||
        !seen.insert(id).second)
      throw std::invalid_argument("invalid or duplicate macro index");
    if (d.node_size_x[id] <= 0 || d.node_size_y[id] <= 0)
      throw std::invalid_argument("macro must have positive area");
    const double hx =
        d.node_halo_x.empty()
            ? 0
            : std::ceil(d.node_halo_x[id] / l.site_width - 1e-9) * l.site_width;
    const double hy =
        d.node_halo_y.empty()
            ? 0
            : std::ceil(d.node_halo_y[id] / l.row_height - 1e-9) * l.row_height;
    Box b{d.node_x[id] - hx, d.node_y[id] - hy,
          d.node_x[id] + d.node_size_x[id] + hx,
          d.node_y[id] + d.node_size_y[id] + hy, id};
    if (!std::isfinite(b.xl) || !std::isfinite(b.yl) || !std::isfinite(b.xh) ||
        !std::isfinite(b.yh))
      throw std::invalid_argument("macro envelope overflow");
    if (b.xl < l.xl - tolerance || b.yl < l.yl - tolerance ||
        b.xh > l.xh + tolerance || b.yh > l.yh + tolerance)
      result.violations.push_back({ViolationKind::bounds, id});
    if (!aligned(d.node_x[id], l.xl, l.site_width) ||
        !aligned(d.node_y[id], l.yl, l.row_height))
      result.violations.push_back({ViolationKind::pitch, id});
    boxes.push_back(b);
  }
  std::sort(boxes.begin(), boxes.end(), [](const Box &a, const Box &b) {
    return a.xl < b.xl || (a.xl == b.xl && a.id < b.id);
  });
  std::vector<Box> active;
  for (const auto &a : boxes) {
    std::erase_if(active,
                  [&](const Box &b) { return b.xh <= a.xl + tolerance; });
    for (const auto &b : active)
      if (a.yl < b.yh - tolerance && b.yl < a.yh - tolerance)
        result.violations.push_back({ViolationKind::overlap,
                                     std::min(a.id, b.id),
                                     std::max(a.id, b.id)});
    active.push_back(a);
  }
  return result;
}
} // namespace DreamPlace::native::macro
