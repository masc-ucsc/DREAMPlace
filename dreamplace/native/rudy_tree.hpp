// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/rudy.hpp"
namespace DreamPlace::native::detail {
struct RouteEdge {
  double xl, yl, xh, yh;
};
// Empty above the cap; caller preserves the weighted bounding-box fallback.
std::vector<RouteEdge>
route_tree(const std::vector<std::pair<double, double>> &, const RudyOptions &);
} // namespace DreamPlace::native::detail
