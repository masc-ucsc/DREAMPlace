// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <cmath>
#include <limits>

#include "dreamplace/native/placement.hpp"
#include "electric_potential/src/electric_density_map.cpp"

namespace DreamPlace::native::detail {
void density_map(const Design& d, int bx, int by, int threads, const std::vector<double>& sx,
                 const std::vector<double>& sy, const std::vector<double>& ox, const std::vector<double>& oy,
                 const std::vector<double>& ratio, std::vector<double>& map) {
  const auto& l = d.layout;
  map.assign(bx * by, 0);
  std::vector<double> cx(bx), cy(by);
  double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  for (int i = 0; i < bx; ++i) cx[i] = l.xl + (i + 0.5) * dx;
  for (int i = 0; i < by; ++i) cy[i] = l.yl + (i + 0.5) * dy;
  computeTriangleDensityMapLauncher<double, AtomicAdd<double>>(
      d.node_x.data(), d.node_y.data(), sx.data(), sy.data(), ox.data(), oy.data(), ratio.data(), cx.data(), cy.data(),
      d.node_x.size(), bx, by, l.xl, l.yl, l.xh, l.yh, dx, dy, threads, AtomicAdd<double>{}, map.data());
}
}  // namespace DreamPlace::native::detail
