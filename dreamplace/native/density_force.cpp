// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <cmath>
#include <limits>

#include "dreamplace/native/placement.hpp"
#include "electric_potential/src/electric_force.cpp"

namespace DreamPlace::native::detail {
Gradient density_force(const Design& d, int bx, int by, int threads, const std::vector<double>& sx,
                       const std::vector<double>& sy, const std::vector<double>& ox, const std::vector<double>& oy,
                       const std::vector<double>& ratio, const std::vector<double>& fx, const std::vector<double>& fy) {
  const auto& l = d.layout;
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  computeElectricForceLauncher<double>(bx, by, bx, by, fx.data(), fy.data(), d.node_x.data(), d.node_y.data(),
                                       sx.data(), sy.data(), ox.data(), oy.data(), ratio.data(), nullptr, nullptr, l.xl,
                                       l.yl, l.xh, l.yh, (l.xh - l.xl) / bx, (l.yh - l.yl) / by, d.node_x.size(),
                                       threads, g.x.data(), g.y.data());
  // The kernel returns electrostatic force; the minimizer needs its negative.
  for (size_t i = 0; i < g.x.size(); ++i) {
    g.x[i] = -g.x[i];
    g.y[i] = -g.y[i];
  }
  return g;
}
}  // namespace DreamPlace::native::detail
