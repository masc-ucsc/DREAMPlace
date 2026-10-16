// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/fillers.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace DreamPlace::native::macro {
double bounded_filler_width(const Design &d, const Options &o,
                            double preferred) {
  if (!o.add_fillers || !o.max_iterations || !d.num_movable)
    return preferred;
  const auto &l = d.layout;
  double movable = 0, fixed = 0;
  for (size_t i = 0; i < d.node_x.size(); ++i) {
    double hx =
        d.node_halo_x.empty()
            ? 0
            : std::ceil(d.node_halo_x[i] / l.site_width - 1e-9) * l.site_width;
    double hy =
        d.node_halo_y.empty()
            ? 0
            : std::ceil(d.node_halo_y[i] / l.row_height - 1e-9) * l.row_height;
    double x = d.node_x[i] - hx, y = d.node_y[i] - hy;
    double w = d.node_size_x[i] + 2 * hx;
    double h = d.node_size_y[i] + 2 * hy;
    if (i < static_cast<size_t>(d.num_movable))
      movable += w * h;
    else
      fixed += std::max(0., std::min(l.xh, x + w) - std::max(l.xl, x)) *
               std::max(0., std::min(l.yh, y + h) - std::max(l.yl, y));
  }
  const double width = l.xh - l.xl;
  const double target = std::max(
      0., o.target_density * (width * (l.yh - l.yl) - fixed) - movable);
  const double initial = std::clamp(
      preferred > 0 ? preferred : movable / d.num_movable / l.row_height,
      l.site_width, width);
  // Match place()'s one-million-node safety limit. Round upward to protect
  // against the final division rounding across the budget boundary.
  return std::clamp(
      std::max(initial,
               std::nextafter(target / (1000000. * l.row_height),
                              std::numeric_limits<double>::infinity())),
      l.site_width, width);
}
} // namespace DreamPlace::native::macro
