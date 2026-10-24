// BSD-3-Clause. See LICENSE.
#include "local_density.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace DreamPlace::native::detail {
double macro_cell_screening_length(const Design &d, double reach) {
  double side = 0;
  for (size_t i = 0; i < d.node_x.size(); ++i) {
    if (i < d.node_macro.size() && d.node_macro[i]) {
      const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
      const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
      side = std::max({side, d.node_size_x[i] - 2 * px,
                       d.node_size_y[i] - 2 * py});
    }
  }
  return reach * side;
}
LocalDensity::LocalDensity(const Design &d, int bins, double range,
                           double target, int num_threads)
    : bins_(bins), physical_nodes_(d.node_x.size()),
      dx_((d.layout.xh - d.layout.xl) / bins),
      dy_((d.layout.yh - d.layout.yl) / bins), target_(target),
      solver_(bins, bins, d.layout.xh - d.layout.xl, d.layout.yh - d.layout.yl,
              true, range * std::sqrt(dx_ * dy_), num_threads) {}
std::vector<double> LocalDensity::cell_density(const Design &d) const {
  std::vector<double> rho(bins_ * bins_);
  // Cloud-in-cell scatter and its transpose gather suppress bin-boundary jumps.
  const auto visit = [&](int i, auto fn) {
    const double px =
        (d.node_x[i] + d.node_size_x[i] / 2 - d.layout.xl) / dx_ - .5;
    const double py =
        (d.node_y[i] + d.node_size_y[i] / 2 - d.layout.yl) / dy_ - .5;
    const int ix = std::floor(px), iy = std::floor(py);
    const double tx = px - ix, ty = py - iy;
    for (int x = 0; x < 2; ++x) {
      for (int y = 0; y < 2; ++y) {
        fn(std::clamp(ix + x, 0, bins_ - 1) * bins_ +
               std::clamp(iy + y, 0, bins_ - 1),
           (x ? tx : 1 - tx) * (y ? ty : 1 - ty));
      }
    }
  };
  const auto cell = [&](int i) {
    return (d.node_macro.empty() || i >= int(d.node_macro.size()) ||
            !d.node_macro[i]) &&
           d.node_size_x[i] > 0 && d.node_size_y[i] > 0;
  };
  // Fillers are excluded: host-controlled local spreading never creates them.
  const int original = std::min(physical_nodes_, int(d.node_x.size()));
  for (int i = 0; i < original; ++i) {
    if (cell(i)) {
      const double area = d.node_size_x[i] * d.node_size_y[i];
      visit(i, [&](int b, double w) { rho[b] += w * area / (dx_ * dy_); });
    }
  }
  return rho;
}
Gradient LocalDensity::gradient(const Design &d, bool macros_only,
                                const std::vector<double> &source,
                                bool ignore_fixed_macros,
                                double fixed_macro_scale) const {
  auto rho = source.empty() ? cell_density(d) : source;
  if (rho.size() != size_t(bins_ * bins_)) {
    throw std::invalid_argument("invalid local density source");
  }
  std::vector<double> phi, fx, fy;
  // Cloud-in-cell scatter and its transpose gather suppress bin-boundary jumps.
  const auto visit = [&](int i, auto fn) {
    const double px =
        (d.node_x[i] + d.node_size_x[i] / 2 - d.layout.xl) / dx_ - .5;
    const double py =
        (d.node_y[i] + d.node_size_y[i] / 2 - d.layout.yl) / dy_ - .5;
    const int ix = std::floor(px), iy = std::floor(py);
    const double tx = px - ix, ty = py - iy;
    for (int x = 0; x < 2; ++x) {
      for (int y = 0; y < 2; ++y) {
        fn(std::clamp(ix + x, 0, bins_ - 1) * bins_ +
               std::clamp(iy + y, 0, bins_ - 1),
           (x ? tx : 1 - tx) * (y ? ty : 1 - ty));
      }
    }
  };
  const auto cell = [&](int i) {
    return (d.node_macro.empty() || i >= int(d.node_macro.size()) ||
            !d.node_macro[i]) &&
           d.node_size_x[i] > 0 && d.node_size_y[i] > 0;
  };
  // Fillers are excluded: host-controlled local spreading never creates them.
  const int original = std::min(physical_nodes_, int(d.node_x.size()));
  // Emit a distance-decaying field from local overpopulation, rather than
  // amplifying arbitrary sparse-region noise. The field extends outside an
  // overfull bin; no physical overlap is needed for a cell to receive force.
  std::vector<double> blocked(rho.size());
  if (!macros_only && target_ > 0 && !d.node_macro.empty()) {
    for (int i = 0; i < std::min(original, int(d.node_macro.size())); ++i) {
      if (d.node_macro[i] && !(ignore_fixed_macros && i >= d.num_movable)) {
        const double scale = i >= d.num_movable ? fixed_macro_scale : 1;
        const double w = d.node_size_x[i] * scale, h = d.node_size_y[i] * scale;
        const double x = d.node_x[i] + (d.node_size_x[i] - w) / 2;
        const double y = d.node_y[i] + (d.node_size_y[i] - h) / 2;
        const int x0 =
            std::clamp(int(std::floor((x - d.layout.xl) / dx_)), 0, bins_);
        const int x1 =
            std::clamp(int(std::ceil((x + w - d.layout.xl) / dx_)), 0, bins_);
        const int y0 =
            std::clamp(int(std::floor((y - d.layout.yl) / dy_)), 0, bins_);
        const int y1 =
            std::clamp(int(std::ceil((y + h - d.layout.yl) / dy_)), 0, bins_);
        for (int a = x0; a < x1; ++a) {
          for (int b = y0; b < y1; ++b) {
            const double ox =
                std::max(0., std::min(x + w, d.layout.xl + (a + 1) * dx_) -
                                 std::max(x, d.layout.xl + a * dx_));
            const double oy =
                std::max(0., std::min(y + h, d.layout.yl + (b + 1) * dy_) -
                                 std::max(y, d.layout.yl + b * dy_));
            blocked[a * bins_ + b] += ox * oy / (dx_ * dy_);
          }
        }
      }
    }
  }
  if (!macros_only) {
    for (size_t b = 0; b < rho.size(); ++b) {
      rho[b] = std::max(0., rho[b] - target_ * std::max(0., 1 - blocked[b]));
    }
  }
  solver_.solve(rho, phi, fx, fy);
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  if (macros_only) {
    // A center-only sample misses cell pressure near a large SRAM's edges.
    // Integrate the field over the entire padded footprint instead. Raw cell
    // area is used here: a macro should prefer a quieter neighborhood even
    // when both sides are below the standard-cell soft density target.
    for (int i = 0; i < std::min(original, d.num_movable); ++i) {
      if (d.node_macro.empty() || !d.node_macro[i]) {
        continue;
      }
      const double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
                   h = d.node_size_y[i];
      const int x0 =
          std::clamp(int(std::floor((x - d.layout.xl) / dx_)), 0, bins_);
      const int x1 =
          std::clamp(int(std::ceil((x + w - d.layout.xl) / dx_)), 0, bins_);
      const int y0 =
          std::clamp(int(std::floor((y - d.layout.yl) / dy_)), 0, bins_);
      const int y1 =
          std::clamp(int(std::ceil((y + h - d.layout.yl) / dy_)), 0, bins_);
      for (int a = x0; a < x1; ++a) {
        for (int b = y0; b < y1; ++b) {
          const double ox =
              std::max(0., std::min(x + w, d.layout.xl + (a + 1) * dx_) -
                               std::max(x, d.layout.xl + a * dx_));
          const double oy =
              std::max(0., std::min(y + h, d.layout.yl + (b + 1) * dy_) -
                               std::max(y, d.layout.yl + b * dy_));
          g.x[i] -= ox * oy * fx[a * bins_ + b];
          g.y[i] -= ox * oy * fy[a * bins_ + b];
        }
      }
    }
    return g;
  }
  for (int i = 0; i < std::min(original, d.num_movable); ++i) {
    if (cell(i)) {
      const double area = d.node_size_x[i] * d.node_size_y[i];
      visit(i, [&](int b, double w) {
        g.x[i] -= w * area * fx[b];
        g.y[i] -= w * area * fy[b];
      });
    }
  }
  return g;
}
Gradient LocalDensity::density_pull(const Design &d,
                                    const std::vector<double> &source,
                                    double range) const {
  if (source.size() != size_t(bins_ * bins_) || !std::isfinite(range) ||
      range <= 0) {
    throw std::invalid_argument("invalid density pull source/range");
  }
  const int stride = bins_ + 1;
  std::vector<double> prefix(stride * stride);
  for (int x = 0; x < bins_; ++x) {
    for (int y = 0; y < bins_; ++y) {
      prefix[(x + 1) * stride + y + 1] =
          source[x * bins_ + y] + prefix[x * stride + y + 1] +
          prefix[(x + 1) * stride + y] - prefix[x * stride + y];
    }
  }
  // Bilinear interpolation of prefix values exactly integrates a piecewise
  // constant density grid, including partially covered boundary bins.
  const auto integral = [&](double x, double y) {
    x = std::clamp((x - d.layout.xl) / dx_, 0., double(bins_));
    y = std::clamp((y - d.layout.yl) / dy_, 0., double(bins_));
    const int ix = std::min(int(x), bins_ - 1),
              iy = std::min(int(y), bins_ - 1);
    const double tx = x - ix, ty = y - iy;
    return (1 - tx) * ((1 - ty) * prefix[ix * stride + iy] +
                       ty * prefix[ix * stride + iy + 1]) +
           tx * ((1 - ty) * prefix[(ix + 1) * stride + iy] +
                 ty * prefix[(ix + 1) * stride + iy + 1]);
  };
  const auto mean = [&](double xl, double yl, double xh, double yh) {
    xl = std::max(xl, d.layout.xl);
    yl = std::max(yl, d.layout.yl);
    xh = std::min(xh, d.layout.xh);
    yh = std::min(yh, d.layout.yh);
    if (xh <= xl || yh <= yl) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return (integral(xh, yh) - integral(xl, yh) - integral(xh, yl) +
            integral(xl, yl)) /
           ((xh - xl) * (yh - yl) / (dx_ * dy_));
  };
  Gradient g;
  g.x.resize(d.node_x.size());
  g.y.resize(d.node_y.size());
  const double reach = range * std::sqrt(dx_ * dy_);
  for (int i = 0; i < d.num_movable; ++i) {
    if (d.node_macro.empty() || !d.node_macro[i]) {
      continue;
    }
    const double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
                 h = d.node_size_y[i];
    if (w <= 0 || h <= 0) {
      continue;
    }
    const double left = mean(x - reach, y - reach, x, y + h + reach);
    const double right = mean(x + w, y - reach, x + w + reach, y + h + reach);
    const double bottom = mean(x - reach, y - reach, x + w + reach, y);
    const double top = mean(x - reach, y + h, x + w + reach, y + h + reach);
    // Missing off-core strips cannot masquerade as empty, attractive space.
    if (std::isfinite(left) && std::isfinite(right)) {
      g.x[i] = w * h * (right - left) / (w + 2 * reach);
    }
    if (std::isfinite(bottom) && std::isfinite(top)) {
      g.y[i] = w * h * (top - bottom) / (h + 2 * reach);
    }
  }
  return g;
}
} // namespace DreamPlace::native::detail
