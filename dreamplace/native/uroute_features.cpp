// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include "dreamplace/native/uroute.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace DreamPlace::native {
namespace {
// Match Python binary64 arithmetic even on targets that contract multiply-add.
// Feature generation is not a neural-network hot loop. Each product rounds
// before the following addition, independent of compiler FMA contraction.
double multiply(double a, double b) {
  volatile double rounded = a * b;
  return rounded;
}
void grid(int x, int y) {
  if (x < 1 || y < 1 || x > 1024 || y > 1024)
    throw std::invalid_argument("URoute grid dimensions must be in 1..1024");
}
using detail::clamped_bin;
} // namespace
std::vector<double> uroute_resample(std::span<const double> input, int fx,
                                    int fy, int tx, int ty) {
  grid(fx, fy);
  grid(tx, ty);
  if (input.size() != size_t(fx) * fy)
    throw std::invalid_argument("URoute resample shape mismatch");
  for (double v : input)
    if (!std::isfinite(v))
      throw std::invalid_argument("nonfinite URoute map");
  if (fx == tx && fy == ty)
    return {input.begin(), input.end()};
  std::vector<double> result(size_t(tx) * ty);
  // Integer overlap numerators avoid platform-dependent boundary rounding.
  // Target tile area in this common coordinate system is fx*fy.
  for (int x = 0; x < tx; ++x) {
    const int x0 = x * fx / tx, x1 = ((x + 1) * fx + tx - 1) / tx;
    for (int y = 0; y < ty; ++y) {
      const int y0 = y * fy / ty, y1 = ((y + 1) * fy + ty - 1) / ty;
      double sum = 0;
      for (int i = x0; i < x1; ++i) {
        const int a =
            std::min((x + 1) * fx, (i + 1) * tx) - std::max(x * fx, i * tx);
        for (int j = y0; j < y1; ++j) {
          const int b =
              std::min((y + 1) * fy, (j + 1) * ty) - std::max(y * fy, j * ty);
          sum += multiply(input[size_t(i) * fy + j], double(a) * b);
        }
      }
      result[size_t(x) * ty + y] = sum / (double(fx) * fy);
    }
  }
  return result;
}

UrouteMaps uroute_features(const Design &d, int nx, int ny, int ox, int oy) {
  d.validate();
  grid(nx, ny);
  grid(ox, oy);
  const auto &l = d.layout;
  const double dx = (l.xh - l.xl) / nx, dy = (l.yh - l.yl) / ny;
  if (!std::isfinite(dx * dy) || dx * dy <= 0)
    throw std::invalid_argument("invalid URoute core area");
  const size_t area = size_t(nx) * ny;
  std::array<std::vector<double>, uroute_feature_channels> raw;
  for (auto &v : raw)
    v.resize(area);
  auto box = [&](double xl, double yl, double xh, double yh, auto visit) {
    if (!std::isfinite(xl) || !std::isfinite(yl) || !std::isfinite(xh) ||
        !std::isfinite(yh))
      throw std::invalid_argument("nonfinite URoute physical box");
    if (xh <= l.xl || yh <= l.yl || xl >= l.xh || yl >= l.yh || xh <= xl ||
        yh <= yl)
      return;
    for (int x = clamped_bin(xl, l.xl, dx, nx);
         x <= clamped_bin(xh, l.xl, dx, nx); ++x) {
      const double a = std::max(0., std::min(xh, l.xl + multiply(x + 1, dx)) -
                                        std::max(xl, l.xl + multiply(x, dx)));
      for (int y = clamped_bin(yl, l.yl, dy, ny);
           y <= clamped_bin(yh, l.yl, dy, ny); ++y) {
        const double b = std::max(0., std::min(yh, l.yl + multiply(y + 1, dy)) -
                                          std::max(yl, l.yl + multiply(y, dy)));
        if (a > 0 && b > 0)
          visit(size_t(x) * ny + y, multiply(a, b));
      }
    }
  };
  for (size_t n = 0; n < d.net_weights.size(); ++n) {
    if (!d.net_mask[n])
      continue;
    const int begin = d.flat_net2pin_start[n],
              end = d.flat_net2pin_start[n + 1];
    if (end - begin < 2)
      continue;
    double xl = std::numeric_limits<double>::infinity(), yl = xl, xh = -xl,
           yh = -xl;
    for (int k = begin; k < end; ++k) {
      const int p = d.flat_net2pin[k], node = d.pin2node[p];
      const double x = d.node_x[node] + d.pin_offset_x[p],
                   y = d.node_y[node] + d.pin_offset_y[p];
      xl = std::min(xl, x);
      xh = std::max(xh, x);
      yl = std::min(yl, y);
      yh = std::max(yh, y);
    }
    const double w = xh - xl, h = yh - yl;
    if (!std::isfinite(w) || !std::isfinite(h))
      throw std::invalid_argument("URoute net bbox overflow");
    if (w == 0 || h == 0)
      continue;
    const double demand = 1 / w + 1 / h;
    box(xl, yl, xh, yh, [&](size_t i, double a) {
      raw[0][i] += multiply(demand, a / (dx * dy));
    });
    for (int k = begin; k < end; ++k) {
      const int p = d.flat_net2pin[k], node = d.pin2node[p];
      const double x = d.node_x[node] + d.pin_offset_x[p],
                   y = d.node_y[node] + d.pin_offset_y[p];
      // Include boundary IOs in the adjacent core tile; exclude outside pins.
      if (x >= l.xl && x <= l.xh && y >= l.yl && y <= l.yh)
        raw[1][size_t(clamped_bin(x, l.xl, dx, nx)) * ny +
               clamped_bin(y, l.yl, dy, ny)] += demand;
    }
  }
  for (size_t i = 0; i < d.node_x.size(); ++i) {
    if (d.placement_mask(i) != PlacementMask::none)
      continue;
    const double px = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
    const double py = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
    const double xl = d.node_x[i] + px, yl = d.node_y[i] + py;
    const double xh = d.node_x[i] + d.node_size_x[i] - px,
                 yh = d.node_y[i] + d.node_size_y[i] - py;
    const bool macro = !d.node_macro.empty() && d.node_macro[i];
    box(xl, yl, xh, yh, [&](size_t j, double a) {
      if (macro)
        raw[2][j] = 1;
      else
        raw[3][j] += a;
    });
  }
  UrouteMaps result{uroute_feature_channels, ox, oy, true, {}, {}};
  result.values.reserve(size_t(uroute_feature_channels) * ox * oy);
  for (const auto &channel : raw) {
    const auto values = uroute_resample(channel, nx, ny, ox, oy);
    const auto [low, high] = std::minmax_element(values.begin(), values.end());
    if (!std::isfinite(*low) || !std::isfinite(*high))
      throw std::invalid_argument("URoute feature accumulation overflow");
    result.ranges.push_back({*low, *high});
    for (double v : values)
      result.values.push_back(
          *high == *low ? 0.f
                        : static_cast<float>((v - *low) / (*high - *low)));
  }
  return result;
}
} // namespace DreamPlace::native
