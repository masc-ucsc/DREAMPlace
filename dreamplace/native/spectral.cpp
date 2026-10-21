// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace DreamPlace::native {
namespace {
using Complex = std::complex<double>;
bool power_two(int n) { return n > 0 && (n & (n - 1)) == 0; }
void fft(std::vector<Complex> &a, bool inverse) {
  const int n = a.size();
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j)
      std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len *= 2) {
    const Complex root =
        std::polar(1., (inverse ? 2 : -2) * std::numbers::pi / len);
    for (int i = 0; i < n; i += len) {
      Complex w = 1;
      for (int j = 0; j < len / 2; ++j) {
        const Complex u = a[i + j], v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= root;
      }
    }
  }
  if (inverse)
    for (auto &x : a)
      x /= n;
}
// Orthonormal DCT-II/III through an even extension. Scratch is per solve.
void cosine(std::vector<double> &line, bool inverse, bool sine = false) {
  const int n = line.size();
  std::vector<Complex> a(2 * n);
  if (!inverse) {
    for (int j = 0; j < n; ++j)
      a[j] = a[2 * n - 1 - j] = line[j];
    fft(a, false);
    for (int k = 0; k < n; ++k)
      line[k] =
          .5 * (a[k] * std::polar(1., -std::numbers::pi * k / (2 * n))).real() *
          std::sqrt((k ? 2. : 1.) / n);
  } else {
    if (sine) {
      std::reverse(line.begin() + 1, line.end());
      line[0] = 0;
    }
    for (int k = 0; k < n; ++k) {
      a[k] = 2 * line[k] / std::sqrt((k ? 2. : 1.) / n) *
             std::polar(1., std::numbers::pi * k / (2 * n));
      if (k)
        a[2 * n - k] = std::conj(a[k]);
    }
    fft(a, true);
    for (int j = 0; j < n; ++j)
      line[j] = a[j].real() * (sine && (j & 1) ? -1 : 1);
  }
}
void fast_solve(int nx, int ny, double width, double height, double screening,
                const std::vector<double> &rho, std::vector<double> &phi,
                std::vector<double> &fx, std::vector<double> &fy) {
  auto coeff = rho;
  std::vector<double> line(nx);
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x)
      line[x] = coeff[x * ny + y];
    cosine(line, false);
    for (int x = 0; x < nx; ++x)
      coeff[x * ny + y] = line[x];
  }
  line.resize(ny);
  const double mass = screening > 0 ? 1 / (screening * screening) : 0;
  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y)
      line[y] = coeff[x * ny + y];
    cosine(line, false);
    const double wx = std::numbers::pi * x / width;
    for (int y = 0; y < ny; ++y) {
      const double wy = std::numbers::pi * y / height;
      coeff[x * ny + y] = (x || y) ? line[y] / (wx * wx + wy * wy + mass) : 0;
    }
  }
  auto temp = coeff, vy = coeff;
  for (int x = 0; x < nx; ++x) {
    for (int y = 0; y < ny; ++y)
      line[y] = coeff[x * ny + y];
    cosine(line, true);
    for (int y = 0; y < ny; ++y)
      temp[x * ny + y] = line[y];
    for (int y = 0; y < ny; ++y)
      line[y] = coeff[x * ny + y] * std::numbers::pi * y / height;
    cosine(line, true, true);
    for (int y = 0; y < ny; ++y)
      vy[x * ny + y] = line[y];
  }
  phi.resize(nx * ny);
  fx.resize(nx * ny);
  fy.resize(nx * ny);
  line.resize(nx);
  for (int y = 0; y < ny; ++y) {
    for (int x = 0; x < nx; ++x)
      line[x] = temp[x * ny + y];
    cosine(line, true);
    for (int x = 0; x < nx; ++x)
      phi[x * ny + y] = line[x];
    for (int x = 0; x < nx; ++x)
      line[x] = temp[x * ny + y] * std::numbers::pi * x / width;
    cosine(line, true, true);
    for (int x = 0; x < nx; ++x)
      fx[x * ny + y] = line[x];
    for (int x = 0; x < nx; ++x)
      line[x] = vy[x * ny + y];
    cosine(line, true);
    for (int x = 0; x < nx; ++x)
      fy[x * ny + y] = line[x];
  }
}
} // namespace
PoissonSolver::PoissonSolver(int nx, int ny, double width, double height,
                             bool fast, double screening_length)
    : nx_(nx), ny_(ny), width_(width), height_(height),
      fast_(fast && power_two(nx) && power_two(ny)),
      screening_(screening_length) {
  if (nx < 2 || ny < 2 || nx > 1024 || ny > 1024 || !std::isfinite(width) ||
      !std::isfinite(height) || width <= 0 || height <= 0)
    throw std::invalid_argument("invalid Poisson grid");
  if (fast_)
    return;
  auto basis = [](int n, auto &c, auto &s) {
    c.resize(n * n);
    s.resize(n * n);
    for (int k = 0; k < n; ++k) {
      double norm = std::sqrt((k ? 2.0 : 1.0) / n);
      for (int x = 0; x < n; ++x) {
        double a = std::numbers::pi * k * (x + 0.5) / n;
        c[k * n + x] = norm * std::cos(a);
        s[k * n + x] = norm * std::sin(a);
      }
    }
  };
  basis(nx, cx_, sx_);
  basis(ny, cy_, sy_);
}

void PoissonSolver::solve(const std::vector<double> &rho,
                          std::vector<double> &phi, std::vector<double> &fx,
                          std::vector<double> &fy) const {
  if (rho.size() != static_cast<size_t>(nx_ * ny_))
    throw std::invalid_argument("density grid size");
  if (fast_) {
    fast_solve(nx_, ny_, width_, height_, screening_, rho, phi, fx, fy);
    return;
  }
  const int nx = nx_, ny = ny_;
  std::vector<double> temp(nx * ny), coeff(nx * ny);
  // Accumulate each output in the original x order, with contiguous y lanes
  // in the inner loop. This exposes independent outputs to SIMD without
  // reassociating a floating-point reduction.
  for (int k = 0; k < nx; ++k)
    for (int x = 0; x < nx; ++x)
      for (int y = 0; y < ny; ++y)
        temp[k * ny + y] += cx_[k * nx + x] * rho[x * ny + y];
  for (int k = 0; k < nx; ++k)
    for (int l = 0; l < ny; ++l) {
      double v = 0;
      for (int y = 0; y < ny; ++y)
        v += cy_[l * ny + y] * temp[k * ny + y];
      double wx = std::numbers::pi * k / width_,
             wy = std::numbers::pi * l / height_;
      coeff[k * ny + l] =
          (k || l) ? v / (wx * wx + wy * wy +
                          (screening_ > 0 ? 1 / (screening_ * screening_) : 0))
                   : 0;
    }
  // Inverse mixed cosine/sine transforms. E = -grad(phi).
  std::vector<double> vy(nx * ny);
  std::fill(temp.begin(), temp.end(), 0);
  for (int k = 0; k < nx; ++k)
    for (int l = 0; l < ny; ++l)
      for (int y = 0; y < ny; ++y) {
        double v = coeff[k * ny + l];
        temp[k * ny + y] += cy_[l * ny + y] * v;
        vy[k * ny + y] += sy_[l * ny + y] * v * std::numbers::pi * l / height_;
      }
  phi.assign(nx * ny, 0);
  fx.assign(nx * ny, 0);
  fy.assign(nx * ny, 0);
  for (int x = 0; x < nx; ++x)
    for (int k = 0; k < nx; ++k)
      for (int y = 0; y < ny; ++y) {
        int i = x * ny + y, j = k * ny + y;
        phi[i] += cx_[k * nx + x] * temp[j];
        fx[i] += sx_[k * nx + x] * temp[j] * std::numbers::pi * k / width_;
        fy[i] += cx_[k * nx + x] * vy[j];
      }
}
} // namespace DreamPlace::native
