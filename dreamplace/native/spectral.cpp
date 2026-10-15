// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace DreamPlace::native {
PoissonSolver::PoissonSolver(int nx, int ny, double width, double height)
    : nx_(nx), ny_(ny), width_(width), height_(height) {
  if (nx < 2 || ny < 2 || nx > 1024 || ny > 1024 || !std::isfinite(width) || !std::isfinite(height) || width <= 0 ||
      height <= 0)
    throw std::invalid_argument("invalid Poisson grid");
  auto basis = [](int n, auto& c, auto& s) {
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

void PoissonSolver::solve(const std::vector<double>& rho, std::vector<double>& phi, std::vector<double>& fx,
                          std::vector<double>& fy) const {
  if (rho.size() != static_cast<size_t>(nx_ * ny_)) throw std::invalid_argument("density grid size");
  const int nx = nx_, ny = ny_;
  std::vector<double> temp(nx * ny), coeff(nx * ny);
  for (int k = 0; k < nx; ++k)
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x) temp[k * ny + y] += cx_[k * nx + x] * rho[x * ny + y];
  for (int k = 0; k < nx; ++k)
    for (int l = 0; l < ny; ++l) {
      double v = 0;
      for (int y = 0; y < ny; ++y) v += cy_[l * ny + y] * temp[k * ny + y];
      double wx = std::numbers::pi * k / width_, wy = std::numbers::pi * l / height_;
      coeff[k * ny + l] = (k || l) ? v / (wx * wx + wy * wy) : 0;
    }
  // Inverse mixed cosine/sine transforms. E = -grad(phi).
  std::vector<double> vx(nx * ny), vy(nx * ny);
  std::fill(temp.begin(), temp.end(), 0);
  for (int k = 0; k < nx; ++k)
    for (int y = 0; y < ny; ++y)
      for (int l = 0; l < ny; ++l) {
        double v = coeff[k * ny + l];
        temp[k * ny + y] += cy_[l * ny + y] * v;
        vy[k * ny + y] += sy_[l * ny + y] * v * std::numbers::pi * l / height_;
      }
  phi.assign(nx * ny, 0);
  fx.assign(nx * ny, 0);
  fy.assign(nx * ny, 0);
  for (int x = 0; x < nx; ++x)
    for (int y = 0; y < ny; ++y)
      for (int k = 0; k < nx; ++k) {
        int i = x * ny + y, j = k * ny + y;
        phi[i] += cx_[k * nx + x] * temp[j];
        fx[i] += sx_[k * nx + x] * temp[j] * std::numbers::pi * k / width_;
        fy[i] += cx_[k * nx + x] * vy[j];
      }
}
}  // namespace DreamPlace::native
