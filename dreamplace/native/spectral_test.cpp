// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/spectral.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <utility>
#include <vector>

// Manufactured Neumann solutions exercise both field directions, the DC
// removal, and nonsquare grids without depending on the transform algorithm.
int main() {
  for (bool fast : {false, true})
    for (const auto &shape : {std::pair{7, 11}, std::pair{16, 16},
                              std::pair{64, 32}, std::pair{128, 64}}) {
      const auto [nx, ny] = shape;
      const double width = 13.7, height = 8.2;
      const double ax = std::numbers::pi / width,
                   ay = 2 * std::numbers::pi / height;
      std::vector<double> rho(nx * ny), phi, fx, fy;
      for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y) {
          const double px = width * (x + .5) / nx, py = height * (y + .5) / ny;
          rho[x * ny + y] =
              (ax * ax + ay * ay) * std::cos(ax * px) * std::cos(ay * py) +
              .3 * 9 * ax * ax * std::cos(3 * ax * px) + 2.1;
        }
      DreamPlace::native::PoissonSolver(nx, ny, width, height, fast)
          .solve(rho, phi, fx, fy);
      for (int x = 0; x < nx; ++x)
        for (int y = 0; y < ny; ++y) {
          const int i = x * ny + y;
          const double px = width * (x + .5) / nx, py = height * (y + .5) / ny;
          const double potential = std::cos(ax * px) * std::cos(ay * py) +
                                   .3 * std::cos(3 * ax * px);
          const double field_x = ax * std::sin(ax * px) * std::cos(ay * py) +
                                 .9 * ax * std::sin(3 * ax * px);
          const double field_y = ay * std::cos(ax * px) * std::sin(ay * py);
          if (!std::isfinite(phi[i]) || !std::isfinite(fx[i]) ||
              !std::isfinite(fy[i]) || std::abs(phi[i] - potential) > 1e-10 ||
              std::abs(fx[i] - field_x) > 1e-10 ||
              std::abs(fy[i] - field_y) > 1e-10) {
            std::cerr << "Poisson solution mismatch on " << nx << 'x' << ny
                      << " at " << x << ',' << y << '\n';
            return 1;
          }
        }
    }
}
