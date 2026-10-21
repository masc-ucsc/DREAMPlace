// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/spectral.hpp"

#include <cmath>
#include <cstring>
#include <iostream>
#include <numbers>
#include <utility>
#include <vector>

// Manufactured Neumann solutions exercise both field directions, the DC
// removal, and nonsquare grids without depending on the transform algorithm.
int main() {
  for (bool fast : {false, true}) {
    std::vector<std::pair<int, int>> shapes{
        {7, 11}, {16, 16}, {64, 32}, {128, 64}, {128, 128}, {128, 256}};
    if (fast) {
      shapes.emplace_back(512, 128);
      shapes.emplace_back(512, 512);
    }
    for (double screening : {0., 2.3}) {
      for (const auto &shape : shapes) {
        const auto [nx, ny] = shape;
        const double width = 13.7, height = 8.2;
        const double mass = screening > 0 ? 1 / (screening * screening) : 0;
        const double ax = std::numbers::pi / width,
                     ay = 2 * std::numbers::pi / height;
        std::vector<double> rho(nx * ny), phi, fx, fy;
        for (int x = 0; x < nx; ++x) {
          for (int y = 0; y < ny; ++y) {
            const double px = width * (x + .5) / nx,
                         py = height * (y + .5) / ny;
            rho[x * ny + y] =
                (ax * ax + ay * ay + mass) * std::cos(ax * px) *
                    std::cos(ay * py) +
                .3 * (9 * ax * ax + mass) * std::cos(3 * ax * px) + 2.1;
          }
        }
        const DreamPlace::native::PoissonSolver solver(nx, ny, width, height,
                                                       fast, screening);
        solver.solve(rho, phi, fx, fy);
        for (int x = 0; x < nx; ++x) {
          for (int y = 0; y < ny; ++y) {
            const int i = x * ny + y;
            const double px = width * (x + .5) / nx,
                         py = height * (y + .5) / ny;
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
        {
          // Independent rows and columns must preserve every output bit across
          // thread counts, including on nonsquare and fallback dense grids.
          std::vector<double> parallel_phi, parallel_fx, parallel_fy;
          const DreamPlace::native::PoissonSolver parallel(
              nx, ny, width, height, fast, screening, 4);
          parallel.solve(rho, parallel_phi, parallel_fx, parallel_fy);
          const auto same_bits = [](const auto &a, const auto &b) {
            return a.size() == b.size() &&
                   std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) ==
                       0;
          };
          if (!same_bits(phi, parallel_phi) || !same_bits(fx, parallel_fx) ||
              !same_bits(fy, parallel_fy)) {
            std::cerr << "Parallel Poisson solve changed output bits\n";
            return 1;
          }
        }
        const auto old_phi = phi, old_fx = fx, old_fy = fy;
        solver.solve(rho, phi, fx, fy);
        if (phi != old_phi || fx != old_fx || fy != old_fy) {
          std::cerr << "Repeated Poisson solve changed its result\n";
          return 1;
        }
      }
    }
  }
}
