// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/spectral.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char **argv) {
  const int loops = argc > 1 ? std::atoi(argv[1]) : 500;
  const int nx = argc > 2 ? std::atoi(argv[2]) : 64,
            ny = argc > 3 ? std::atoi(argv[3]) : nx;
  if (loops < 1)
    return 2;
  DreamPlace::native::PoissonSolver solver(nx, ny, 101.7, 83.3);
  std::vector<double> rho(nx * ny), phi, fx, fy;
  for (int i = 0; i < nx * ny; ++i)
    rho[i] = std::sin(i * .73) + std::cos(i * .17);
  const auto start = std::chrono::steady_clock::now();
  double checksum = 0;
  for (int i = 0; i < loops; ++i) {
    solver.solve(rho, phi, fx, fy);
    checksum += phi[i % phi.size()] + fx[i % fx.size()] + fy[i % fy.size()];
  }
  const auto seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
  std::cout << std::setprecision(17) << "seconds=" << seconds
            << " checksum=" << checksum << '\n';
  if (argc > 4) {
    std::ofstream output(argv[4], std::ios::binary);
    for (const auto *values : {&phi, &fx, &fy})
      output.write(reinterpret_cast<const char *>(values->data()),
                   values->size() * sizeof(double));
    if (!output)
      return 1;
  }
}
