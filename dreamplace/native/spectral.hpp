// BSD-3-Clause. See LICENSE.
#pragma once
#include <vector>

namespace DreamPlace::native {
// Separable orthonormal DCT-II / DCT-III Poisson solve with Neumann boundaries.
// Plain CPU transform, O(nx*ny*(nx+ny)); caches bases, no FFT/torch dependency.
class PoissonSolver {
 public:
  PoissonSolver(int nx, int ny, double width, double height);
  void solve(const std::vector<double>& density, std::vector<double>& potential, std::vector<double>& field_x,
             std::vector<double>& field_y) const;

 private:
  int nx_, ny_;
  double width_, height_;
  std::vector<double> cx_, cy_, sx_, sy_;
};
}  // namespace DreamPlace::native
