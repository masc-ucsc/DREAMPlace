// BSD-3-Clause. See LICENSE.
#pragma once
#include <vector>

namespace DreamPlace::native {
// Separable orthonormal DCT-II / DCT-III Poisson solve with Neumann boundaries.
// Default: dense CPU transform O(nx*ny*(nx+ny)), preserving reduction order.
// Opt-in fast radix-2 DCT uses O(nx*ny*log(nx*ny)); no external FFT dependency.
class PoissonSolver {
public:
  PoissonSolver(int nx, int ny, double width, double height, bool fast = false,
                double screening_length = 0);
  void solve(const std::vector<double> &density, std::vector<double> &potential,
             std::vector<double> &field_x, std::vector<double> &field_y) const;

private:
  int nx_, ny_;
  double width_, height_;
  bool fast_;
  double screening_;
  std::vector<double> cx_, cy_, sx_, sy_;
};
} // namespace DreamPlace::native
