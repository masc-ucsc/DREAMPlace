// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"
#include "dreamplace/native/spectral.hpp"
namespace DreamPlace::native::detail {
// Physical body dimensions exclude the private halo padding applied by place.
double macro_cell_screening_length(const Design &d, double reach);
// Screened Poisson particle-mesh field for real standard cells only. No nets,
// artificial particles or physical footprint inflation. Macro receivers are
// opt-in.
class LocalDensity {
public:
  LocalDensity(const Design &d, int bins, double range, double target = 0,
               int num_threads = 1);
  // macros_only integrates the raw standard-cell field over each movable
  // macro footprint (already halo padded by place). Macros never emit charge.
  std::vector<double> cell_density(const Design &d) const;
  Gradient gradient(const Design &d, bool macros_only = false,
                    const std::vector<double> &source = {}) const;
  // Prefix-sum local slope: no Poisson solve, cells are the only sources.
  Gradient density_pull(const Design &d, const std::vector<double> &source,
                        double range) const;

private:
  int bins_, physical_nodes_;
  double dx_, dy_, target_;
  PoissonSolver solver_;
};
} // namespace DreamPlace::native::detail
