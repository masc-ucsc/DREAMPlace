// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>

#include "dreamplace/native/placement.hpp"
#include "electric_potential/src/electric_density_map.cpp"

namespace DreamPlace::native::detail {
namespace {
// Accumulates into a per-worker copy of the map. Floating-point atomic adds into
// one shared map make the sums depend on thread arrival order, so identical
// runs would drift apart through the GP iterations. With a static schedule each
// worker sees a fixed node range, and the copies are merged in worker order.
struct WorkerSlotAdd {
  using type = double;
  double* base;
  std::vector<std::vector<double>>* slots;
  template <typename V> void operator()(type* dst, V v) const {
    int worker = 0;
#ifdef _OPENMP
    worker = omp_get_thread_num();
#endif
    (*slots)[worker][dst - base] += v;
  }
};
}  // namespace

void density_map(const Design& d, int bx, int by, int threads, const std::vector<double>& sx,
                 const std::vector<double>& sy, const std::vector<double>& ox, const std::vector<double>& oy,
                 const std::vector<double>& ratio, std::vector<double>& map) {
  const auto& l = d.layout;
  map.assign(bx * by, 0);
  std::vector<double> cx(bx), cy(by);
  double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  for (int i = 0; i < bx; ++i) cx[i] = l.xl + (i + 0.5) * dx;
  for (int i = 0; i < by; ++i) cy[i] = l.yl + (i + 0.5) * dy;
#ifdef _OPENMP
  // Bound the private copies (256 MiB) for very fine grids and many workers.
  const size_t bins = static_cast<size_t>(bx) * by;
  const int workers = std::min<size_t>(threads, std::max<size_t>(1, (size_t{1} << 25) / bins));
  if (workers > 1) {
    std::vector<std::vector<double>> slots(workers, std::vector<double>(bins, 0.));
    computeTriangleDensityMapLauncher<double, WorkerSlotAdd>(
        d.node_x.data(), d.node_y.data(), sx.data(), sy.data(), ox.data(), oy.data(), ratio.data(), cx.data(),
        cy.data(), d.node_x.size(), bx, by, l.xl, l.yl, l.xh, l.yh, dx, dy, workers,
        WorkerSlotAdd{map.data(), &slots}, map.data());
    for (const auto& slot : slots)
      for (size_t i = 0; i < bins; ++i) map[i] += slot[i];
    return;
  }
  threads = workers;
#endif
  computeTriangleDensityMapLauncher<double, AtomicAdd<double>>(
      d.node_x.data(), d.node_y.data(), sx.data(), sy.data(), ox.data(), oy.data(), ratio.data(), cx.data(), cy.data(),
      d.node_x.size(), bx, by, l.xl, l.yl, l.xh, l.yh, dx, dy, threads, AtomicAdd<double>{}, map.data());
}
}  // namespace DreamPlace::native::detail
