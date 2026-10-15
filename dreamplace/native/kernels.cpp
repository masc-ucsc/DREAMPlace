// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include "abacus_legalize/src/abacus_legalize_cpu.h"
#include "dreamplace/native/placement.hpp"
#include "greedy_legalize/src/function_cpu.h"
#include "utility/src/detailed_place_db.h"
#include "utility/src/legalization_db.h"

namespace DreamPlace {
template <typename T>
int computeWeightedAverageWirelengthMergedLauncher(const T*, const T*, const int*, const int*, const unsigned char*,
                                                   int, const T*, T*, T*, T*, int);
template <typename T>
int kreorderCPULauncher(DetailedPlaceDB<T>&, int, int, int);
template <typename T>
int globalSwapCPULauncher(DetailedPlaceDB<T>, int);
template <typename T>
void independentSetMatchingCPULauncher(DetailedPlaceDB<T>, int, int, int, int);
template <typename T>
bool macroLegalizationLauncher(LegalizationDB<T>);
}  // namespace DreamPlace
namespace DreamPlace::native::detail {
bool kernel_is_legal(const Design& d) {
  const auto& l = d.layout;
  std::vector<int> regions(d.node_x.size(), 0);
  return legalityCheckKernelCPU(d.node_x.data(), d.node_y.data(), d.node_size_x.data(), d.node_size_y.data(),
                                static_cast<const double*>(nullptr), static_cast<const int*>(nullptr), regions.data(),
                                l.xl, l.yl, l.xh, l.yh, l.site_width, l.row_height, 1.0, d.node_x.size(), d.num_movable,
                                0);
}

Gradient wirelength(const Design& d, double gamma, int threads) {
  const int n = d.node_x.size(), p = d.pin2node.size(), nets = d.net_weights.size();
  std::vector<double> x(p), y(p), gx(p), gy(p), wl(nets);
  for (int i = 0; i < p; ++i) {
    x[i] = d.node_x[d.pin2node[i]] + d.pin_offset_x[i];
    y[i] = d.node_y[d.pin2node[i]] + d.pin_offset_y[i];
  }
  double inv = 1 / gamma;
  computeWeightedAverageWirelengthMergedLauncher<double>(x.data(), y.data(), d.flat_net2pin.data(),
                                                         d.flat_net2pin_start.data(), d.net_mask.data(), nets, &inv,
                                                         wl.data(), gx.data(), gy.data(), threads);
  Gradient g;
  g.x.resize(n);
  g.y.resize(n);
  for (int i = 0; i < nets; ++i) g.value += wl[i] * d.net_weights[i];
  for (int i = 0; i < p; ++i) {
    double w = d.net_weights[d.pin2net[i]];
    g.x[d.pin2node[i]] += w * gx[i];
    g.y[d.pin2node[i]] += w * gy[i];
  }
  return g;
}

void legalize(Design& d, int bx, int by) {
  auto ix = d.node_x, iy = d.node_y;
  std::vector<double> weights(ix.size(), 1);
  std::vector<int> regions(ix.size(), 0);
  const auto& l = d.layout;
  LegalizationDB<double> db{};
  db.init_x = ix.data();
  db.init_y = iy.data();
  db.x = d.node_x.data();
  db.y = d.node_y.data();
  db.node_size_x = d.node_size_x.data();
  db.node_size_y = d.node_size_y.data();
  db.node_weights = weights.data();
  db.node2fence_region_map = regions.data();
  db.xl = l.xl;
  db.yl = l.yl;
  db.xh = l.xh;
  db.yh = l.yh;
  db.site_width = l.site_width;
  db.row_height = l.row_height;
  db.num_bins_x = bx;
  db.num_bins_y = by;
  db.bin_size_x = (l.xh - l.xl) / bx;
  db.bin_size_y = (l.yh - l.yl) / by;
  db.num_sites_x = std::round((l.xh - l.xl) / l.site_width);
  db.num_sites_y = std::round((l.yh - l.yl) / l.row_height);
  db.num_nodes = ix.size();
  db.num_movable_nodes = d.num_movable;
  if (!macroLegalizationLauncher(db)) throw std::runtime_error("DREAMPlace macro legalization failed");
  greedyLegalizationCPU(db, db.init_x, db.init_y, db.node_size_x, db.node_size_y, db.x, db.y, db.xl, db.yl, db.xh,
                        db.yh, db.site_width, db.row_height, bx, by, db.num_nodes, db.num_movable_nodes);
  // Abacus operates within the rows established by greedy legalization.
  if (is_legal(d)) {
    auto gx = d.node_x, gy = d.node_y;
    abacusLegalizationCPU(ix.data(), iy.data(), db.node_size_x, db.node_size_y, weights.data(), db.x, db.y, db.xl,
                          db.yl, db.xh, db.yh, db.site_width, db.row_height, 1, by, db.num_nodes, db.num_movable_nodes);
    if (!is_legal(d)) {
      d.node_x = std::move(gx);
      d.node_y = std::move(gy);
    }
  }
}

void detailed_place(Design& d, int bx, int by, int threads) {
  // Legacy detailed kernels optimize unweighted HPWL. Avoid undoing the
  // timing-weighted objective until weighted detailed placement is supported.
  if (std::any_of(d.net_weights.begin(), d.net_weights.end(), [](double w) { return w != 1; })) return;
  auto ix = d.node_x, iy = d.node_y;
  std::vector<int> regions(ix.size(), 0);
  const auto& l = d.layout;
  DetailedPlaceDB<double> db{};
  db.init_x = ix.data();
  db.init_y = iy.data();
  db.x = d.node_x.data();
  db.y = d.node_y.data();
  db.node_size_x = d.node_size_x.data();
  db.node_size_y = d.node_size_y.data();
  db.node2fence_region_map = regions.data();
  db.flat_net2pin_map = d.flat_net2pin.data();
  db.flat_net2pin_start_map = d.flat_net2pin_start.data();
  db.flat_node2pin_map = d.flat_node2pin.data();
  db.flat_node2pin_start_map = d.flat_node2pin_start.data();
  db.pin2node_map = d.pin2node.data();
  db.pin2net_map = d.pin2net.data();
  db.pin_offset_x = d.pin_offset_x.data();
  db.pin_offset_y = d.pin_offset_y.data();
  db.net_mask = d.net_mask.data();
  db.xl = l.xl;
  db.yl = l.yl;
  db.xh = l.xh;
  db.yh = l.yh;
  db.site_width = l.site_width;
  db.row_height = l.row_height;
  db.num_bins_x = bx;
  db.num_bins_y = by;
  db.bin_size_x = (l.xh - l.xl) / bx;
  db.bin_size_y = (l.yh - l.yl) / by;
  db.num_sites_x = std::round((l.xh - l.xl) / l.site_width);
  db.num_sites_y = std::round((l.yh - l.yl) / l.row_height);
  db.num_nodes = ix.size();
  db.num_movable_nodes = d.num_movable;
  db.num_pins = d.pin2node.size();
  db.num_nets = d.net_weights.size();
  // Each pass must preserve legality and not increase HPWL.
  auto run = [&](auto pass) {
    auto x = d.node_x, y = d.node_y;
    double before = hpwl(d);
    pass();
    if (!is_legal(d) || hpwl(d) > before + 1e-8) {
      // Copy into existing arrays: db contains non-owning pointers.
      std::copy(x.begin(), x.end(), d.node_x.begin());
      std::copy(y.begin(), y.end(), d.node_y.begin());
    }
  };
  run([&] { kreorderCPULauncher(db, 4, 2, threads); });
  run([&] { independentSetMatchingCPULauncher(db, 16, 16, 2, threads); });
  run([&] { globalSwapCPULauncher(db, 2); });
}
}  // namespace DreamPlace::native::detail
