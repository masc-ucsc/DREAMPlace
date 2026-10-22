// BSD-3-Clause. See LICENSE.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

#include "abacus_legalize/src/abacus_legalize_cpu.h"
#include "dreamplace/native/placement.hpp"
#include "greedy_legalize/src/function_cpu.h"
#include "utility/src/detailed_place_db.h"
#include "utility/src/legalization_db.h"

namespace DreamPlace {
template <typename T>
int computeWeightedAverageWirelengthMergedLauncher(const T *, const T *,
                                                   const int *, const int *,
                                                   const unsigned char *, int,
                                                   const T *, T *, T *, T *,
                                                   int);
template <typename T>
int kreorderCPULauncher(DetailedPlaceDB<T> &, int, int, int);
template <typename T> int globalSwapCPULauncher(DetailedPlaceDB<T>, int);
template <typename T>
void independentSetMatchingCPULauncher(DetailedPlaceDB<T>, int, int, int, int);
template <typename T>
bool macroLegalizationLauncher(
    LegalizationDB<T>, const std::vector<MacroAlignConstraint> * = nullptr);
} // namespace DreamPlace
namespace DreamPlace::native::detail {
// Legacy row kernels infer movable identity from an index prefix. Present only
// standard cells in that prefix; even one-row macros remain exact obstacles.
struct CellMovableView {
  Design design;
  std::vector<int> original;
  explicit CellMovableView(const Design &source) {
    design.layout = source.layout;
    for (int i = 0; i < source.num_movable; ++i)
      if (!source.node_macro[i])
        original.push_back(i);
    design.num_movable = original.size();
    for (int i = 0; i < source.num_movable; ++i)
      if (source.node_macro[i])
        original.push_back(i);
    for (size_t i = source.num_movable; i < source.node_x.size(); ++i)
      original.push_back(i);
    const auto copy = [&](auto &to, const auto &from) {
      if (!from.empty())
        for (int i : original)
          to.push_back(from[i]);
    };
    copy(design.node_x, source.node_x);
    copy(design.node_y, source.node_y);
    copy(design.node_size_x, source.node_size_x);
    copy(design.node_size_y, source.node_size_y);
    copy(design.node_macro, source.node_macro);
    copy(design.node_legalize_weights, source.node_legalize_weights);
    copy(design.node_legalize_last, source.node_legalize_last);
    copy(design.node_legalize_anchor_x, source.node_legalize_anchor_x);
    copy(design.node_legalize_anchor_y, source.node_legalize_anchor_y);
    design.pin2node = source.pin2node;
    std::vector<int> inverse(original.size());
    for (size_t i = 0; i < original.size(); ++i)
      inverse[original[i]] = i;
    for (auto &node : design.pin2node)
      node = inverse[node];
    design.pin_offset_x = source.pin_offset_x;
    design.pin_offset_y = source.pin_offset_y;
    design.flat_net2pin = source.flat_net2pin;
    design.flat_net2pin_start = source.flat_net2pin_start;
    design.net_weights = source.net_weights;
    design.net_mask = source.net_mask;
    design.build_maps();
  }
  void writeback(Design &target) const {
    for (int i = 0; i < design.num_movable; ++i) {
      target.node_x[original[i]] = design.node_x[i];
      target.node_y[original[i]] = design.node_y[i];
    }
  }
};
bool has_movable_macro(const Design &d) {
  return !d.node_macro.empty() &&
         std::any_of(d.node_macro.begin(), d.node_macro.begin() + d.num_movable,
                     [](auto value) { return value != 0; });
}

bool kernel_is_legal(const Design &d) {
  const auto &l = d.layout;
  std::vector<int> regions(d.node_x.size(), 0);
  if (!d.node_macro.empty() &&
      std::any_of(d.node_macro.begin(), d.node_macro.begin() + d.num_movable,
                  [](auto macro) { return macro != 0; })) {
    // The legacy row checker has no macro identity. Present movable macros as
    // fixed obstacles, retaining its independent checks on every CORE row/site.
    // Native is_legal separately checks movable macro bounds and exact overlap.
    std::vector<double> x, y, w, h;
    int core_movable = 0;
    const auto append = [&](int i) {
      x.push_back(d.node_x[i]);
      y.push_back(d.node_y[i]);
      w.push_back(d.node_size_x[i]);
      h.push_back(d.node_size_y[i]);
    };
    for (int i = 0; i < d.num_movable; ++i)
      if (!d.node_macro[i]) {
        append(i);
        ++core_movable;
      }
    for (int i = 0; i < d.num_movable; ++i)
      if (d.node_macro[i])
        append(i);
    for (size_t i = d.num_movable; i < d.node_x.size(); ++i)
      append(i);
    return legalityCheckKernelCPU(
        x.data(), y.data(), w.data(), h.data(),
        static_cast<const double *>(nullptr), static_cast<const int *>(nullptr),
        regions.data(), l.xl, l.yl, l.xh, l.yh, l.site_width, l.row_height, 1.0,
        d.node_x.size(), core_movable, 0);
  }
  return legalityCheckKernelCPU(
      d.node_x.data(), d.node_y.data(), d.node_size_x.data(),
      d.node_size_y.data(), static_cast<const double *>(nullptr),
      static_cast<const int *>(nullptr), regions.data(), l.xl, l.yl, l.xh, l.yh,
      l.site_width, l.row_height, 1.0, d.node_x.size(), d.num_movable, 0);
}

Gradient wirelength(const Design &d, double gamma, int threads) {
  const int n = d.node_x.size(), p = d.pin2node.size(),
            nets = d.net_weights.size();
  std::vector<double> x(p), y(p), gx(p), gy(p), wl(nets);
  for (int i = 0; i < p; ++i) {
    x[i] = d.node_x[d.pin2node[i]] + d.pin_offset_x[i];
    y[i] = d.node_y[d.pin2node[i]] + d.pin_offset_y[i];
  }
  double inv = 1 / gamma;
  computeWeightedAverageWirelengthMergedLauncher<double>(
      x.data(), y.data(), d.flat_net2pin.data(), d.flat_net2pin_start.data(),
      d.net_mask.data(), nets, &inv, wl.data(), gx.data(), gy.data(), threads);
  Gradient g;
  g.x.resize(n);
  g.y.resize(n);
  for (int i = 0; i < nets; ++i)
    g.value += wl[i] * d.net_weights[i];
  for (int i = 0; i < p; ++i) {
    double w = d.net_weights[d.pin2net[i]];
    g.x[d.pin2node[i]] += w * gx[i];
    g.y[d.pin2node[i]] += w * gy[i];
  }
  return g;
}

// Priority packing uses ordered free row intervals. Search effort is bounded
// per cell; unsupported or fragmented rows fall back atomically to the original
// legalizer, which still consumes the timing displacement weights.
bool priority_pack(Design &d, int row_budget = 64, int gap_budget = 128) {
  const auto &l = d.layout;
  const int rows = std::llround((l.yh - l.yl) / l.row_height);
  const int64_t sites = std::floor((l.xh - l.xl) / l.site_width + 1e-8);
  if (rows <= 0 || sites <= 0)
    return false;
  using Row = std::map<int64_t, int64_t>;
  std::vector<Row> free(rows);
  for (auto &row : free)
    row.emplace(0, sites);
  auto reserve = [&](int row, int64_t lo, int64_t hi) {
    auto &gaps = free[row];
    auto it = gaps.upper_bound(lo);
    if (it != gaps.begin())
      --it;
    while (it != gaps.end() && it->first < hi) {
      const auto left = it->first, right = it->second;
      if (right <= lo) {
        ++it;
        continue;
      }
      it = gaps.erase(it);
      if (left < lo)
        gaps.emplace(left, lo);
      if (right > hi) {
        gaps.emplace(hi, right);
        break;
      }
    }
  };
  std::vector<int> order;
  for (int i = 0; i < int(d.node_x.size()); ++i) {
    if (i >= d.num_movable || (!d.node_macro.empty() && d.node_macro[i]) ||
        d.node_size_y[i] > l.row_height + 1e-8) {
      // Clamp in floating point: fixed terminals may lie far outside the core.
      const auto grid = [](double value, double limit) {
        return std::clamp(value, 0.0, limit);
      };
      const int lo = int(
          grid(std::floor((d.node_y[i] - l.yl) / l.row_height + 1e-8), rows));
      const int hi = int(grid(
          std::ceil((d.node_y[i] + d.node_size_y[i] - l.yl) / l.row_height -
                    1e-8),
          rows));
      const int64_t left =
          int64_t(grid(std::floor((d.node_x[i] - l.xl) / l.site_width + 1e-8),
                       double(sites)));
      const int64_t right = int64_t(grid(
          std::ceil((d.node_x[i] + d.node_size_x[i] - l.xl) / l.site_width -
                    1e-8),
          double(sites)));
      if (right > left)
        for (int row = lo; row < hi; ++row)
          reserve(row, left, right);
    } else if (d.node_size_x[i] > 0 && d.node_size_y[i] > 0)
      order.push_back(i);
  }
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    const bool last_a =
        !d.node_legalize_last.empty() && d.node_legalize_last[a];
    const bool last_b =
        !d.node_legalize_last.empty() && d.node_legalize_last[b];
    if (last_a != last_b)
      return !last_a;
    return !d.node_legalize_weights.empty() &&
           d.node_legalize_weights[a] > d.node_legalize_weights[b];
  });
  for (int i : order) {
    double best = std::numeric_limits<double>::infinity();
    int best_row = -1;
    int64_t best_site = 0;
    const int64_t width = std::ceil(d.node_size_x[i] / l.site_width - 1e-8);
    const int64_t wanted = std::llround((d.node_x[i] - l.xl) / l.site_width);
    const int near = std::clamp(
        int(std::llround((d.node_y[i] - l.yl) / l.row_height)), 0, rows - 1);
    int rows_visited = 0, gaps_visited = 0;
    for (int delta = 0;
         delta < rows && rows_visited < row_budget && gaps_visited < gap_budget;
         ++delta) {
      if (delta > 1 && (delta - 1) * l.row_height > best)
        break;
      for (int side = 0; side < (delta ? 2 : 1); ++side) {
        const int row = near + (side ? -delta : delta);
        if (row < 0 || row >= rows || rows_visited++ >= row_budget)
          continue;
        const double dy = std::abs(l.yl + row * l.row_height - d.node_y[i]);
        if (dy > best)
          continue;
        const auto &gaps = free[row];
        auto right = gaps.lower_bound(wanted);
        auto left = right;
        bool has_left = left != gaps.begin();
        if (has_left)
          --left;
        while ((has_left || right != gaps.end()) && gaps_visited < gap_budget) {
          const double left_bound =
              has_left ? std::max(0.0, (wanted - left->second) * l.site_width)
                       : std::numeric_limits<double>::infinity();
          const double right_bound =
              right != gaps.end()
                  ? std::max(0.0, (right->first - wanted) * l.site_width)
                  : std::numeric_limits<double>::infinity();
          if (dy + std::min(left_bound, right_bound) > best)
            break;
          const bool take_left = left_bound <= right_bound;
          const auto gap = take_left ? *left : *right;
          ++gaps_visited;
          if (gap.second - gap.first >= width) {
            const auto site = std::clamp(wanted, gap.first, gap.second - width);
            const double cost =
                dy + std::abs(l.xl + site * l.site_width - d.node_x[i]);
            if (cost < best) {
              best = cost;
              best_site = site;
              best_row = row;
            }
          }
          if (take_left) {
            if (left == gaps.begin())
              has_left = false;
            else
              --left;
          } else
            ++right;
        }
      }
    }
    if (best_row < 0)
      return false;
    d.node_x[i] = l.xl + best_site * l.site_width;
    d.node_y[i] = l.yl + best_row * l.row_height;
    reserve(best_row, best_site, best_site + width);
  }
  return is_legal(d);
}

double protected_displacement(const Design &d) {
  if (d.node_legalize_weights.empty() || d.node_legalize_anchor_x.empty())
    return 0;
  double cost = 0;
  for (int i = 0; i < d.num_movable; ++i)
    cost += (d.node_legalize_weights[i] - 1) *
            (std::abs(d.node_x[i] - d.node_legalize_anchor_x[i]) +
             std::abs(d.node_y[i] - d.node_legalize_anchor_y[i]));
  return cost;
}

void legalize(Design &d, int bx, int by, const AlignmentSnap *snap,
              AlignmentReport *report,
              const std::function<void(Design &)> &after_macro_legalize,
              bool macro_only) {
  auto ix = d.node_x, iy = d.node_y;
  if (snap)
    for (const auto &group : snap->groups) {
      auto &anchors = group.axis == Axis::x ? ix : iy;
      for (int i : group.nodes)
        anchors[i] = group.target;
    }
  if (snap)
    for (const auto &target : snap->targets) {
      if (target.node < 0 || target.node >= d.num_movable ||
          d.node_macro.empty() || !d.node_macro[target.node] ||
          !std::isfinite(target.x) || !std::isfinite(target.y))
        throw std::invalid_argument("snap targets must name movable macros");
      ix[target.node] = target.x;
      iy[target.node] = target.y;
    }
  auto weights = d.node_legalize_weights;
  if (weights.empty())
    weights.assign(ix.size(), 1);
  std::vector<int> regions(ix.size(), 0);
  const auto &l = d.layout;
  LegalizationDB<double> db{};
  db.init_x = ix.data();
  db.init_y = iy.data();
  db.x = d.node_x.data();
  db.y = d.node_y.data();
  db.node_size_x = d.node_size_x.data();
  db.node_size_y = d.node_size_y.data();
  std::vector<double> macro_sizes_x, macro_sizes_y;
  if (macro_only) {
    macro_sizes_x = d.node_size_x;
    macro_sizes_y = d.node_size_y;
    for (int i = 0; i < d.num_movable; ++i)
      if (d.node_macro.empty() || !d.node_macro[i])
        macro_sizes_x[i] = macro_sizes_y[i] = 0;
    // The legacy height classifier also recognizes tall standard cells.
    // Explicit macro-only mode must leave those cells untouched as well.
    db.node_size_x = macro_sizes_x.data();
    db.node_size_y = macro_sizes_y.data();
  }
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
  if ((snap || after_macro_legalize || macro_only) && !d.node_macro.empty())
    db.node_macro_mask = d.node_macro.data();
  const auto macro_start = std::chrono::steady_clock::now();
  if (!macroLegalizationLauncher(db, snap ? &snap->constraints : nullptr))
    throw std::runtime_error("DREAMPlace macro legalization failed");
  if (after_macro_legalize) {
    after_macro_legalize(d);
    d.validate();
    // A halo adapter can replace the coordinate/size arrays while translating
    // between physical and padded geometry. Refresh all database pointers.
    db.x = d.node_x.data();
    db.y = d.node_y.data();
    db.node_size_x = d.node_size_x.data();
    db.node_size_y = d.node_size_y.data();
    db.node_macro_mask = d.node_macro.empty() ? nullptr : d.node_macro.data();
  }
  const auto cell_start = std::chrono::steady_clock::now();
  if (report)
    report->macro_legalization_seconds =
        std::chrono::duration<double>(cell_start - macro_start).count();
  if (macro_only)
    return;
  // Greedy's fixed-obstacle bins and row blanks are built from init_x/y.
  // Macros have just moved: reserve their legalized footprints, while keeping
  // the original standard-cell anchors for the displacement objective.
  for (int i = 0; i < d.num_movable; ++i) {
    if (db.is_dummy_fixed(i)) {
      ix[i] = d.node_x[i];
      iy[i] = d.node_y[i];
    }
  }
  if ((snap || after_macro_legalize) && has_movable_macro(d)) {
    CellMovableView cells(d);
    if (cells.design.num_movable)
      legalize(cells.design, bx, by);
    cells.writeback(d);
    if (report)
      report->cell_legalization_seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        cell_start)
              .count();
    return;
  }
  bool packed = false;
  const bool deferred =
      std::any_of(d.node_legalize_last.begin(), d.node_legalize_last.end(),
                  [](auto last) { return last != 0; });
  if (!d.node_legalize_weights.empty() || deferred) {
    const auto x = d.node_x, y = d.node_y;
    packed = priority_pack(d);
    if (!packed) {
      std::copy(x.begin(), x.end(), d.node_x.begin());
      std::copy(y.begin(), y.end(), d.node_y.begin());
      // A bounded wider retry retains functional-first processing when the
      // nearest rows are full. Infeasible fragmentation still uses the original
      // legalizer atomically, before the host's final local tie-site repair.
      if (deferred)
        packed = priority_pack(d, 256, 512);
      if (!packed) {
        std::copy(x.begin(), x.end(), d.node_x.begin());
        std::copy(y.begin(), y.end(), d.node_y.begin());
      }
    }
  }
  if (!packed) {
    greedyLegalizationCPU(db, db.init_x, db.init_y, db.node_size_x,
                          db.node_size_y, db.x, db.y, db.xl, db.yl, db.xh,
                          db.yh, db.site_width, db.row_height, bx, by,
                          db.num_nodes, db.num_movable_nodes);
  }
  // Abacus operates within the rows established by greedy legalization.
  if (is_legal(d)) {
    auto gx = d.node_x, gy = d.node_y;
    const double protected_before = protected_displacement(d);
    abacusLegalizationCPU(
        ix.data(), iy.data(), db.node_size_x, db.node_size_y,
        d.node_legalize_weights.empty() ? nullptr : weights.data(), db.x, db.y,
        db.xl, db.yl, db.xh, db.yh, db.site_width, db.row_height, 1, by,
        db.num_nodes, db.num_movable_nodes);
    bool moved_macro = false;
    if ((snap || after_macro_legalize) && !d.node_macro.empty())
      for (int i = 0; i < d.num_movable; ++i)
        moved_macro |=
            d.node_macro[i] && (d.node_x[i] != gx[i] || d.node_y[i] != gy[i]);
    if (!is_legal(d) || moved_macro ||
        protected_displacement(d) > protected_before + 1e-8) {
      d.node_x = std::move(gx);
      d.node_y = std::move(gy);
    }
  }
  if (report)
    report->cell_legalization_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      cell_start)
            .count();
}

void detailed_place(Design &d, int bx, int by, int threads,
                    bool preserve_macros, int iterations,
                    std::vector<DetailedSample> *history) {
  if (preserve_macros && has_movable_macro(d)) {
    CellMovableView cells(d);
    if (cells.design.num_movable)
      detailed_place(cells.design, bx, by, threads, false, iterations, history);
    cells.writeback(d);
    return;
  }
  auto ix = d.node_x, iy = d.node_y;
  std::vector<int> regions(ix.size(), 0);
  const auto &l = d.layout;
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
  db.net_weights = d.net_weights.data();
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
  auto run = [&](const char *kernel, auto pass) {
    auto x = d.node_x, y = d.node_y;
    double before = hpwl(d);
    const double protected_before = protected_displacement(d);
    pass();
    bool moved_macro = false;
    if (preserve_macros && !d.node_macro.empty())
      for (int i = 0; i < d.num_movable; ++i)
        moved_macro |=
            d.node_macro[i] && (d.node_x[i] != x[i] || d.node_y[i] != y[i]);
    std::string legality_issue;
    const bool legal = is_legal(d, history ? &legality_issue : nullptr);
    const double after = hpwl(d);
    const bool rejected = !legal || after > before + 1e-8 || moved_macro ||
                          protected_displacement(d) > protected_before + 1e-8;
    if (history)
      history->push_back({kernel, before, after, legal, !moved_macro, !rejected,
                          legality_issue});
    if (rejected) {
      // Copy into existing arrays: db contains non-owning pointers.
      std::copy(x.begin(), x.end(), d.node_x.begin());
      std::copy(y.begin(), y.end(), d.node_y.begin());
    }
  };
  run("k_reorder", [&] { kreorderCPULauncher(db, 4, iterations, threads); });
  run("independent_set_matching", [&] {
    independentSetMatchingCPULauncher(db, 16, 16, iterations, threads);
  });
  run("global_swap", [&] { globalSwapCPULauncher(db, iterations); });
}
} // namespace DreamPlace::native::detail
