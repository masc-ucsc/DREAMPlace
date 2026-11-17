/**
 * @file   sp_legalize.h
 * @brief  Sequence-pair macro legalization: the fallback when the Hannan grid
 * search and the separation LP cannot find a legal macro pose.
 *
 * The Hannan search places macros greedily one at a time and the LP keeps the
 * left/right/above/below relations of the global placement, so a dense macro
 * set (padded macros filling 80%+ of the core) can be packable and still fail
 * both. Here the relations themselves are searched, as floorplanners do
 * (Murata's sequence pair, the annealing of Parquet and OpenROAD's mpl):
 *
 *  - the pair (G+, G-) starts from the global placement (G- by x+y, G+ by x-y),
 *    so the first relations are the placer's;
 *  - a pair is scored by longest paths on both axes inside the core: fixed
 *    obstacles are pinned nodes, the overflow (lower bound above upper bound)
 *    is the infeasibility, and each macro then goes to its global-placement
 *    target clamped into its feasible window (a feasible pair always yields a
 *    legal, grid-aligned pose);
 *  - annealing swaps entries of G+, G- or both and minimizes overflow first,
 *    then the weighted displacement from the global placement.
 *
 * Sizes are the launcher's (site/row-rounded, already padded by any halo the
 * caller added), so touching padded boxes keep each macro's own halo.
 */
#ifndef DREAMPLACE_MACRO_LEGALIZE_SP_LEGALIZE_H
#define DREAMPLACE_MACRO_LEGALIZE_SP_LEGALIZE_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "utility/src/legalization_db.h"

DREAMPLACE_BEGIN_NAMESPACE

struct SequencePairNode {
  double w = 0, h = 0;        // size (padded, rounded)
  double tx = 0, ty = 0;      // target origin (global placement) or pinned origin
  double weight = 1;
  bool pinned = false;        // fixed obstacle inside the core
};

struct SequencePairResult {
  double overflow = 0;        // total lower-bound-over-upper-bound on both axes, in DBU
  double displacement = 0;    // weighted |x - tx| + |y - ty| over movable nodes
  std::vector<double> x, y;   // origins
  double cost() const { return overflow * 1000 + displacement; }
};

/// Score one sequence pair and compute its pose. `gp`/`gm` are permutations of
/// node indices; `lo`/`hi` are the core bounds; `grid` the site/row pitch from `lo`.
inline SequencePairResult sequencePairEvaluate(const std::vector<SequencePairNode>& nodes,
                                               const std::vector<int>& gp, const std::vector<int>& gm,
                                               double xl, double yl, double xh, double yh,
                                               double site, double row, bool want_pose = true) {
  const int n = int(nodes.size());
  std::vector<int> rank_p(n), rank_m(n);
  for (int k = 0; k < n; ++k) {
    rank_p[gp[k]] = k;
    rank_m[gm[k]] = k;
  }
  SequencePairResult out;
  out.x.assign(n, 0);
  out.y.assign(n, 0);
  std::vector<double> lb(n), ub(n);
  for (int axis = 0; axis < 2; ++axis) {
    const double lo = axis == 0 ? xl : yl, hi = axis == 0 ? xh : yh, g = axis == 0 ? site : row;
    const auto size = [&](int i) { return axis == 0 ? nodes[i].w : nodes[i].h; };
    const auto target = [&](int i) { return axis == 0 ? nodes[i].tx : nodes[i].ty; };
    const auto ceil_g = [&](double v) { return g > 0 ? lo + std::ceil((v - lo) / g - 1e-9) * g : v; };
    const auto floor_g = [&](double v) { return g > 0 ? lo + std::floor((v - lo) / g + 1e-9) * g : v; };
    const auto round_g = [&](double v) { return g > 0 ? lo + std::round((v - lo) / g) * g : v; };
    // j precedes i on this axis: left of (before in both) or below (after in G+, before in G-)
    const auto before = [&](int j, int i) {
      return rank_m[j] < rank_m[i] && (axis == 0 ? rank_p[j] < rank_p[i] : rank_p[j] > rank_p[i]);
    };
    // G- is a topological order of both relations
    for (int k = 0; k < n; ++k) {
      const int i = gm[k];
      double v = lo;
      for (int q = 0; q < k; ++q) {
        const int j = gm[q];
        if (before(j, i)) v = std::max(v, lb[j] + size(j));
      }
      if (nodes[i].pinned) {
        out.overflow += std::max(0.0, v - target(i));
        lb[i] = target(i);
      } else {
        lb[i] = ceil_g(v);
      }
    }
    for (int k = n - 1; k >= 0; --k) {
      const int i = gm[k];
      double v = hi - size(i);
      for (int q = k + 1; q < n; ++q) {
        const int j = gm[q];
        if (before(i, j)) v = std::min(v, ub[j] - size(i));
      }
      if (nodes[i].pinned) {
        out.overflow += std::max(0.0, target(i) - v);
        ub[i] = target(i);
      } else {
        ub[i] = floor_g(v);
        out.overflow += std::max(0.0, lb[i] - ub[i]);
      }
    }
    if (!want_pose && out.overflow > 0) continue;
    auto& pos = axis == 0 ? out.x : out.y;
    for (int k = 0; k < n; ++k) {
      const int i = gm[k];
      if (nodes[i].pinned) {
        pos[i] = target(i);
        continue;
      }
      double v = lo;
      for (int q = 0; q < k; ++q) {
        const int j = gm[q];
        if (before(j, i)) v = std::max(v, pos[j] + size(j));
      }
      v = ceil_g(v);
      // a feasible pair keeps v <= ub[i]: every predecessor sat at or below its own upper bound
      pos[i] = std::clamp(round_g(target(i)), v, std::max(v, ub[i]));
      out.displacement += nodes[i].weight * std::abs(pos[i] - target(i));
    }
  }
  return out;
}

/// Anneal a sequence pair for `nodes`. Returns the best result found (overflow 0 means legal).
inline SequencePairResult sequencePairAnneal(const std::vector<SequencePairNode>& nodes, double xl, double yl,
                                             double xh, double yh, double site, double row,
                                             double seconds = 3.0, unsigned seed = 1) {
  const int n = int(nodes.size());
  std::vector<int> gp(n), gm(n);
  std::iota(gp.begin(), gp.end(), 0);
  std::iota(gm.begin(), gm.end(), 0);
  const auto cx = [&](int i) { return nodes[i].tx + nodes[i].w / 2; };
  const auto cy = [&](int i) { return nodes[i].ty + nodes[i].h / 2; };
  std::stable_sort(gm.begin(), gm.end(), [&](int a, int b) { return cx(a) + cy(a) < cx(b) + cy(b); });
  std::stable_sort(gp.begin(), gp.end(), [&](int a, int b) { return cx(a) - cy(a) < cx(b) - cy(b); });
  auto best = sequencePairEvaluate(nodes, gp, gm, xl, yl, xh, yh, site, row);
  if (n < 2) return best;
  auto best_gp = gp, best_gm = gm;
  double cur = best.cost();
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> pick(0, n - 1);
  std::uniform_real_distribution<double> unit(0, 1);
  const auto move = [&](std::vector<int>& p, std::vector<int>& m) {
    int a = pick(rng), b = pick(rng);
    while (b == a) b = pick(rng);
    const double r = unit(rng);
    if (r < .4) {
      std::swap(p[a], p[b]);
    } else if (r < .8) {
      std::swap(m[a], m[b]);
    } else {  // swap the same two nodes in both sequences
      const int u = p[a], v = p[b];
      std::swap(p[a], p[b]);
      const auto iu = std::find(m.begin(), m.end(), u), iv = std::find(m.begin(), m.end(), v);
      std::iter_swap(iu, iv);
    }
  };
  // initial temperature: the mean uphill step of a few random moves
  double uphill = 0;
  int ups = 0;
  for (int s = 0; s < 64; ++s) {
    auto p = gp, m = gm;
    move(p, m);
    const double d = sequencePairEvaluate(nodes, p, m, xl, yl, xh, yh, site, row).cost() - cur;
    if (d > 0) {
      uphill += d;
      ++ups;
    }
  }
  const double t0 = ups ? uphill / ups : 1.0, t_end = t0 * 1e-4;
  const long moves = std::clamp<long>(long(n) * 4000, 20000, 400000);
  const auto start = std::chrono::steady_clock::now();
  double t = t0;
  const double cool = std::pow(t_end / t0, 1.0 / double(moves));
  for (long it = 0; it < moves; ++it, t *= cool) {
    if ((it & 255) == 0 &&
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() > seconds)
      break;
    auto p = gp, m = gm;
    move(p, m);
    const auto r = sequencePairEvaluate(nodes, p, m, xl, yl, xh, yh, site, row);
    const double c = r.cost(), d = c - cur;
    if (d <= 0 || unit(rng) < std::exp(-d / t)) {
      gp = std::move(p);
      gm = std::move(m);
      cur = c;
      if (c < best.cost()) {
        best = r;
        best_gp = gp;
        best_gm = gm;
      }
    }
  }
  return best;
}

/// Legalize `macros` (movable) around `fixed_macros` with a sequence pair.
/// Writes db.x/db.y of the movable macros only when a legal pose was found.
template <typename T>
bool spLegalizeLauncher(const LegalizationDB<T>& db, const std::vector<int>& macros,
                        const std::vector<int>& fixed_macros, double seconds = 3.0) {
  if (macros.empty()) return true;
  if (macros.size() > 400) return false;  // O(n^2) scoring: dense macro arrays need another packer
  std::vector<SequencePairNode> nodes;
  nodes.reserve(macros.size() + fixed_macros.size());
  for (int id : macros) {
    SequencePairNode s;
    s.w = db.node_size_x[id];
    s.h = db.node_size_y[id];
    s.tx = std::clamp<double>(db.init_x[id], db.xl, std::max<double>(db.xl, db.xh - s.w));
    s.ty = std::clamp<double>(db.init_y[id], db.yl, std::max<double>(db.yl, db.yh - s.h));
    s.weight = db.node_weights ? std::max<double>(db.node_weights[id], 1e-6) : 1;
    nodes.push_back(s);
  }
  for (int id : fixed_macros) {
    // only obstacles inside the core constrain the packing; clip them to it
    const double fx = std::max<double>(db.init_x[id], db.xl), fy = std::max<double>(db.init_y[id], db.yl);
    const double fxh = std::min<double>(db.init_x[id] + db.node_size_x[id], db.xh);
    const double fyh = std::min<double>(db.init_y[id] + db.node_size_y[id], db.yh);
    if (fxh <= fx || fyh <= fy) continue;
    SequencePairNode s;
    s.w = fxh - fx;
    s.h = fyh - fy;
    s.tx = fx;
    s.ty = fy;
    s.pinned = true;
    nodes.push_back(s);
  }
  const auto r = sequencePairAnneal(nodes, db.xl, db.yl, db.xh, db.yh, db.site_width, db.row_height, seconds);
  // a warning: the fallback ran because the Hannan search and the LPs failed
  dreamplacePrint(r.overflow > 0 ? kERROR : kWARN,
                  "Sequence-pair macro legalization (fallback): %lu macros, %lu obstacles, overflow %g, displacement %g\n",
                  macros.size(), nodes.size() - macros.size(), r.overflow, r.displacement);
  if (r.overflow > 0) return false;
  for (size_t k = 0; k < macros.size(); ++k) {
    db.x[macros[k]] = r.x[k];
    db.y[macros[k]] = r.y[k];
  }
  return true;
}

DREAMPLACE_END_NAMESPACE

#endif
