// BSD-3-Clause. See LICENSE.
#pragma once
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "utility/src/namespace.h"

DREAMPLACE_BEGIN_NAMESPACE

// The macro LP is a difference-constrained L1 problem, not a general LP:
//   minimize sum_i |x_i - anchor_i| subject to x_v - x_u <= bound.
// Its dual is a min-cost circulation. Each difference constraint is an
// uncapacitated arc u -> v with cost=bound. Two anchor arcs to/from an origin
// have capacity 1 and costs +/-anchor. Optimal node potentials recover x.
// This uses LEMON network simplex directly, avoiding Limbo's parser/model/LP
// wrapper and any GLPK backend. Coordinates and costs are integer DBU.
class MacroDifferenceLP {
 public:
  using Int = std::int64_t;
  struct Constraint {
    int u, v;
    Int bound;
  };
  std::vector<Int> lower, upper, anchor;
  std::vector<Constraint> constraints;

  static Int coordinate(double value) {
    if (!std::isfinite(value) || std::abs(value) > 1e12)
      throw std::invalid_argument("macro LP coordinates must be finite and within 1e12 DBU");
    return static_cast<Int>(value);
  }

  std::vector<Int> solve() const {
    const int n = anchor.size();
    if (anchor.size() > 1000000 || lower.size() != anchor.size() || upper.size() != anchor.size())
      throw std::invalid_argument("invalid macro LP dimensions");
    lemon::ListDigraph graph;
    std::vector<lemon::ListDigraph::Node> nodes;
    for (int i = 0; i <= n; ++i) nodes.push_back(graph.addNode());
    lemon::ListDigraph::ArcMap<Int> costs(graph), caps(graph);
    using Solver = lemon::NetworkSimplex<lemon::ListDigraph, Int, Int>;
    const Int inf = std::numeric_limits<Int>::max();
    auto arc = [&](int u, int v, Int cost, Int cap) {
      auto a = graph.addArc(nodes[u], nodes[v]);
      costs[a] = cost;
      caps[a] = cap;
    };
    for (int i = 0; i < n; ++i) {
      if (lower[i] > upper[i]) throw std::runtime_error("infeasible macro bounds");
      arc(n, i, upper[i], inf);
      arc(i, n, -lower[i], inf);
      arc(n, i, anchor[i], 1);
      arc(i, n, -anchor[i], 1);
    }
    for (auto c : constraints) {
      if (c.u < 0 || c.u >= n || c.v < 0 || c.v >= n) throw std::invalid_argument("invalid macro constraint");
      arc(c.u, c.v, c.bound, inf);
    }
    Solver solver(graph);
    if (solver.upperMap(caps).costMap(costs).run() != Solver::OPTIMAL)
      throw std::runtime_error("macro separation constraints are infeasible");
    std::vector<Int> x(n);
    for (int i = 0; i < n; ++i) {
      x[i] = solver.potential(nodes[i]) - solver.potential(nodes[n]);
      if (x[i] < lower[i] || x[i] > upper[i]) throw std::runtime_error("invalid macro LP potential");
    }
    for (auto c : constraints)
      if (x[c.v] - x[c.u] > c.bound) throw std::runtime_error("macro LP constraint violation");
    return x;
  }
};
DREAMPLACE_END_NAMESPACE
