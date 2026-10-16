// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/diff_lp.hpp"
#include <algorithm>
#include <cmath>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <limits>
#include <numeric>
namespace DreamPlace::native::macro {
namespace {
std::string message(const std::vector<std::string> &ids) {
  std::string out = "infeasible H2D hard constraints:";
  for (const auto &id : ids)
    out += " " + id;
  return out;
}
void require(bool ok) {
  if (!ok)
    throw std::invalid_argument(
        "invalid or numerically unsafe H2D difference LP");
}
} // namespace
Infeasible::Infeasible(std::vector<std::string> ids)
    : std::runtime_error(message(ids)), cycle(std::move(ids)) {}
int DifferenceLP::variable(Int lo, Int hi, double coefficient, std::string id) {
  int n = variables.size();
  variables.push_back({lo, hi, coefficient, std::move(id)});
  return n;
}
void DifferenceLP::constrain(int u, int v, Int bound, double penalty,
                             std::string id) {
  arcs.push_back({u, v, bound, penalty, std::move(id)});
}
void DifferenceLP::anchor(int v, Int x, double w) {
  displacement.push_back({v, x, w});
}
DifferenceLP::Solution DifferenceLP::solve() const {
  const int n = variables.size(), origin = n;
  require(n <= 100000 && arcs.size() <= 10000000);
  auto edges = arcs;
  double finite_capacity = 0, max_cost = 1;
  auto coordinate = [&](Int x) {
    require(x >= -1000000000LL && x <= 1000000000LL);
    max_cost = std::max(max_cost, std::abs(double(x)));
  };
  for (int i = 0; i < n; ++i) {
    const auto &v = variables[i];
    coordinate(v.lower);
    coordinate(v.upper);
    require(v.lower <= v.upper && std::isfinite(v.coefficient) &&
            std::abs(v.coefficient) <= 1e9);
    finite_capacity += std::abs(v.coefficient);
    edges.push_back({-1, i, v.upper, 0, v.id + ".upper"});
    edges.push_back({i, -1, -v.lower, 0, v.id + ".lower"});
  }
  for (const auto &a : edges) {
    require(a.u >= -1 && a.u < n && a.v >= -1 && a.v < n &&
            std::isfinite(a.penalty) && a.penalty >= 0 && a.penalty <= 1e9);
    coordinate(a.bound);
    finite_capacity += a.penalty;
  }
  for (const auto &a : displacement) {
    require(a.variable >= 0 && a.variable < n && std::isfinite(a.weight) &&
            a.weight >= 0 && a.weight <= 1e9);
    coordinate(a.anchor);
    finite_capacity += 2 * a.weight;
  }
  // Diagnose hard negative cycles, including fixed/core bounds, before solving
  // a dual circulation. Soft arcs cannot make the primal infeasible.
  std::vector<Int> distance(n + 1);
  std::vector<int> previous(n + 1, -1);
  int changed = -1;
  for (int iter = 0; iter <= n; ++iter) {
    changed = -1;
    for (int e = 0; e < int(edges.size()); ++e) {
      const auto &a = edges[e];
      if (a.penalty)
        continue;
      int u = a.u < 0 ? origin : a.u, v = a.v < 0 ? origin : a.v;
      if (distance[v] > distance[u] + a.bound) {
        distance[v] = distance[u] + a.bound;
        previous[v] = e;
        changed = v;
      }
    }
    if (changed < 0)
      break;
  }
  if (changed >= 0) {
    int v = changed;
    for (int i = 0; i <= n; ++i) {
      int e = previous[v];
      require(e >= 0);
      v = edges[e].u < 0 ? origin : edges[e].u;
    }
    int first = v;
    std::vector<std::string> ids;
    do {
      int e = previous[v];
      ids.push_back(edges[e].id);
      v = edges[e].u < 0 ? origin : edges[e].u;
    } while (v != first && ids.size() <= size_t(n + 2));
    throw Infeasible(std::move(ids));
  }
  const double safe = double(std::numeric_limits<Int>::max() / 64);
  require(std::isfinite(finite_capacity) && finite_capacity * max_cost < safe);
  const double scale = std::max(
      1., std::min(1024., std::floor(safe / std::max(1., finite_capacity *
                                                             max_cost))));
  auto quantize = [&](double w) -> Int {
    if (w == 0)
      return 0;
    double q = std::round(std::abs(w) * scale);
    require(std::isfinite(q) && q < safe);
    return (w < 0 ? -1 : 1) * static_cast<Int>(std::max(1., q));
  };
  lemon::ListDigraph graph;
  std::vector<lemon::ListDigraph::Node> nodes;
  for (int i = 0; i <= n; ++i)
    nodes.push_back(graph.addNode());
  lemon::ListDigraph::ArcMap<Int> costs(graph), caps(graph);
  lemon::ListDigraph::NodeMap<Int> supply(graph);
  auto arc = [&](int u, int v, Int cost, Int capacity) {
    auto e = graph.addArc(nodes[u < 0 ? origin : u], nodes[v < 0 ? origin : v]);
    costs[e] = cost;
    caps[e] = capacity;
  };
  for (const auto &a : edges)
    arc(a.u, a.v, a.bound,
        a.penalty ? quantize(a.penalty) : std::numeric_limits<Int>::max());
  for (const auto &a : displacement)
    if (a.weight) {
      arc(-1, a.variable, a.anchor, quantize(a.weight));
      arc(a.variable, -1, -a.anchor, quantize(a.weight));
    }
  Int total = 0;
  for (int i = 0; i < n; ++i) {
    Int q = quantize(variables[i].coefficient);
    require(std::abs(double(total)) + std::abs(double(q)) < safe);
    supply[nodes[i]] = q;
    total += q;
  }
  supply[nodes[origin]] = -total;
  using Solver = lemon::NetworkSimplex<lemon::ListDigraph, Int, Int>;
  Solver solver(graph);
  if (solver.upperMap(caps).costMap(costs).supplyMap(supply).run() !=
      Solver::OPTIMAL)
    throw std::runtime_error("H2D network simplex did not find an optimum");
  Solution out;
  out.x.resize(n);
  for (int i = 0; i < n; ++i) {
    out.x[i] = solver.potential(nodes[i]) - solver.potential(nodes[origin]);
    require(out.x[i] >= variables[i].lower && out.x[i] <= variables[i].upper);
    out.objective += variables[i].coefficient * out.x[i];
  }
  for (const auto &a : displacement)
    out.objective += a.weight * std::abs(double(out.x[a.variable] - a.anchor));
  for (const auto &a : arcs) {
    Int u = a.u < 0 ? 0 : out.x[a.u], v = a.v < 0 ? 0 : out.x[a.v];
    Int amount = std::max<Int>(0, v - u - a.bound);
    if (amount) {
      if (!a.penalty)
        throw std::runtime_error("H2D hard LP constraint violated");
      out.slack.push_back({a.id, amount, a.penalty});
      out.slack_cost += a.penalty * amount;
    }
  }
  out.objective += out.slack_cost;
  return out;
}
} // namespace DreamPlace::native::macro
