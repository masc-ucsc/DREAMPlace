// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/toposearch.hpp"
#include "dreamplace/native/rudy.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <stdexcept>
namespace DreamPlace::native::macro {
SearchResult toposearch(const Design &src, std::span<const NodeInfo> info,
                        const Realization &incumbent, const Parameters &p,
                        const SequencePair &previous,
                        std::function<double(const Design &)> route_estimator) {
  p.validate();
  SearchResult out;
  out.best = incumbent;
  out.sequence = sequence_pair(incumbent.topology);
  if (!p.toposearch || p.search == "none" || incumbent.topology.groups.empty())
    return out;
  struct Costs {
    double wire = 0, displace = 0, order = 0, anchor = 0, slack = 0, route = 0;
  };
  auto physical = [&](const Realization &r) {
    auto a = abstract_design(src, info, r.topology, p);
    Design d = src;
    expand(a, r.topology, d);
    return d;
  };

  const auto &layout = src.layout;
  const double caph = p.route_cap_h > 0 ? p.route_cap_h : 1 / layout.site_width,
               capv = p.route_cap_v > 0 ? p.route_cap_v : 1 / layout.row_height;
  std::map<std::string, std::pair<double, double>> seed;
  for (const auto &g : incumbent.topology.groups)
    seed[g.id] = {g.x, g.y};
  auto cheap = [&](const Realization &r) {
    Costs c;
    c.wire = hpwl(physical(r));
    c.order = r.order_violation;
    c.anchor = r.anchor_deviation;
    c.slack = r.slack_amount;
    for (const auto &g : r.topology.groups)
      c.displace += std::abs(g.x - seed.at(g.id).first) +
                    std::abs(g.y - seed.at(g.id).second);
    return c;
  };
  auto route = [&](const Realization &r) {
    auto design = physical(r);
    return route_estimator
               ? route_estimator(design)
               : rudy(design, {p.bins, p.bins, caph, capv}).overflow;
  };
  Costs scale = cheap(incumbent);
  scale.route = route(incumbent);
  if (!p.normalize_at_seed)
    scale = {1, 1, 1, 1, 1, 1};
  else {
    scale.wire = std::max(1., scale.wire);
    scale.displace =
        std::max(1., (layout.xh - layout.xl) + (layout.yh - layout.yl));
    scale.order = std::max(1., scale.order);
    scale.anchor = std::max(1., scale.anchor);
    scale.slack = 1;
    scale.route = std::max(1., scale.route);
  }
  auto cheap_cost = [&](const Costs &c) {
    return c.wire / scale.wire + p.mu * c.displace / scale.displace +
           p.rho * c.order / scale.order + p.kappa * c.anchor / scale.anchor +
           p.M * c.slack;
  };
  const double basecheap = cheap_cost(cheap(incumbent)),
               baseroute = p.alpha * route(incumbent) / scale.route;
  out.cost = basecheap + baseroute;
  out.seed_cost = out.cost;
  out.rudy_overflow = route(incumbent);
  // Frozen groups must keep their relation to every other group relative to the
  // incumbent sequence pair. Rank each permutation once: the incumbent never
  // changes, so only the proposal needs ranking per evaluation.
  const auto initial = sequence_pair(incumbent.topology);
  std::map<std::string, int> initial_plus_rank, initial_minus_rank;
  for (int k = 0; k < int(initial.plus.size()); ++k)
    initial_plus_rank[initial.plus[k]] = k;
  for (int k = 0; k < int(initial.minus.size()); ++k)
    initial_minus_rank[initial.minus[k]] = k;
  std::vector<std::string> frozen_ids;
  for (const auto &g : incumbent.topology.groups)
    if (g.freeze_topology)
      frozen_ids.push_back(g.id);
  for (int chain = 0; chain < (p.search == "sa" ? p.chains : 1); ++chain) {
    std::mt19937 rng(unsigned(p.seed) + 0x9e3779b9U * unsigned(chain + 1));
    std::uniform_real_distribution<double> unit(0, 1);
    Topology current = incumbent.topology;
    SequencePair currentS = sequence_pair(current);
    double currentcheap = basecheap, currentroute = baseroute;
    auto evaluate = [&](Topology candidate, SequencePair proposed,
                        double temperature, bool greedy) {
      ++out.evaluations;
      // A swap of two other groups can cross a frozen group. Preserve its
      // pairwise relations in both permutations, not just its array index.
      if (!frozen_ids.empty()) {
        const auto rank = [](const std::vector<std::string> &order) {
          std::map<std::string, int> out;
          for (int k = 0; k < int(order.size()); ++k)
            out[order[k]] = k;
          return out;
        };
        const auto plus_rank = rank(proposed.plus),
                   minus_rank = rank(proposed.minus);
        auto preserves = [&](const std::vector<std::string> &before,
                             const std::map<std::string, int> &before_rank,
                             const std::map<std::string, int> &after_rank,
                             const std::string &id) {
          const auto a = before_rank.find(id), b = after_rank.find(id);
          if (a == before_rank.end() || b == after_rank.end())
            return false;
          for (const auto &other : before)
            if (other != id) {
              const auto x = before_rank.find(other),
                         y = after_rank.find(other);
              if (y == after_rank.end() ||
                  ((x->second < a->second) != (y->second < b->second)))
                return false;
            }
          return true;
        };
        for (const auto &id : frozen_ids)
          if (!preserves(initial.plus, initial_plus_rank, plus_rank, id) ||
              !preserves(initial.minus, initial_minus_rank, minus_rank, id))
            return false;
      }
      Realization realized;
      try {
        auto a = abstract_design(src, info, candidate, p);
        auto graph = sequence_graph(candidate, proposed);
        // O(G^2) longest-path packing screens the sequence pair before an LP.
        // This is a necessary condition only; the LP handles obstacles/anchors.
        for (Axis axis : {Axis::x, Axis::y}) {
          const bool x = axis == Axis::x;
          const int n = candidate.groups.size();
          std::map<std::string, int> ids;
          std::vector<std::vector<int>> edges(n);
          std::vector<int> degree(n);
          std::vector<double> pos(n), upper(n), size(n), halo(n);
          for (int i = 0; i < n; ++i) {
            const auto &g = candidate.groups[i];
            ids[g.id] = i;
            size[i] = x ? g.width : g.height;
            halo[i] = x ? g.halo_x : g.halo_y;
            pos[i] = (x ? layout.xl : layout.yl) + halo[i];
            upper[i] = (x ? layout.xh : layout.yh) - size[i] - halo[i];
            if (g.fixed)
              pos[i] = upper[i] = x ? g.x : g.y;
          }
          for (const auto &e : graph)
            if (e.axis == axis) {
              int u = ids.at(e.from), v = ids.at(e.to);
              edges[u].push_back(v);
              ++degree[v];
            }
          std::vector<int> ready;
          for (int i = 0; i < n; ++i)
            if (!degree[i])
              ready.push_back(i);
          int visited = 0;
          while (!ready.empty()) {
            int u = ready.back();
            ready.pop_back();
            ++visited;
            if (pos[u] > upper[u] + 1e-6) {
              ++out.screened;
              return false;
            }
            for (int v : edges[u]) {
              pos[v] = std::max(pos[v], pos[u] + size[u] + halo[u] + halo[v]);
              if (!--degree[v])
                ready.push_back(v);
            }
          }
          if (visited != n) {
            ++out.screened;
            return false;
          }
        }
        realized = compact(a, candidate, p, std::move(graph));
        out.lp_solves += realized.lp_solves;
      } catch (const Infeasible &) {
        return false;
      } catch (const std::invalid_argument &) {
        return false;
      }
      Costs c = cheap(realized);
      double cc = cheap_cost(c), delta = cc - currentcheap;
      // Delayed acceptance: first the cheap target, then the conditional RUDY
      // target. Never choose a best iterate using a partially evaluated cost.
      if (!greedy && delta > 0 && unit(rng) >= std::exp(-delta / temperature))
        return false;
      double rr = p.alpha * route(realized) / scale.route, full = cc + rr;
      bool accept =
          greedy ? full < currentcheap + currentroute
                 : rr <= currentroute ||
                       unit(rng) < std::exp((currentroute - rr) / temperature);
      if (accept) {
        current = realized.topology;
        currentS = proposed;
        currentcheap = cc;
        currentroute = rr;
        ++out.accepted;
      }
      if (realized.legal && full < out.cost) {
        realized.topology.realized = proposed;
        out.best = std::move(realized);
        out.sequence = std::move(proposed);
        out.cost = full;
        out.rudy_overflow =
            p.alpha ? rr * scale.route / p.alpha : route(out.best);
      }
      return accept;
    };
    if (!previous.plus.empty() && previous.plus.size() == current.groups.size())
      evaluate(current, previous, p.T0, true);
    const int n = current.groups.size();
    for (int eval = 0; eval < p.evals; ++eval) {
      Topology candidate = current;
      SequencePair proposal = currentS;
      int i = rng() % n, j = n > 1 ? rng() % n : i;
      if (n > 1 && i == j)
        j = (j + 1) % n;
      int move;
      if (p.search == "lite")
        move = eval % 4;
      else {
        double draw = unit(rng) *
                      (p.move_plus + p.move_minus + p.move_both + p.move_flip);
        move = draw < p.move_plus                                ? 0
               : draw < p.move_plus + p.move_minus               ? 1
               : draw < p.move_plus + p.move_minus + p.move_both ? 2
                                                                 : 3;
      }
      if (move == 3) {
        if (!p.orient_moves)
          continue;
        auto &g = candidate.groups[i];
        if (g.fixed || g.freeze_topology)
          continue;
        std::vector<Orientation> allowed;
        for (unsigned o = 0; o < 4; ++o)
          if ((g.allowed_orientations & (1U << o)) &&
              o != unsigned(g.orientation))
            allowed.push_back(static_cast<Orientation>(o));
        if (allowed.empty())
          continue;
        g.orientation = allowed[rng() % allowed.size()];
      } else {
        if (n < 2)
          continue;
        const auto &ga = candidate.groups[i];
        const auto &gb = candidate.groups[j];
        if (ga.fixed || gb.fixed || ga.freeze_topology || gb.freeze_topology)
          continue;
        auto swap = [&](auto &v) {
          auto a = std::find(v.begin(), v.end(), ga.id),
               b = std::find(v.begin(), v.end(), gb.id);
          std::iter_swap(a, b);
        };
        if (move == 0 || move == 2)
          swap(proposal.plus);
        if (move == 1 || move == 2)
          swap(proposal.minus);
      }
      evaluate(std::move(candidate), std::move(proposal),
               std::max(1e-12, p.T0 * std::pow(p.cooling, eval)),
               p.search == "lite");
    }
  }

  return out;
}
} // namespace DreamPlace::native::macro
