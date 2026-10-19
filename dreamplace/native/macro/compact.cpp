// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/compact.hpp"
#include "dreamplace/native/macro/check.hpp"
#include "dreamplace/native/macro/packing.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_set>
namespace DreamPlace::native::macro {
namespace {
struct Object {
  std::string id;
  int node;
  bool fixed;
};
std::vector<Object> objects(const Abstracted &a, const Topology &t) {
  std::vector<Object> out;
  std::set<int> used;
  for (const auto &g : t.groups) {
    int i = a.group_node.at(g.id);
    used.insert(i);
    out.push_back({g.id, i, i >= a.design.num_movable});
  }
  for (int i = a.design.num_movable; i < int(a.design.node_x.size()); ++i)
    if (!used.contains(i) && a.design.node_size_x[i] > 0 &&
        a.design.node_size_y[i] > 0)
      out.push_back({a.info[i].name, i, true});
  return out;
}
double hx(const Design &d, int i) {
  return d.node_halo_x.empty()
             ? 0
             : std::ceil(d.node_halo_x[i] / d.layout.site_width - 1e-9) *
                   d.layout.site_width;
}
double hy(const Design &d, int i) {
  return d.node_halo_y.empty()
             ? 0
             : std::ceil(d.node_halo_y[i] / d.layout.row_height - 1e-9) *
                   d.layout.row_height;
}
Separation relation(const Abstracted &a, const Object &u, const Object &v) {
  const auto &d = a.design;
  int i = u.node, j = v.node;
  double cx = d.node_x[j] + d.node_size_x[j] / 2 - d.node_x[i] -
              d.node_size_x[i] / 2,
         cy = d.node_y[j] + d.node_size_y[j] / 2 - d.node_y[i] -
              d.node_size_y[i] / 2;
  double wx = (d.node_size_x[i] + d.node_size_x[j]) / 2 + hx(d, i) + hx(d, j),
         wy = (d.node_size_y[i] + d.node_size_y[j]) / 2 + hy(d, i) + hy(d, j);
  bool horizontal =
      std::abs(cx) - wx >= 0
          ? (std::abs(cy) - wy < 0 || std::abs(cx) / wx >= std::abs(cy) / wy)
          : (std::abs(cy) - wy < 0 && wx - std::abs(cx) <= wy - std::abs(cy));
  bool forward = horizontal ? cx > 0 : cy > 0;
  if ((horizontal ? cx : cy) == 0)
    forward = u.id < v.id;
  return {forward ? u.id : v.id, forward ? v.id : u.id,
          "separate:" + u.id + ":" + v.id, horizontal ? Axis::x : Axis::y};
}
long long coordinate(double v) {
  if (!std::isfinite(v) || std::abs(v) > 1e9)
    throw std::invalid_argument("H2D pitch coordinate out of range");
  return std::llround(v);
}
} // namespace
SequencePair sequence_pair(const Topology &t) {
  SequencePair s;
  std::map<std::string, std::pair<double, double>> c;
  for (const auto &g : t.groups) {
    s.plus.push_back(g.id);
    s.minus.push_back(g.id);
    c[g.id] = {g.x + g.width / 2, g.y + g.height / 2};
  }
  auto sort = [&](auto &v, bool plus) {
    std::sort(v.begin(), v.end(), [&](const auto &a, const auto &b) {
      double x = c[a].first + (plus ? -1 : 1) * c[a].second,
             y = c[b].first + (plus ? -1 : 1) * c[b].second;
      return x < y || (x == y && a < b);
    });
  };
  sort(s.plus, true);
  sort(s.minus, false);
  return s;
}
std::vector<Separation> sequence_graph(const Topology &t,
                                       const SequencePair &s) {
  if (s.plus.size() != t.groups.size() || s.minus.size() != t.groups.size())
    throw std::invalid_argument("sequence pair size mismatch");
  std::set<std::string> names;
  for (const auto &g : t.groups)
    names.insert(g.id);
  if (std::set<std::string>(s.plus.begin(), s.plus.end()) != names ||
      std::set<std::string>(s.minus.begin(), s.minus.end()) != names)
    throw std::invalid_argument("sequence pair is not a group permutation");
  std::map<std::string, int> minus;
  for (int i = 0; i < int(s.minus.size()); ++i)
    minus[s.minus[i]] = i;
  std::vector<Separation> out;
  for (int i = 0; i < int(s.plus.size()); ++i)
    for (int j = i + 1; j < int(s.plus.size()); ++j) {
      const auto &u = s.plus[i];
      const auto &v = s.plus[j];
      bool horizontal = minus[u] < minus[v];
      out.push_back({horizontal ? u : v, horizontal ? v : u,
                     "sp:" + u + ":" + v, horizontal ? Axis::x : Axis::y});
    }
  return out;
}
std::vector<Separation> position_graph(const Abstracted &a, const Topology &t) {
  auto obj = objects(a, t);
  std::vector<Separation> out;
  for (size_t i = 0; i < obj.size(); ++i)
    for (size_t j = i + 1; j < obj.size(); ++j)
      if (!obj[i].fixed || !obj[j].fixed)
        out.push_back(relation(a, obj[i], obj[j]));
  return out;
}
double anchor_deviation(const Topology &t, const Layout &l) {
  std::map<std::string, const Group *> groups;
  for (const auto &g : t.groups)
    groups[g.id] = &g;
  double total = 0;
  for (const auto &a : t.anchors) {
    const auto &g = *groups.at(a.group);
    double distance = 0;
    switch (a.side) {
    case Side::left:
      distance = std::abs(g.x - g.halo_x - l.xl);
      break;
    case Side::right:
      distance = std::abs(l.xh - g.x - g.width - g.halo_x);
      break;
    case Side::bottom:
      distance = std::abs(g.y - g.halo_y - l.yl);
      break;
    case Side::top:
      distance = std::abs(l.yh - g.y - g.height - g.halo_y);
      break;
    }
    total += a.strength.weight * std::max(0., distance - a.tolerance);
  }
  return total;
}
double order_violation(const Topology &t) {
  std::map<std::string, const Group *> groups;
  for (const auto &g : t.groups)
    groups[g.id] = &g;
  double total = 0;
  for (const auto &e : t.order) {
    const auto &a = *groups.at(e.from);
    const auto &b = *groups.at(e.to);
    total += e.strength.weight * std::max(0., e.axis == Axis::x
                                                  ? a.x + a.width - b.x
                                                  : a.y + a.height - b.y);
  }
  return total;
}
Realization compact(const Abstracted &a, const Topology &t, const Parameters &p,
                    std::vector<Separation> graph, bool displacement_only) {
  p.validate();
  a.design.validate();
  Realization out;
  out.topology = t;
  const auto &d = a.design;
  const auto obj = objects(a, t);
  std::map<std::string, int> index;
  std::map<int, int> physical;
  for (int i = 0; i < int(obj.size()); ++i) {
    index[obj[i].id] = i;
    physical[obj[i].node] = i;
  }
  // Edges are matched by id or by unordered endpoint pair. Hash those keys once
  // instead of rescanning the whole (quadratic-size) graph per candidate edge.
  using PairKey = std::pair<std::string, std::string>;
  const auto pair_key = [](const std::string &a, const std::string &b) {
    return a < b ? PairKey{a, b} : PairKey{b, a};
  };
  // Retain obstacle separation independently of a movable sequence-pair
  // proposal.
  auto obstacle = position_graph(a, t);
  std::unordered_set<std::string> present;
  for (const auto &g : graph)
    present.insert(g.id);
  for (const auto &e : obstacle) {
    int u = index.at(e.from), v = index.at(e.to);
    if (obj[u].fixed || obj[v].fixed)
      if (present.insert(e.id).second)
        graph.push_back(e);
  }
  std::set<std::string> frozen;
  for (const auto &g : t.groups)
    if (g.freeze_topology)
      frozen.insert(g.id);
  if (!frozen.empty()) {
    // Each pinned realized edge replaces every graph edge on its pair, and a
    // later pinned edge on the same pair replaces an earlier one.
    std::vector<const Separation *> pinned;
    std::map<PairKey, size_t> last;
    for (const auto &e : t.realized_graph)
      if ((frozen.contains(e.from) || frozen.contains(e.to)) &&
          index.contains(e.from) && index.contains(e.to)) {
        last[pair_key(e.from, e.to)] = pinned.size();
        pinned.push_back(&e);
      }
    if (!pinned.empty()) {
      std::erase_if(graph, [&](const Separation &v) {
        return last.contains(pair_key(v.from, v.to));
      });
      for (size_t k = 0; k < pinned.size(); ++k)
        if (last.at(pair_key(pinned[k]->from, pinned[k]->to)) == k)
          graph.push_back(*pinned[k]);
    }
  }
  std::vector<OrderEdge> forced;
  for (const auto &e : t.order)
    if (p.order_in_lp || e.strength.hard)
      forced.push_back(e);
  // A requested order owns its pair's separating direction. The geometric
  // prototype must not inject the opposite relation for that same pair.
  {
    std::set<PairKey> owned;
    for (const auto &e : forced)
      owned.insert(pair_key(e.from, e.to));
    if (!owned.empty())
      std::erase_if(graph, [&](const Separation &s) {
        return owned.contains(pair_key(s.from, s.to));
      });
  }
  std::vector<PositionBox> positions(obj.size());
  std::set<std::string> relaxed;
  for (int axis = 0; axis < 2; ++axis) {
    const double pitch =
        axis ? (p.pitch_y > 0 ? p.pitch_y : d.layout.row_height)
             : (p.pitch_x > 0 ? p.pitch_x : d.layout.site_width);
    const double site = axis ? d.layout.row_height : d.layout.site_width;
    if (std::abs(pitch / site - std::round(pitch / site)) > 1e-6)
      throw std::invalid_argument(
          "H2D snap pitch must be a whole site/row multiple");
    const double origin = axis ? d.layout.yl : d.layout.xl,
                 extent = (axis ? d.layout.yh : d.layout.xh) - origin;
    double span = extent / pitch + 1;
    for (const auto &o : obj)
      span +=
          (axis ? d.node_size_y[o.node] : d.node_size_x[o.node]) / pitch + 1;
    long long limit = coordinate(std::ceil(span) * 2 + 100);
    while (true) {
      DifferenceLP lp;
      std::vector<int> vars;
      std::vector<double> residual;
      for (const auto &o : obj) {
        int i = o.node;
        double pos = axis ? d.node_y[i] : d.node_x[i],
               size = axis ? d.node_size_y[i] : d.node_size_x[i],
               halo = axis ? hy(d, i) : hx(d, i);
        long long here =
            coordinate(o.fixed ? std::floor((pos - origin) / pitch)
                               : std::round((pos - origin) / pitch));
        const double offset = o.fixed ? pos - origin - here * pitch : 0;
        residual.push_back(offset);
        int v = lp.variable(o.fixed ? here : -limit, o.fixed ? here : limit, 0,
                            o.id);
        vars.push_back(v);
        lp.constrain(-1, v,
                     coordinate(std::floor(
                         (extent - size - halo - offset) / pitch + 1e-9)),
                     p.M, "outline:" + o.id + ":upper");
        lp.constrain(v, -1,
                     -coordinate(std::ceil((halo - offset) / pitch - 1e-9)),
                     p.M, "outline:" + o.id + ":lower");
        if (displacement_only && !o.fixed)
          lp.anchor(v, here, 1);
      }
      for (const auto &e : graph)
        if ((e.axis == Axis::y) == bool(axis)) {
          int u = index.at(e.from), v = index.at(e.to), i = obj[u].node,
              j = obj[v].node;
          double gap = (axis ? d.node_size_y[i] + hy(d, i) + hy(d, j)
                             : d.node_size_x[i] + hx(d, i) + hx(d, j));
          lp.constrain(vars[v], vars[u],
                       -coordinate(std::ceil(
                           (gap + residual[u] - residual[v]) / pitch - 1e-9)),
                       displacement_only ? 0 : p.M, e.id);
        }
      for (const auto &e : forced)
        if (!relaxed.contains(e.id) && (e.axis == Axis::y) == bool(axis)) {
          int u = index.at(e.from), v = index.at(e.to), i = obj[u].node;
          int j = obj[v].node;
          double size = axis ? d.node_size_y[i] + hy(d, i) + hy(d, j)
                             : d.node_size_x[i] + hx(d, i) + hx(d, j);
          lp.constrain(vars[v], vars[u],
                       -coordinate(std::ceil(
                           (size + residual[u] - residual[v]) / pitch - 1e-9)),
                       0, e.id);
        }
      for (const auto &anchor : t.anchors) {
        if (!anchor.strength.hard && anchor.strength.weight == 0)
          continue;
        bool horizontal =
            anchor.side == Side::left || anchor.side == Side::right;
        if (horizontal == bool(axis))
          continue;
        int i = index.at(anchor.group), node = obj[i].node;
        double size = axis ? d.node_size_y[node] : d.node_size_x[node],
               halo = axis ? hy(d, node) : hx(d, node);
        bool lower = anchor.side == Side::left || anchor.side == Side::bottom;
        double target = lower ? halo : extent - size - halo;
        long long at = coordinate(std::round((target - residual[i]) / pitch));
        double penalty =
            anchor.strength.hard
                ? 0
                : (p.anchor_mode == "exact" ? p.M * anchor.strength.weight
                                            : anchor.strength.weight);
        lp.constrain(-1, vars[i], at, penalty, anchor.id + ":upper");
        lp.constrain(vars[i], -1, -at, penalty, anchor.id + ":lower");
      }
      if (!displacement_only)
        for (int net = 0; net < int(d.net_weights.size()); ++net) {
          if (!d.net_mask[net] || !d.net_weights[net])
            continue;
          std::map<int, std::pair<double, double>> pins;
          double fixedmin = 1e100, fixedmax = -1e100;
          for (int k = d.flat_net2pin_start[net];
               k < d.flat_net2pin_start[net + 1]; ++k) {
            int pin = d.flat_net2pin[k], node = d.pin2node[pin];
            double off = axis ? d.pin_offset_y[pin] : d.pin_offset_x[pin];
            auto it = physical.find(node);
            if (it != physical.end() && !obj[it->second].fixed) {
              int v = it->second;
              auto f = pins.find(v);
              if (f == pins.end())
                pins[v] = {off, off};
              else {
                f->second.first = std::min(f->second.first, off);
                f->second.second = std::max(f->second.second, off);
              }
            } else {
              double pos =
                  (axis ? d.node_y[node] : d.node_x[node]) + off - origin;
              fixedmin = std::min(fixedmin, pos);
              fixedmax = std::max(fixedmax, pos);
            }
          }
          double weight = d.net_weights[net] * pitch;
          if (pins.size() == 1 && fixedmax >= fixedmin) {
            auto [v, off] = *pins.begin();
            lp.constrain(
                -1, vars[v],
                coordinate(std::floor((fixedmax - off.second) / pitch)), weight,
                "net:" + std::to_string(net) + ":max");
            lp.constrain(vars[v], -1,
                         coordinate(std::floor((off.first - fixedmin) / pitch)),
                         weight, "net:" + std::to_string(net) + ":min");
          } else if (pins.size() >= 2) {
            int lo = lp.variable(-limit, limit, -weight, "net.lo"),
                hi = lp.variable(-limit, limit, weight, "net.hi");
            for (auto [v, off] : pins) {
              lp.constrain(vars[v], lo,
                           coordinate(std::floor(off.first / pitch)), 0,
                           "pin.lo");
              lp.constrain(hi, vars[v],
                           -coordinate(std::ceil(off.second / pitch)), 0,
                           "pin.hi");
            }
            if (fixedmax >= fixedmin) {
              lp.constrain(-1, lo, coordinate(std::floor(fixedmin / pitch)), 0,
                           "fixed.pin.lo");
              lp.constrain(hi, -1, -coordinate(std::ceil(fixedmax / pitch)), 0,
                           "fixed.pin.hi");
            }
          }
        }
      try {
        auto solved = lp.solve();
        ++out.lp_solves;
        for (int i = 0; i < int(obj.size()); ++i) {
          double value = obj[i].fixed ? (axis ? d.node_y[obj[i].node]
                                              : d.node_x[obj[i].node])
                                      : origin + pitch * solved.x[vars[i]];
          if (axis)
            positions[i].yl = value;
          else
            positions[i].xl = value;
        }
        for (const auto &s :
             solved.slack) { // HPWL hinge arcs are objective terms, not
                             // feasibility violations.
          if (s.id.starts_with("net:"))
            continue;
          auto physical_slack = s;
          physical_slack.amount *= coordinate(pitch);
          out.slack_amount += s.amount * pitch;
          out.slack.push_back(std::move(physical_slack));
        }
        break;
      } catch (const Infeasible &e) {
        const OrderEdge *drop = nullptr;
        for (const auto &r : forced)
          if (!r.strength.hard && !relaxed.contains(r.id) &&
              std::find(e.cycle.begin(), e.cycle.end(), r.id) !=
                  e.cycle.end() &&
              (drop == nullptr || r.strength.weight < drop->strength.weight ||
               (r.strength.weight == drop->strength.weight && r.id < drop->id)))
            drop = &r;
        if (!drop)
          throw;
        relaxed.insert(drop->id);
        out.relaxed.push_back(drop->id);
      }
    }
  }
  auto placed = a;
  for (int i = 0; i < int(obj.size()); ++i) {
    placed.design.node_x[obj[i].node] = positions[i].xl;
    placed.design.node_y[obj[i].node] = positions[i].yl;
  }
  update_groups(placed, out.topology);
  std::vector<int> nodes;
  for (const auto &o : obj)
    nodes.push_back(o.node);
  out.legal = check(placed.design, nodes).legal();
  out.hpwl = hpwl(placed.design);
  out.anchor_deviation = anchor_deviation(out.topology, d.layout);
  out.order_violation = order_violation(out.topology);
  for (const auto &e : forced)
    if (!relaxed.contains(e.id))
      graph.push_back({e.from, e.to, e.id, e.axis});
  out.graph = std::move(graph);
  out.topology.realized_graph = out.graph;
  return out;
}
Realization group_legalize(const Abstracted &a, const Topology &t,
                           const Parameters &p) {
  auto graph = position_graph(a, t);
  Realization result;
  try {
    result = compact(a, t, p, graph, true);
    if (result.legal)
      return result;
  } catch (const Infeasible &e) {
    if (std::none_of(e.cycle.begin(), e.cycle.end(), [](const auto &id) {
          return id.starts_with("separate:");
        }))
      throw;
  }
  // Reuse the existing macro heuristic/LP/Hannan launcher on group rectangles.
  // Movable standard cells are omitted; fixed objects remain exact obstacles.
  auto obj = objects(a, t);
  std::stable_sort(
      obj.begin(), obj.end(),
      [](const Object &a, const Object &b) { return a.fixed < b.fixed; });
  Design packed;
  packed.layout = a.design.layout;
  for (const auto &o : obj) {
    int i = o.node;
    if (!o.fixed)
      ++packed.num_movable;
    packed.node_x.push_back(a.design.node_x[i]);
    packed.node_y.push_back(a.design.node_y[i]);
    packed.node_size_x.push_back(a.design.node_size_x[i]);
    packed.node_size_y.push_back(a.design.node_size_y[i]);
    packed.node_halo_x.push_back(hx(a.design, i));
    packed.node_halo_y.push_back(hy(a.design, i));
  }
  packed.build_maps();
  Options control;
  control.max_iterations = 0;
  control.detailed_place = false;
  bool packed_fallback = false;
  try {
    place(packed, control);
  } catch (const std::runtime_error &) {
    if (!pack_rectangles(packed))
      throw;
    packed_fallback = true;
  }
  auto fallback = a;
  for (int i = 0; i < int(obj.size()); ++i) {
    fallback.design.node_x[obj[i].node] = packed.node_x[i];
    fallback.design.node_y[obj[i].node] = packed.node_y[i];
  }
  Topology legal = t;
  if (packed_fallback)
    legal.diagnostics.push_back("legacy group legalization failed; used "
                                "deterministic rectangle packing");
  update_groups(fallback, legal);
  auto legalgraph = position_graph(fallback, legal);
  auto refined = compact(fallback, legal, p, legalgraph, true);
  if (refined.legal)
    return refined;
  // The heuristic cannot certify hard guidance. Do not return an apparently
  // legal packing after the constrained refinement rejected it.
  if (std::any_of(t.order.begin(), t.order.end(),
                  [](const auto &e) { return e.strength.hard; }) ||
      std::any_of(t.anchors.begin(), t.anchors.end(),
                  [](const auto &a) { return a.strength.hard; }))
    throw Infeasible({"hard H2D constraints prevent a legal packing"});
  result.topology = std::move(legal);
  result.graph = std::move(legalgraph);
  result.legal = true;
  result.hpwl = hpwl(fallback.design);
  result.anchor_deviation = anchor_deviation(result.topology, a.design.layout);
  result.order_violation = order_violation(result.topology);
  result.slack = std::move(refined.slack);
  result.slack_amount = refined.slack_amount;
  result.lp_solves += refined.lp_solves;
  return result;
}
} // namespace DreamPlace::native::macro
