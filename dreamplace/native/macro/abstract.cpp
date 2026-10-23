// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/abstract.hpp"
#include "dreamplace/native/rigid.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <tuple>
namespace DreamPlace::native::macro {
Abstracted abstract_design(const Design &src, std::span<const NodeInfo> info,
                           const Topology &t, const Parameters &p,
                           bool freeze) {
  validate(t, src, info);
  p.validate();
  Abstracted a;
  auto &d = a.design;
  d.layout = src.layout;
  struct Item {
    int source = -1, group = -1, blockage = -1;
    bool movable = false;
  };
  std::vector<Item> items;
  std::map<std::string, int> names;
  for (int i = 0; i < int(info.size()); ++i)
    names[info[i].name] = i;
  for (const auto &n : info)
    a.source_names.push_back(n.name);
  a.source_to_node.assign(info.size(), -1);
  std::set<int> grouped;
  for (int g = 0; g < int(t.groups.size()); ++g) {
    const auto &group = t.groups[g];
    bool movable = !freeze && !group.fixed;
    for (const auto &m : group.members) {
      int i = names.at(m.name);
      grouped.insert(i);
      movable &= i < src.num_movable;
    }
    items.push_back({-1, g, -1, movable});
  }
  for (int i = 0; i < int(info.size()); ++i)
    if (!grouped.contains(i))
      items.push_back({i, -1, -1, i < src.num_movable});
  for (int b = 0; b < int(t.blockages.size()); ++b)
    items.push_back({-1, -1, b, false});
  std::stable_sort(
      items.begin(), items.end(),
      [](const Item &a, const Item &b) { return a.movable > b.movable; });
  for (const auto &item : items) {
    int index = d.node_x.size();
    double x, y, w, h, hx = 0, hy = 0;
    Orientation o = Orientation::N;
    NodeInfo ni;
    PositionBox box{};
    int source = item.source;
    if (item.group >= 0) {
      const auto &g = t.groups[item.group];
      a.group_node[g.id] = index;
      x = g.x;
      y = g.y;
      w = g.width;
      h = g.height;
      hx = g.halo_x;
      hy = g.halo_y;
      o = g.orientation;
      ni = {"@group:" + g.id,      {}, g.id, true, false, false, true,
            g.allowed_orientations};
      for (const auto &m : g.members)
        a.source_to_node[names.at(m.name)] = index;
    } else if (item.blockage >= 0) {
      const auto &b = t.blockages[item.blockage];
      x = b.rect.xl;
      y = b.rect.yl;
      w = b.rect.xh - x;
      h = b.rect.yh - y;
      ni = {"@blockage:" + b.id, {}, {}, true};
    } else {
      x = src.node_x[source];
      y = src.node_y[source];
      w = src.node_size_x[source];
      h = src.node_size_y[source];
      hx = src.node_halo_x.empty() ? 0 : src.node_halo_x[source];
      hy = src.node_halo_y.empty() ? 0 : src.node_halo_y[source];
      o = src.node_orientation.empty() ? Orientation::N
                                       : src.node_orientation[source];
      ni = info[source];
      a.source_to_node[source] = index;
    }
    box = {d.layout.xl, d.layout.yl, d.layout.xh - w, d.layout.yh - h};
    if (source >= 0 && !src.node_boxes.empty())
      box = src.node_boxes[source];
    if (item.group >= 0 && !freeze && p.anchor_mode == "exact")
      for (const auto &anchor : t.anchors)
        if (anchor.group == t.groups[item.group].id &&
            (anchor.strength.hard || anchor.strength.weight > 0)) {
          auto candidate = box;
          const auto &l = d.layout;
          switch (anchor.side) {
          case Side::left:
            candidate.xl = candidate.xh = l.xl + hx;
            break;
          case Side::right:
            candidate.xl = candidate.xh = l.xh - w - hx;
            break;
          case Side::bottom:
            candidate.yl = candidate.yh = l.yl + hy;
            break;
          case Side::top:
            candidate.yl = candidate.yh = l.yh - h - hy;
            break;
          }
          // Keep all compatible exact projections. In a tight/conflicting case,
          // Compact's explicit slack reports the deviation instead of an empty
          // clamp.
          candidate.xl = std::max(candidate.xl, box.xl);
          candidate.xh = std::min(candidate.xh, box.xh);
          candidate.yl = std::max(candidate.yl, box.yl);
          candidate.yh = std::min(candidate.yh, box.yh);
          if (candidate.xl <= candidate.xh && candidate.yl <= candidate.yh)
            box = candidate;
          else if (anchor.strength.hard)
            throw std::invalid_argument("conflicting hard anchors: " +
                                        anchor.id);
        }
    if (item.movable)
      ++d.num_movable;
    d.node_x.push_back(x);
    d.node_y.push_back(y);
    d.node_size_x.push_back(w);
    d.node_size_y.push_back(h);
    d.node_halo_x.push_back(hx);
    d.node_halo_y.push_back(hy);
    d.node_orientation.push_back(o);
    d.node_boxes.push_back(box);
    d.node_inflatable.push_back(item.movable && !ni.macro);
    a.info.push_back(std::move(ni));
    a.source_of_node.push_back(source);
  }
  d.pin2node = src.pin2node;
  d.pin_offset_x = src.pin_offset_x;
  d.pin_offset_y = src.pin_offset_y;
  d.flat_net2pin = src.flat_net2pin;
  d.flat_net2pin_start = src.flat_net2pin_start;
  d.net_weights = src.net_weights;
  d.net_mask = src.net_mask;
  struct Membership {
    const Group *group;
    const Member *member;
  };
  std::map<int, Membership> membership;
  for (const auto &g : t.groups)
    for (const auto &m : g.members)
      membership[names.at(m.name)] = {&g, &m};
  std::vector<RigidPinFrame> frames(info.size());
  for (int node = 0; node < int(info.size()); ++node) {
    auto &frame = frames[node];
    frame.node = a.source_to_node[node];
    auto it = membership.find(node);
    if (it != membership.end()) {
      const auto &[g, m] = it->second;
      auto xy =
          member_position(*g, *m, src.node_size_x[node], src.node_size_y[node]);
      Orientation old = src.node_orientation.empty()
                            ? Orientation::N
                            : src.node_orientation[node];
      frame.x = xy.first - g->x;
      frame.y = xy.second - g->y;
      frame.transform = compose(old, compose(g->orientation, m->orientation));
    }
  }
  rebind_rigid_pins(src, d, frames);
  for (int net = 0; net < int(d.net_weights.size()); ++net) {
    int first = -1;
    bool one = true;
    for (int k = d.flat_net2pin_start[net]; k < d.flat_net2pin_start[net + 1];
         ++k) {
      int n = d.pin2node[d.flat_net2pin[k]];
      if (first < 0)
        first = n;
      else
        one &= n == first;
    }
    if (one)
      d.net_mask[net] = 0;
  }
  for (const auto &e : t.order) {
    auto u = a.group_node.at(e.from), v = a.group_node.at(e.to);
    d.order_edges.push_back({u, v, e.axis, 0, e.strength, e.id});
  }
  d.build_maps();
  d.validate();
  return a;
}
void update_groups(const Abstracted &a, Topology &t) {
  for (auto &g : t.groups) {
    int i = a.group_node.at(g.id);
    g.x = a.design.node_x[i];
    g.y = a.design.node_y[i];
    if (!a.design.node_orientation.empty())
      g.orientation = a.design.node_orientation[i];
  }
}
void expand(const Abstracted &a, const Topology &topology, Design &src) {
  auto t = topology;
  update_groups(a, t);
  Design out = src;
  if (out.node_orientation.empty())
    out.node_orientation.assign(out.node_x.size(), Orientation::N);
  for (int i = 0; i < int(a.source_of_node.size()); ++i) {
    int source = a.source_of_node[i];
    if (source >= 0) {
      out.node_x[source] = a.design.node_x[i];
      out.node_y[source] = a.design.node_y[i];
    }
  }
  std::map<std::string, int> source_ids;
  for (int i = 0; i < int(a.source_names.size()); ++i)
    source_ids[a.source_names[i]] = i;
  for (const auto &g : t.groups)
    for (const auto &m : g.members) {
      auto found = source_ids.find(m.name);
      int source = found == source_ids.end() ? -1 : found->second;
      if (source < 0)
        throw std::invalid_argument("missing expanded member");
      auto xy = member_position(g, m, src.node_size_x[source],
                                src.node_size_y[source]);
      out.node_x[source] = xy.first;
      out.node_y[source] = xy.second;
      const auto orientation = compose(g.orientation, m.orientation);
      const auto old = src.node_orientation.empty()
                           ? Orientation::N
                           : src.node_orientation[source];
      for (int k = src.flat_node2pin_start[source];
           k < src.flat_node2pin_start[source + 1]; ++k) {
        int pin = src.flat_node2pin[k];
        auto offset =
            orient_pin(src.pin_offset_x[pin], src.pin_offset_y[pin],
                       src.node_size_x[source], src.node_size_y[source],
                       compose(old, orientation));
        out.pin_offset_x[pin] = offset.first;
        out.pin_offset_y[pin] = offset.second;
      }
      out.node_orientation[source] = orientation;
    }
  for (int i = src.num_movable; i < int(src.node_x.size()); ++i)
    if (out.node_x[i] != src.node_x[i] || out.node_y[i] != src.node_y[i] ||
        out.node_orientation[i] != (src.node_orientation.empty()
                                        ? Orientation::N
                                        : src.node_orientation[i]))
      throw std::runtime_error("H2D moved a fixed input");
  out.validate();
  src = std::move(out);
}
void add_pseudonets(Abstracted &a, const Topology &t, const Parameters &p) {
  if (!p.pseudonets || !p.beta || t.dataflow.empty())
    return;
  auto &d = a.design;
  std::map<std::string, int> endpoint;
  std::map<std::string, std::pair<std::pair<double, double>, int>> clusters;
  for (const auto &g : t.groups)
    for (const auto &m : g.members)
      endpoint[m.name] = a.group_node.at(g.id);
  for (int i = 0; i < int(a.info.size()); ++i) {
    const auto &n = a.info[i];
    if (n.io)
      endpoint[n.name] = i;
    else if (!n.macro) {
      auto &c = clusters[n.parent];
      c.first.first += d.node_x[i] + d.node_size_x[i] / 2;
      c.first.second += d.node_y[i] + d.node_size_y[i] / 2;
      ++c.second;
    }
  }
  std::vector<double> weights;
  for (const auto &f : t.dataflow)
    weights.push_back(f.weight);
  std::sort(weights.begin(), weights.end());
  double threshold = weights[std::min(
      weights.size() - 1, size_t(p.theta_q * (weights.size() - 1)))];
  auto resolve = [&](const std::string &name) -> int {
    auto it = endpoint.find(name);
    if (it != endpoint.end())
      return it->second;
    std::string parent;
    if (name.starts_with("@cluster:"))
      parent = name.substr(9);
    else {
      auto original =
          std::find(a.source_names.begin(), a.source_names.end(), name);
      if (original == a.source_names.end())
        return -1;
      const int source = original - a.source_names.begin();
      parent = a.info[a.source_to_node[source]].parent;
    }
    auto c = clusters.find(parent);
    if (c == clusters.end())
      return -1;
    std::string key = "@cluster:" + parent;
    auto found = endpoint.find(key);
    int id;
    if (found != endpoint.end())
      id = found->second;
    else {
      id = d.node_x.size();
      endpoint[key] = id;
      d.node_x.push_back(c->second.first.first / c->second.second);
      d.node_y.push_back(c->second.first.second / c->second.second);
      d.node_size_x.push_back(0);
      d.node_size_y.push_back(0);
      d.node_halo_x.push_back(0);
      d.node_halo_y.push_back(0);
      d.node_orientation.push_back(Orientation::N);
      d.node_boxes.push_back(
          {d.layout.xl, d.layout.yl, d.layout.xh, d.layout.yh});
      d.node_inflatable.push_back(0);
      a.info.push_back({key, parent, {}, false, false, true});
      a.source_of_node.push_back(-1);
    }
    return id;
  };
  for (const auto &f : t.dataflow)
    if (f.weight >= threshold) {
      int u = resolve(f.from), v = resolve(f.to);
      if (u < 0 || v < 0 || u == v)
        continue;
      for (int node : {u, v}) {
        d.flat_net2pin.push_back(d.pin2node.size());
        d.pin2node.push_back(node);
        d.pin_offset_x.push_back(d.node_size_x[node] / 2);
        d.pin_offset_y.push_back(d.node_size_y[node] / 2);
      }
      d.flat_net2pin_start.push_back(d.flat_net2pin.size());
      d.net_weights.push_back(p.beta * f.weight);
      d.net_mask.push_back(1);
    }
  d.build_maps();
  d.validate();
}
} // namespace DreamPlace::native::macro
