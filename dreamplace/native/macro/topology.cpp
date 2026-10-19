// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/topology.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
namespace DreamPlace::native::macro {
namespace {
void need(bool ok, std::string message) {
  if (!ok)
    throw std::invalid_argument(std::move(message));
}
std::map<std::string, int> index(std::span<const NodeInfo> info,
                                 const Design &d) {
  need(info.size() == d.node_x.size(), "H2D metadata size mismatch");
  std::map<std::string, int> out;
  for (int i = 0; i < int(info.size()); ++i)
    need(!info[i].name.empty() && out.emplace(info[i].name, i).second,
         "duplicate/empty physical name");
  return out;
}
} // namespace
std::string fingerprint(const Design &d, std::span<const NodeInfo> info,
                        bool full) {
  d.validate();
  const auto ids = index(info, d);
  std::uint64_t h = 14695981039346656037ULL;
  auto hash = [&](std::string_view text) {
    for (unsigned char c : text)
      h = (h ^ c) * 1099511628211ULL;
    h = (h ^ 0xff) * 1099511628211ULL;
  };
  auto number = [&](double v) {
    hash(std::to_string(std::bit_cast<std::uint64_t>(v)));
  };
  for (double v : {d.layout.xl, d.layout.yl, d.layout.xh, d.layout.yh,
                   d.layout.site_width, d.layout.row_height})
    number(v);
  for (const auto &[name, i] : ids)
    if (full || info[i].macro) {
      hash(name);
      hash(info[i].parent);
      hash(info[i].master);
      hash(std::to_string(info[i].allowed_orientations));
      number(d.node_size_x[i]);
      number(d.node_size_y[i]);
      number(d.node_halo_x.empty() ? 0 : d.node_halo_x[i]);
      number(d.node_halo_y.empty() ? 0 : d.node_halo_y[i]);
      hash(i < d.num_movable ? "movable" : "fixed");
      if (i >= d.num_movable) {
        number(d.node_x[i]);
        number(d.node_y[i]);
        hash(std::to_string(
            d.node_orientation.empty() ? 0 : int(d.node_orientation[i])));
      }
      std::vector<std::pair<double, double>> offsets;
      for (int k = d.flat_node2pin_start[i]; k < d.flat_node2pin_start[i + 1];
           ++k) {
        int p = d.flat_node2pin[k];
        auto xy =
            orient_pin(d.pin_offset_x[p], d.pin_offset_y[p], d.node_size_x[i],
                       d.node_size_y[i],
                       d.node_orientation.empty() ? Orientation::N
                                                  : d.node_orientation[i]);
        offsets.push_back(xy);
      }
      std::sort(offsets.begin(), offsets.end());
      for (auto [x, y] : offsets) {
        number(x);
        number(y);
      }
    }
  if (full)
    for (size_t n = 0; n < d.net_weights.size(); ++n) {
      number(d.net_weights[n]);
      hash(std::to_string(d.net_mask[n]));
      for (int k = d.flat_net2pin_start[n]; k < d.flat_net2pin_start[n + 1];
           ++k)
        hash(info[d.pin2node[d.flat_net2pin[k]]].name);
    }
  std::ostringstream s;
  s << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << h;
  return s.str();
}
std::string topology_key(const Topology &t) {
  std::ostringstream s;
  s << std::setprecision(17);
  auto str = [&](const std::string &v) { s << v.size() << ":" << v << ";"; };
  auto strength = [&](Strength v) { s << v.hard << "," << v.weight << ";"; };
  str(t.design);
  for (const auto &g : t.groups) {
    str(g.id);
    s << g.cols << "," << g.rows << "," << g.x << "," << g.y << "," << g.width
      << "," << g.height << "," << g.halo_x << "," << g.halo_y << ","
      << int(g.orientation) << "," << g.allowed_orientations << "," << g.fixed
      << "," << g.freeze_topology << ";";
    for (const auto &m : g.members) {
      str(m.name);
      s << m.dx << "," << m.dy << "," << int(m.orientation) << ";";
      str(m.key_parent);
      str(m.key_signature);
      s << m.key_ordinal << ";";
    }
  }
  for (const auto &a : t.anchors) {
    str(a.id);
    str(a.group);
    s << int(a.side) << "," << a.tolerance << ";";
    strength(a.strength);
  }
  for (const auto &e : t.order) {
    str(e.id);
    str(e.from);
    str(e.to);
    s << int(e.axis) << ";";
    strength(e.strength);
  }
  for (const auto &f : t.fixed) {
    str(f.name);
    s << f.x << "," << f.y << "," << int(f.orientation) << ";";
    strength(f.strength);
  }
  for (const auto &b : t.blockages) {
    str(b.id);
    s << b.rect.xl << "," << b.rect.yl << "," << b.rect.xh << "," << b.rect.yh
      << ";";
    strength(b.strength);
  }
  for (const auto &f : t.dataflow) {
    str(f.from);
    str(f.to);
    s << f.weight << ";";
  }
  for (const auto &v : t.realized.plus)
    str(v);
  s << "|";
  for (const auto &v : t.realized.minus)
    str(v);
  for (const auto &e : t.realized_graph) {
    str(e.id);
    str(e.from);
    str(e.to);
    s << int(e.axis) << ";";
  }
  std::uint64_t h = 14695981039346656037ULL;
  for (unsigned char c : s.str())
    h = (h ^ c) * 1099511628211ULL;
  std::ostringstream out;
  out << "fnv1a64:" << std::hex << h;
  return out.str();
}
Topology remap_unstable(const Topology &t, const Design &d,
                        std::span<const NodeInfo> info) {
  d.validate();
  const auto ids = index(info, d);
  using Key = std::pair<std::string, std::string>;
  std::map<Key, std::vector<std::string>> candidates;
  std::map<std::string, Key> keys;
  for (const auto &[name, i] : ids)
    if (info[i].macro && !info[i].stable) {
      std::ostringstream sig;
      sig << std::setprecision(17) << info[i].master << ":" << d.node_size_x[i]
          << "," << d.node_size_y[i] << ":";
      std::vector<std::pair<double, double>> pins;
      for (int k = d.flat_node2pin_start[i]; k < d.flat_node2pin_start[i + 1];
           ++k) {
        int p = d.flat_node2pin[k];
        auto o =
            d.node_orientation.empty() ? Orientation::N : d.node_orientation[i];
        auto xy = orient_pin(d.pin_offset_x[p], d.pin_offset_y[p],
                             d.node_size_x[i], d.node_size_y[i], o);
        pins.push_back(xy);
      }
      std::sort(pins.begin(), pins.end());
      for (auto [x, y] : pins)
        sig << x << "," << y << ";";
      Key key{info[i].parent, sig.str()};
      keys[name] = key;
      candidates[key].push_back(name);
    }
  Topology out = t;
  std::map<std::string, std::string> aliases;
  for (auto &g : out.groups)
    for (auto &m : g.members) {
      if (!ids.contains(m.name) && m.key_ordinal >= 0) {
        auto found = candidates.find({m.key_parent, m.key_signature});
        if (found != candidates.end() &&
            size_t(m.key_ordinal) < found->second.size()) {
          auto name = found->second[m.key_ordinal];
          aliases[m.name] = name;
          out.diagnostics.push_back(
              "unstable macro matched by parent/signature/ordinal: " + m.name +
              " -> " + name);
          m.name = name;
        }
      }
      if (auto key = keys.find(m.name); key != keys.end()) {
        m.key_parent = key->second.first;
        m.key_signature = key->second.second;
        const auto &v = candidates.at(key->second);
        m.key_ordinal = std::find(v.begin(), v.end(), m.name) - v.begin();
      }
    }
  for (auto &f : out.fixed)
    if (aliases.contains(f.name))
      f.name = aliases.at(f.name);
  for (auto &f : out.dataflow) {
    if (aliases.contains(f.from))
      f.from = aliases.at(f.from);
    if (aliases.contains(f.to))
      f.to = aliases.at(f.to);
  }
  return out;
}
std::pair<double, double> member_position(const Group &g, const Member &m,
                                          double w, double h) {
  const unsigned bits = unsigned(g.orientation);
  return {g.x + ((bits & 1) ? g.width - m.dx - w : m.dx),
          g.y + ((bits & 2) ? g.height - m.dy - h : m.dy)};
}
unsigned group_orientations(const Group &g, std::span<const NodeInfo> info,
                            const std::map<std::string, int> &names) {
  unsigned result = 0;
  for (unsigned o = 0; o < 4; ++o) {
    bool ok = true;
    for (const auto &m : g.members)
      ok &= bool(info[names.at(m.name)].allowed_orientations &
                 (1U << (o ^ unsigned(m.orientation))));
    if (ok)
      result |= 1U << o;
  }
  return result;
}
unsigned group_orientations(const Group &g, std::span<const NodeInfo> info) {
  std::map<std::string, int> names;
  for (int i = 0; i < int(info.size()); ++i)
    names[info[i].name] = i;
  return group_orientations(g, info, names);
}
void validate(const Topology &t, const Design &d,
              std::span<const NodeInfo> info) {
  d.validate();
  const auto ids = index(info, d);
  need(t.schema_version == 1, "unsupported topology schema");
  std::set<std::string> groups, members, keys;
  for (const auto &g : t.groups) {
    need(!g.id.empty() && groups.insert(g.id).second,
         "duplicate/empty group id");
    need(!g.members.empty() && g.cols > 0 && g.rows > 0, "empty group/tiling");
    for (double v : {g.x, g.y, g.width, g.height, g.halo_x, g.halo_y})
      need(std::isfinite(v), "nonfinite group geometry");
    need(g.width > 0 && g.height > 0 && g.halo_x >= 0 && g.halo_y >= 0,
         "invalid group dimensions");
    need(g.allowed_orientations > 0 && g.allowed_orientations <= 15 &&
             valid_orientation(g.orientation) &&
             (g.allowed_orientations & (1U << unsigned(g.orientation))),
         "group orientation disallowed");
    for (size_t a = 0; a < g.members.size(); ++a) {
      const auto &m = g.members[a];
      need(ids.contains(m.name) && info[ids.at(m.name)].macro,
           "unknown macro member: " + m.name);
      need(members.insert(m.name).second, "macro belongs to multiple groups");
      int i = ids.at(m.name);
      need(std::isfinite(m.dx) && std::isfinite(m.dy) && m.dx >= 0 &&
               m.dy >= 0 && m.dx + d.node_size_x[i] <= g.width + 1e-6 &&
               m.dy + d.node_size_y[i] <= g.height + 1e-6,
           "member outside group");
      need(std::abs(m.dx / d.layout.site_width -
                    std::round(m.dx / d.layout.site_width)) < 1e-6 &&
               std::abs(m.dy / d.layout.row_height -
                        std::round(m.dy / d.layout.row_height)) < 1e-6,
           "member offset off pitch");
      need(valid_orientation(m.orientation) &&
               (info[i].allowed_orientations &
                (1U << unsigned(compose(g.orientation, m.orientation)))),
           "member orientation disallowed");
      need(g.members.size() == 1 || i < d.num_movable,
           "fixed macro cannot be in a movable array");
      auto hx = [&](int n) {
        return d.node_halo_x.empty()
                   ? 0
                   : std::ceil(d.node_halo_x[n] / d.layout.site_width - 1e-9) *
                         d.layout.site_width;
      };
      auto hy = [&](int n) {
        return d.node_halo_y.empty()
                   ? 0
                   : std::ceil(d.node_halo_y[n] / d.layout.row_height - 1e-9) *
                         d.layout.row_height;
      };
      for (size_t b = 0; b < a; ++b) {
        const auto &n = g.members[b];
        int j = ids.at(n.name);
        need(m.dx + d.node_size_x[i] + hx(i) + hx(j) <= n.dx + 1e-6 ||
                 n.dx + d.node_size_x[j] + hx(i) + hx(j) <= m.dx + 1e-6 ||
                 m.dy + d.node_size_y[i] + hy(i) + hy(j) <= n.dy + 1e-6 ||
                 n.dy + d.node_size_y[j] + hy(i) + hy(j) <= m.dy + 1e-6,
             "overlapping haloed members");
      }
    }
  }
  for (const auto &a : t.anchors) {
    need(int(a.side) >= 0 && int(a.side) <= 3, "invalid anchor side");
    need(groups.contains(a.group), "anchor names unknown group");
    need(std::isfinite(a.tolerance) && a.tolerance >= 0 &&
             std::isfinite(a.strength.weight) && a.strength.weight >= 0,
         "invalid anchor");
    need(!a.id.empty() && keys.insert(a.id).second, "duplicate constraint id");
  }
  for (const auto &e : t.order) {
    need(e.axis == Axis::x || e.axis == Axis::y, "invalid order axis");
    need(groups.contains(e.from) && groups.contains(e.to) && e.from != e.to,
         "invalid order endpoints");
    need(std::isfinite(e.strength.weight) && e.strength.weight >= 0,
         "invalid order weight");
    need(!e.id.empty() && keys.insert(e.id).second, "duplicate constraint id");
  }
  std::set<std::string> fixed_names, blockage_ids;
  for (const auto &f : t.fixed)
    need(fixed_names.insert(f.name).second && ids.contains(f.name) &&
             info[ids.at(f.name)].macro && std::isfinite(f.x) &&
             std::isfinite(f.y) && valid_orientation(f.orientation) &&
             std::isfinite(f.strength.weight) && f.strength.weight >= 0,
         "invalid fixed macro");
  for (const auto &b : t.blockages)
    need(blockage_ids.insert(b.id).second && std::isfinite(b.rect.xl) &&
             std::isfinite(b.rect.yl) && std::isfinite(b.rect.xh) &&
             std::isfinite(b.rect.yh) && b.rect.xl < b.rect.xh &&
             b.rect.yl < b.rect.yh && !b.id.empty() &&
             std::isfinite(b.strength.weight) && b.strength.weight >= 0,
         "invalid blockage");
  for (const auto &f : t.dataflow)
    need(!f.from.empty() && !f.to.empty() && std::isfinite(f.weight) &&
             f.weight >= 0,
         "invalid flow weight/endpoint");
}
Topology merge(const Topology &base, const Topology &requested, const Design &d,
               std::span<const NodeInfo> info) {
  need(requested.schema_version == 1, "unsupported topology schema");
  const auto overlay = remap_unstable(requested, d, info);
  Topology out = base;
  out.diagnostics.insert(out.diagnostics.end(), overlay.diagnostics.begin(),
                         overlay.diagnostics.end());
  auto ids = index(info, d);
  auto replace = [](auto &dst, const auto &src, auto key) {
    for (const auto &v : src) {
      auto it = std::find_if(dst.begin(), dst.end(),
                             [&](const auto &a) { return key(a) == key(v); });
      if (it == dst.end())
        dst.push_back(v);
      else
        *it = v;
    }
  };
  replace(out.groups, overlay.groups, [](const Group &g) { return g.id; });
  replace(out.anchors, overlay.anchors, [](const Anchor &a) { return a.id; });
  replace(out.order, overlay.order, [](const OrderEdge &e) { return e.id; });
  replace(out.fixed, overlay.fixed, [](const Fixed &f) { return f.name; });
  replace(out.blockages, overlay.blockages,
          [](const Blockage &b) { return b.id; });
  for (auto &g : out.groups)
    std::erase_if(g.members, [&](const Member &m) {
      bool gone = !ids.contains(m.name) || !info[ids.at(m.name)].macro;
      if (gone)
        out.diagnostics.push_back("dropped stale macro: " + m.name);
      return gone;
    });
  std::erase_if(out.groups, [](const Group &g) { return g.members.empty(); });
  std::set<std::string> groupids;
  for (const auto &g : out.groups)
    groupids.insert(g.id);
  std::erase_if(out.anchors,
                [&](const Anchor &a) { return !groupids.contains(a.group); });
  std::erase_if(out.order, [&](const OrderEdge &e) {
    return !groupids.contains(e.from) || !groupids.contains(e.to);
  });
  std::erase_if(out.fixed,
                [&](const Fixed &f) { return !ids.contains(f.name); });
  if (!overlay.realized_graph.empty())
    out.realized_graph = overlay.realized_graph;
  if (!overlay.realized.plus.empty())
    out.realized = overlay.realized;
  if (!overlay.dataflow.empty())
    out.dataflow = overlay.dataflow;
  out.macro_fingerprint = fingerprint(d, info);
  out.netlist_hash = fingerprint(d, info, true);
  return remap_unstable(out, d, info);
}
} // namespace DreamPlace::native::macro
