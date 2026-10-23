// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rigid.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace DreamPlace::native {
namespace {
void require(bool ok, const char *reason) {
  if (!ok)
    throw std::invalid_argument(reason);
}
double halo(const std::vector<double> &v, int i, double pitch) {
  return v.empty() ? 0 : std::ceil(v[i] / pitch - 1e-9) * pitch;
}
bool close(double a, double b) {
  return std::abs(a - b) <= 1e-9 * std::max({1., std::abs(a), std::abs(b)});
}
struct Box {
  double xl, yl, xh, yh;
  double area() const { return (xh - xl) * (yh - yl); }
};
bool intersects(const Box &a, const Box &b) {
  return a.xl < b.xh - 1e-9 && b.xl < a.xh - 1e-9 && a.yl < b.yh - 1e-9 &&
         b.yl < a.yh - 1e-9;
}
Box padded_box(const Design &d, int i, double x, double y) {
  const auto hx = halo(d.node_halo_x, i, d.layout.site_width);
  const auto hy = halo(d.node_halo_y, i, d.layout.row_height);
  return {x - hx, y - hy, x + d.node_size_x[i] + hx, y + d.node_size_y[i] + hy};
}
template <class T>
void select(std::vector<T> &out, const std::vector<T> &in,
            const std::vector<int> &prototype) {
  if (in.empty()) {
    out.clear();
    return;
  }
  out.clear();
  out.reserve(prototype.size());
  for (int i : prototype)
    out.push_back(in.at(i));
}
void select_nodes(Design &out, const Design &in,
                  const std::vector<int> &prototype) {
  select(out.node_x, in.node_x, prototype);
  select(out.node_y, in.node_y, prototype);
  select(out.node_size_x, in.node_size_x, prototype);
  select(out.node_size_y, in.node_size_y, prototype);
  select(out.node_macro, in.node_macro, prototype);
  select(out.node_align_class, in.node_align_class, prototype);
  select(out.node_snap_class, in.node_snap_class, prototype);
  select(out.node_blocked_capacity_h, in.node_blocked_capacity_h, prototype);
  select(out.node_blocked_capacity_v, in.node_blocked_capacity_v, prototype);
  select(out.node_legalize_weights, in.node_legalize_weights, prototype);
  select(out.node_legalize_last, in.node_legalize_last, prototype);
  select(out.node_legalize_anchor_x, in.node_legalize_anchor_x, prototype);
  select(out.node_legalize_anchor_y, in.node_legalize_anchor_y, prototype);
  select(out.node_halo_x, in.node_halo_x, prototype);
  select(out.node_halo_y, in.node_halo_y, prototype);
  select(out.node_orientation, in.node_orientation, prototype);
  select(out.node_orientation_fixed, in.node_orientation_fixed, prototype);
  select(out.node_group_id, in.node_group_id, prototype);
  select(out.node_boxes, in.node_boxes, prototype);
  select(out.node_inflatable, in.node_inflatable, prototype);
  select(out.node_centroid_color, in.node_centroid_color, prototype);
  select(out.centroid_node_area, in.centroid_node_area, prototype);
}
int next_class(const std::vector<int> &values, const RigidMacroMap &map,
               bool snap) {
  int result =
      values.empty() ? -1 : *std::max_element(values.begin(), values.end());
  for (const auto &[shape, classes] : map.shape_classes) {
    (void)shape;
    result = std::max(result, snap ? classes.second : classes.first);
  }
  require(result < std::numeric_limits<int>::max(),
          "macro shape class overflow");
  return result + 1;
}
} // namespace

void rebind_rigid_pins(const Design &source, Design &destination,
                       std::span<const RigidPinFrame> frames) {
  source.validate();
  require(frames.size() == source.node_x.size(),
          "rigid pin frame size mismatch");
  auto nodes = source.pin2node;
  auto x = source.pin_offset_x, y = source.pin_offset_y;
  for (size_t i = 0; i < frames.size(); ++i) {
    const auto &f = frames[i];
    require(f.node >= 0 && f.node < int(destination.node_x.size()) &&
                std::isfinite(f.x) && std::isfinite(f.y) &&
                valid_orientation(f.transform),
            "invalid rigid pin frame");
    for (int k = source.flat_node2pin_start[i];
         k < source.flat_node2pin_start[i + 1]; ++k) {
      const int p = source.flat_node2pin[k];
      auto offset = orient_pin(x[p], y[p], source.node_size_x[i],
                               source.node_size_y[i], f.transform);
      nodes[p] = f.node;
      x[p] = f.x + offset.first;
      y[p] = f.y + offset.second;
      require(std::isfinite(x[p]) && std::isfinite(y[p]),
              "rigid pin offset overflow");
    }
  }
  Design candidate = destination;
  candidate.pin2node = std::move(nodes);
  candidate.pin_offset_x = std::move(x);
  candidate.pin_offset_y = std::move(y);
  candidate.flat_net2pin = source.flat_net2pin;
  candidate.flat_net2pin_start = source.flat_net2pin_start;
  candidate.net_weights = source.net_weights;
  candidate.net_mask = source.net_mask;
  candidate.build_maps();
  candidate.validate();
  destination = std::move(candidate);
}

MacroMergeResult merge_macros(Design &design, RigidMacroMap &mapping,
                              std::span<const MacroMerge> requests) {
  MacroMergeResult result;
  try {
    design.validate();
    require(!requests.empty(), "empty macro merge batch");
    require(design.node_pad_x.empty() && design.node_pad_y.empty(),
            "macro merge requires physical, unpadded geometry");
    auto map = mapping;
    const int n = design.node_x.size();
    if (map.source_to_node.empty()) {
      require(map.groups.empty(), "missing rigid source mapping");
      map.source_to_node.resize(n);
      std::iota(map.source_to_node.begin(), map.source_to_node.end(), 0);
    }
    for (int i : map.source_to_node)
      require(i >= 0 && i < n, "invalid rigid source mapping");
    std::set<int> old_group_nodes;
    for (const auto &group : map.groups) {
      require(group.node >= 0 && group.node < design.num_movable &&
                  old_group_nodes.insert(group.node).second && group.id > 0 &&
                  (design.node_orientation.empty() ||
                   design.node_orientation[group.node] == Orientation::N),
              "invalid or reoriented rigid group");
      for (const auto &m : group.members)
        require(m.source >= 0 && m.source < int(map.source_to_node.size()) &&
                    map.source_to_node[m.source] == group.node &&
                    std::isfinite(m.x) && std::isfinite(m.y),
                "invalid rigid group member mapping");
    }
    std::vector<int> owner(n, -1), remap(n, -1), prototype;
    std::vector<Box> unions;
    std::vector<double> halos_x, halos_y;
    for (int g = 0; g < int(requests.size()); ++g) {
      const auto &request = requests[g];
      require(request.members.size() >= 2 &&
                  request.poses.size() == request.members.size(),
              "macro merge needs members and matching poses");
      Box outer{INFINITY, INFINITY, -INFINITY, -INFINITY};
      double area = 0, hx = 0, hy = 0;
      std::vector<Box> boxes;
      for (size_t k = 0; k < request.members.size(); ++k) {
        const int i = request.members[k];
        require(i >= 0 && i < design.num_movable &&
                    !design.node_macro.empty() && design.node_macro[i],
                "macro merge member is not a movable macro");
        require(owner[i] < 0, "macro merge members are not independent");
        owner[i] = g;
        const auto [x, y] = request.poses[k];
        require(std::isfinite(x) && std::isfinite(y),
                "nonfinite macro merge pose");
        const auto bx = halo(design.node_halo_x, i, design.layout.site_width);
        const auto by = halo(design.node_halo_y, i, design.layout.row_height);
        if (k == 0) {
          hx = bx;
          hy = by;
        }
        require(close(hx, bx) && close(hy, by),
                "macro merge has incompatible outside halos");
        const auto box = padded_box(design, i, x, y);
        for (const auto &other : boxes)
          require(!intersects(box, other),
                  "macro merge violates member halo clearance");
        boxes.push_back(box);
        area += box.area();
        outer.xl = std::min(outer.xl, box.xl);
        outer.yl = std::min(outer.yl, box.yl);
        outer.xh = std::max(outer.xh, box.xh);
        outer.yh = std::max(outer.yh, box.yh);
      }
      require(close(area, outer.area()),
              "macro merge outer rectangle has holes");
      require(outer.xl >= design.layout.xl && outer.yl >= design.layout.yl &&
                  outer.xh <= design.layout.xh && outer.yh <= design.layout.yh,
              "macro merge leaves the core");
      for (int i = design.num_movable; i < n; ++i)
        if (design.node_size_x[i] > 0 && design.node_size_y[i] > 0)
          require(!intersects(outer, padded_box(design, i, design.node_x[i],
                                                design.node_y[i])),
                  "macro merge overlaps a fixed obstacle");
      for (const auto &other : unions)
        require(!intersects(outer, other),
                "macro merge batch rectangles overlap");
      unions.push_back(outer);
      halos_x.push_back(hx);
      halos_y.push_back(hy);
    }
    // Preserve source order for all untouched movable nodes and all fixed
    // nodes. Group nodes precede them, so the movable/fixed partition survives.
    for (const auto &r : requests)
      prototype.push_back(r.members.front());
    for (int i = 0; i < design.num_movable; ++i)
      if (owner[i] < 0)
        prototype.push_back(i);
    const int movable = prototype.size();
    for (int i = design.num_movable; i < n; ++i)
      prototype.push_back(i);
    for (int i = 0; i < int(prototype.size()); ++i)
      remap[prototype[i]] = i;
    for (int g = 0; g < int(requests.size()); ++g)
      for (int i : requests[g].members)
        remap[i] = g;
    Design out = design;
    select_nodes(out, design, prototype);
    out.num_movable = movable;
    out.node_orientation.resize(prototype.size(), Orientation::N);
    out.node_orientation_fixed.resize(prototype.size(), 0);
    out.node_group_id.resize(prototype.size(), 0);
    out.node_align_class.resize(prototype.size(), -1);
    out.node_snap_class.resize(prototype.size(), -1);
    std::vector<RigidPinFrame> frames(n);
    for (int i = 0; i < n; ++i)
      frames[i].node = remap[i];
    std::vector<RigidGroup> groups;
    for (const auto &group : map.groups)
      if (owner.at(group.node) < 0) {
        auto kept = group;
        kept.node = remap[group.node];
        groups.push_back(std::move(kept));
      }
    for (int g = 0; g < int(requests.size()); ++g) {
      require(map.next_id > 0 &&
                  map.next_id < std::numeric_limits<std::uint64_t>::max(),
              "rigid group identity overflow");
      RigidGroup group{map.next_id++, g, {}};
      const auto &box = unions[g];
      const int first = requests[g].members.front();
      const double tile_w = design.node_size_x[first] + 2 * halos_x[g];
      const double tile_h = design.node_size_y[first] + 2 * halos_y[g];
      const double columns = std::round((box.xh - box.xl) / tile_w);
      const double rows = std::round((box.yh - box.yl) / tile_h);
      require(columns >= 1 && rows >= 1 &&
                  close(columns * tile_w, box.xh - box.xl) &&
                  close(rows * tile_h, box.yh - box.yl),
              "macro merge is not a congruent rectangular tiling");
      out.node_x[g] = box.xl + halos_x[g];
      out.node_y[g] = box.yl + halos_y[g];
      out.node_size_x[g] = columns * tile_w - 2 * halos_x[g];
      out.node_size_y[g] = rows * tile_h - 2 * halos_y[g];
      out.node_orientation[g] = Orientation::N;
      out.node_orientation_fixed[g] = 1;
      out.node_group_id[g] = group.id;
      const auto shape = std::pair{columns * tile_w, rows * tile_h};
      if (!map.shape_classes.contains(shape))
        map.shape_classes.emplace(
            shape, std::pair{next_class(design.node_align_class, map, false),
                             next_class(design.node_snap_class, map, true)});
      out.node_align_class[g] = map.shape_classes.at(shape).first;
      out.node_snap_class[g] = map.shape_classes.at(shape).second;
      if (!out.node_inflatable.empty())
        out.node_inflatable[g] = 0;
      if (!out.node_legalize_anchor_x.empty()) {
        out.node_legalize_anchor_x[g] = out.node_x[g];
        out.node_legalize_anchor_y[g] = out.node_y[g];
      }
      if (!out.node_centroid_color.empty())
        out.node_centroid_color[g] = -1;
      if (!out.centroid_node_area.empty())
        out.centroid_node_area[g] = 0;
      PositionBox bounds{out.layout.xl + halos_x[g], out.layout.yl + halos_y[g],
                         out.layout.xh - out.node_size_x[g] - halos_x[g],
                         out.layout.yh - out.node_size_y[g] - halos_y[g]};
      for (size_t k = 0; k < requests[g].members.size(); ++k) {
        const int i = requests[g].members[k];
        const auto [x, y] = requests[g].poses[k];
        auto &f = frames[i];
        require(close(design.node_size_x[i] + 2 * halos_x[g], tile_w) &&
                    close(design.node_size_y[i] + 2 * halos_y[g], tile_h),
                "macro merge requires equal padded shapes");
        f.x = std::round((x - out.node_x[g]) / tile_w) * tile_w;
        f.y = std::round((y - out.node_y[g]) / tile_h) * tile_h;
        require(close(x, out.node_x[g] + f.x) && close(y, out.node_y[g] + f.y),
                "macro merge member is off the rigid tiling");
        if (!out.centroid_node_area.empty())
          out.centroid_node_area[g] += design.centroid_node_area[i];
        if (!out.node_legalize_weights.empty())
          out.node_legalize_weights[g] = std::max(
              out.node_legalize_weights[g], design.node_legalize_weights[i]);
        if (!design.node_boxes.empty()) {
          const auto &b = design.node_boxes[i];
          bounds.xl = std::max(bounds.xl, b.xl - f.x);
          bounds.yl = std::max(bounds.yl, b.yl - f.y);
          bounds.xh = std::min(bounds.xh, b.xh - f.x);
          bounds.yh = std::min(bounds.yh, b.yh - f.y);
        }
        const auto old =
            std::find_if(map.groups.begin(), map.groups.end(),
                         [&](const auto &v) { return v.node == i; });
        if (old != map.groups.end()) {
          for (const auto &m : old->members)
            group.members.push_back({m.source, f.x + m.x, f.y + m.y});
        } else {
          auto source = std::find(map.source_to_node.begin(),
                                  map.source_to_node.end(), i);
          require(source != map.source_to_node.end(),
                  "macro merge has no source identity");
          group.members.push_back(
              {int(source - map.source_to_node.begin()), f.x, f.y});
        }
      }
      require(bounds.xl <= bounds.xh && bounds.yl <= bounds.yh &&
                  out.node_x[g] >= bounds.xl && out.node_x[g] <= bounds.xh &&
                  out.node_y[g] >= bounds.yl && out.node_y[g] <= bounds.yh,
              "macro merge violates member position bounds");
      if (!out.node_boxes.empty())
        out.node_boxes[g] = bounds;
      std::sort(
          group.members.begin(), group.members.end(),
          [](const auto &a, const auto &b) { return a.source < b.source; });
      result.group_ids.push_back(group.id);
      groups.push_back(std::move(group));
    }
    // Preserve each logical cohesion member as an offset sample. A group may
    // span several source colors; this keeps their energies and accumulates
    // all member forces on the single rigid node.
    for (auto &group : out.cohesion_groups) {
      std::vector<std::pair<double, double>> offsets;
      std::vector<double> areas;
      for (size_t k = 0; k < group.nodes.size(); ++k) {
        const int i = group.nodes[k];
        const int j = remap[i];
        const auto old =
            group.offsets.empty() ? std::pair{0., 0.} : group.offsets[k];
        offsets.emplace_back(frames[i].x + design.node_size_x[i] / 2 -
                                 out.node_size_x[j] / 2 + old.first,
                             frames[i].y + design.node_size_y[i] / 2 -
                                 out.node_size_y[j] / 2 + old.second);
        areas.push_back(group.areas.empty()
                            ? design.node_size_x[i] * design.node_size_y[i]
                            : group.areas[k]);
        group.nodes[k] = j;
      }
      group.offsets = std::move(offsets);
      group.areas = std::move(areas);
    }
    for (auto &group : out.centroid_groups) {
      // Resolve colors to logical samples before contraction. A hierarchy
      // group may own only part of a rigid rectangle: neither its full area
      // nor the union center can stand in for those original member moments.
      const size_t explicit_count = group.nodes.size();
      if (group.areas.empty())
        for (int i : group.nodes)
          group.areas.push_back(design.centroid_node_area[i]);
      if (group.offsets.empty())
        group.offsets.resize(explicit_count, {0, 0});
      for (int color : group.colors)
        for (int i = 0; i < n; ++i)
          if (design.node_centroid_color[i] == color) {
            group.nodes.push_back(i);
            group.areas.push_back(design.centroid_node_area[i]);
            group.offsets.emplace_back(0, 0);
          }
      group.colors.clear();
      for (size_t k = 0; k < group.nodes.size(); ++k) {
        const int i = group.nodes[k], j = remap[i];
        group.offsets[k].first +=
            frames[i].x + (design.node_size_x[i] - out.node_size_x[j]) / 2;
        group.offsets[k].second +=
            frames[i].y + (design.node_size_y[i] - out.node_size_y[j]) / 2;
        group.nodes[k] = j;
      }
    }
    out.order_edges.clear();
    for (auto edge : design.order_edges) {
      const auto &a = frames[edge.from], &b = frames[edge.to];
      const double offset_a = edge.axis == Axis::x ? a.x : a.y;
      const double offset_b = edge.axis == Axis::x ? b.x : b.y;
      const double size = edge.axis == Axis::x ? design.node_size_x[edge.from]
                                               : design.node_size_y[edge.from];
      const double group_size = edge.axis == Axis::x ? out.node_size_x[a.node]
                                                     : out.node_size_y[a.node];
      const double separation = size + edge.separation + offset_a - offset_b;
      if (a.node == b.node) {
        require(!edge.strength.hard || separation <= 1e-9,
                "macro merge violates internal hard order");
        continue;
      }
      edge.from = a.node;
      edge.to = b.node;
      edge.separation = separation - group_size;
      out.order_edges.push_back(std::move(edge));
    }
    rebind_rigid_pins(design, out, frames);
    for (int net = 0; net < int(out.net_weights.size()); ++net) {
      int node = -1;
      bool internal = true;
      for (int k = out.flat_net2pin_start[net];
           k < out.flat_net2pin_start[net + 1]; ++k) {
        const int i = out.pin2node[out.flat_net2pin[k]];
        if (node < 0)
          node = i;
        internal &= node == i;
      }
      if (internal)
        out.net_mask[net] = 0;
    }
    out.validate();
    for (auto &i : map.source_to_node)
      i = remap[i];
    map.groups = std::move(groups);
    result.previous_node = prototype;
    std::fill_n(result.previous_node.begin(), requests.size(), -1);
    design = std::move(out);
    mapping = std::move(map);
    result.applied = true;
  } catch (const std::invalid_argument &error) {
    result.reason = error.what();
    result.previous_node.clear();
    result.group_ids.clear();
  }
  return result;
}

Design expand_macros(const Design &grouped, const Design &source,
                     const RigidMacroMap &map) {
  grouped.validate();
  source.validate();
  require(map.source_to_node.size() == source.node_x.size(),
          "rigid expansion source size mismatch");
  require(grouped.pin2node.size() == source.pin2node.size() &&
              grouped.flat_net2pin == source.flat_net2pin &&
              grouped.flat_net2pin_start == source.flat_net2pin_start,
          "rigid expansion netlist identity mismatch");
  Design out = source;
  out.node_group_id.assign(source.node_x.size(), 0);
  std::vector<double> offsets_x(source.node_x.size()),
      offsets_y(source.node_x.size());
  std::vector<bool> member(source.node_x.size());
  for (const auto &group : map.groups) {
    require(group.node >= 0 && group.node < grouped.num_movable &&
                (grouped.node_orientation.empty() ||
                 grouped.node_orientation[group.node] == Orientation::N),
            "rigid group orientation is frozen");
    for (const auto &m : group.members) {
      require(m.source >= 0 && m.source < source.num_movable &&
                  !member[m.source] &&
                  map.source_to_node[m.source] == group.node,
              "invalid rigid expansion member");
      member[m.source] = true;
      offsets_x[m.source] = m.x;
      offsets_y[m.source] = m.y;
      out.node_group_id[m.source] = group.id;
    }
  }
  if (!source.node_orientation.empty() || !grouped.node_orientation.empty())
    out.node_orientation.resize(source.node_x.size(), Orientation::N);
  for (size_t i = 0; i < source.node_x.size(); ++i) {
    const int node = map.source_to_node[i];
    require(node >= 0 && node < int(grouped.node_x.size()),
            "invalid rigid expansion node");
    out.node_x[i] = grouped.node_x[node] + offsets_x[i];
    out.node_y[i] = grouped.node_y[node] + offsets_y[i];
    if (!member[i]) {
      out.node_size_x[i] = grouped.node_size_x[node];
      out.node_size_y[i] = grouped.node_size_y[node];
      if (!out.node_orientation.empty())
        out.node_orientation[i] = grouped.node_orientation.empty()
                                      ? Orientation::N
                                      : grouped.node_orientation[node];
    }
    if (i >= size_t(source.num_movable) && source.node_size_x[i] > 0 &&
        source.node_size_y[i] > 0)
      require(out.node_x[i] == source.node_x[i] &&
                  out.node_y[i] == source.node_y[i],
              "rigid expansion moved a fixed node");
  }
  for (size_t p = 0; p < source.pin2node.size(); ++p) {
    const int i = source.pin2node[p];
    require(grouped.pin2node[p] == map.source_to_node[i],
            "rigid expansion changed pin ownership");
    out.pin_offset_x[p] = grouped.pin_offset_x[p] - offsets_x[i];
    out.pin_offset_y[p] = grouped.pin_offset_y[p] - offsets_y[i];
    if (member[i])
      require(close(out.pin_offset_x[p], source.pin_offset_x[p]) &&
                  close(out.pin_offset_y[p], source.pin_offset_y[p]),
              "rigid group changed a frozen member pin");
  }
  out.net_weights = grouped.net_weights;
  out.build_maps();
  out.validate();
  return out;
}
Design contract_macros(const Design &expanded, const Design &grouped,
                       const Design &source, const RigidMacroMap &map) {
  expanded.validate();
  const auto before = expand_macros(grouped, source, map);
  require(expanded.num_movable == before.num_movable &&
              expanded.node_x.size() == before.node_x.size() &&
              expanded.pin2node == before.pin2node &&
              expanded.flat_net2pin == before.flat_net2pin &&
              expanded.flat_net2pin_start == before.flat_net2pin_start,
          "macro magnet host update changed netlist topology");
  Design out = grouped;
  std::vector<bool> member(expanded.node_x.size());
  for (const auto &group : map.groups) {
    require(!group.members.empty(), "empty rigid group");
    const auto &first = group.members.front();
    out.node_x[group.node] = expanded.node_x[first.source] - first.x;
    out.node_y[group.node] = expanded.node_y[first.source] - first.y;
    for (const auto &m : group.members) {
      const int i = m.source;
      member[i] = true;
      require(close(expanded.node_x[i], out.node_x[group.node] + m.x) &&
                  close(expanded.node_y[i], out.node_y[group.node] + m.y) &&
                  expanded.node_size_x[i] == before.node_size_x[i] &&
                  expanded.node_size_y[i] == before.node_size_y[i] &&
                  (expanded.node_orientation.empty()
                       ? Orientation::N
                       : expanded.node_orientation[i]) ==
                      (before.node_orientation.empty()
                           ? Orientation::N
                           : before.node_orientation[i]),
              "host update broke a rigid macro group");
    }
  }
  for (size_t i = 0; i < expanded.node_x.size(); ++i) {
    if (member[i])
      continue;
    const int node = map.source_to_node[i];
    out.node_x[node] = expanded.node_x[i];
    out.node_y[node] = expanded.node_y[i];
    out.node_size_x[node] = expanded.node_size_x[i];
    out.node_size_y[node] = expanded.node_size_y[i];
    if (!out.node_orientation.empty())
      out.node_orientation[node] = expanded.node_orientation.empty()
                                       ? Orientation::N
                                       : expanded.node_orientation[i];
  }
  for (size_t p = 0; p < expanded.pin2node.size(); ++p) {
    const int i = expanded.pin2node[p];
    if (member[i]) {
      require(expanded.pin_offset_x[p] == before.pin_offset_x[p] &&
                  expanded.pin_offset_y[p] == before.pin_offset_y[p],
              "host update changed a frozen member pin");
    } else {
      out.pin_offset_x[p] = expanded.pin_offset_x[p];
      out.pin_offset_y[p] = expanded.pin_offset_y[p];
    }
  }
  out.net_weights = expanded.net_weights;
  out.validate();
  // Also enforce fixed obstacles and the inverse pin transform before commit.
  (void)expand_macros(out, source, map);
  return out;
}
} // namespace DreamPlace::native
