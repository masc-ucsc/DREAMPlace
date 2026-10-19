// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/parameters.hpp"
#include "dreamplace/native/placement.hpp"
#include <map>
#include <optional>
#include <span>
#include <string>
namespace DreamPlace::native::macro {
struct NodeInfo {
  std::string name, parent, master;
  bool macro = false, sequential = false, io = false, stable = true;
  unsigned allowed_orientations = 15;
};
struct Member {
  std::string name;
  double dx = 0, dy = 0;
  Orientation orientation = Orientation::N;
  // Diagnostic fallback only for unnamed instances: parent, pin/size signature,
  // and name-sorted ordinal within the equivalent siblings.
  std::string key_parent, key_signature;
  int key_ordinal = -1;
};
struct Group {
  std::string id;
  std::vector<Member> members;
  int cols = 1, rows = 1;
  double x = 0, y = 0, width = 0, height = 0, halo_x = 0, halo_y = 0;
  Orientation orientation = Orientation::N;
  unsigned allowed_orientations = 15;
  bool fixed = false, freeze_topology = false;
};
struct Anchor {
  std::string id, group;
  Side side = Side::left;
  double tolerance = 0;
  Strength strength;
};
struct OrderEdge {
  std::string id, from, to;
  Axis axis = Axis::x;
  Strength strength;
};
struct Fixed {
  std::string name;
  double x = 0, y = 0;
  Orientation orientation = Orientation::N;
  Strength strength{true, 1};
};
struct Blockage {
  std::string id;
  PositionBox rect;
  Strength strength{true, 1};
};
struct Flow {
  std::string from, to;
  double weight = 0;
};
struct SequencePair {
  std::vector<std::string> plus, minus;
};
struct Separation {
  std::string from, to, id;
  Axis axis = Axis::x;
};
struct Topology {
  int schema_version = 1;
  std::string design, macro_fingerprint, netlist_hash,
      provider = "native_hierarchy";
  std::vector<Group> groups;
  std::vector<Anchor> anchors;
  std::vector<OrderEdge> order;
  std::vector<Fixed> fixed;
  std::vector<Blockage> blockages;
  std::vector<Flow> dataflow;
  SequencePair realized;
  std::vector<Separation> realized_graph;
  std::vector<std::string> diagnostics;
};
std::string fingerprint(const Design &, std::span<const NodeInfo>,
                        bool full_netlist = false);
std::string topology_key(const Topology &);
Topology remap_unstable(const Topology &, const Design &,
                        std::span<const NodeInfo>);
void validate(const Topology &, const Design &, std::span<const NodeInfo>);
// Overlay entries replace matching stable ids. Missing macros are pruned with
// diagnostics; native seeding fills newly introduced leaves afterwards.
Topology merge(const Topology &base, const Topology &overlay, const Design &,
               std::span<const NodeInfo>);
std::pair<double, double> member_position(const Group &, const Member &,
                                          double width, double height);
unsigned group_orientations(const Group &, std::span<const NodeInfo>);
// Same, for callers that already hold a name -> info index map. The overload
// above rebuilds that map over every node, which is O(nodes) per group.
unsigned group_orientations(const Group &, std::span<const NodeInfo>,
                            const std::map<std::string, int> &names);
} // namespace DreamPlace::native::macro
