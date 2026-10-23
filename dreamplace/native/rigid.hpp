// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"
#include <map>
#include <span>

namespace DreamPlace::native {
struct RigidPinFrame {
  int node = -1;
  double x = 0, y = 0; // member origin relative to the destination node
  Orientation transform = Orientation::N; // transform already-oriented pins
};
// Shared by plain magnetization and H2D abstraction. Pins keep their identities
// and net CSR; only their owner and member-to-group coordinate frame change.
// A malformed frame leaves the destination untouched.
void rebind_rigid_pins(const Design &source, Design &destination,
                       std::span<const RigidPinFrame> frames);

struct RigidMember {
  int source = -1;
  double x = 0, y = 0;
};
struct RigidGroup {
  std::uint64_t id = 0;
  int node = -1;
  std::vector<RigidMember> members; // flat source identities, never nested
};
struct RigidMacroMap {
  std::vector<int> source_to_node;
  std::vector<RigidGroup> groups;
  std::uint64_t next_id = 1;
  std::map<std::pair<double, double>, std::pair<int, int>> shape_classes;
};
struct MacroMerge {
  std::vector<int> members;
  std::vector<std::pair<double, double>> poses; // physical lower-left origins
};
struct MacroMergeResult {
  bool applied = false;
  std::string reason;
  std::vector<int> previous_node; // new node -> old node; -1 for a new group
  std::vector<std::uint64_t> group_ids;
};
// Atomic, physical-coordinate replacement. All members must be movable macros
// and their equally padded rectangles must tile the merged outer rectangle.
// Cells in interior channels are left in place for the following GP relaxation.
MacroMergeResult merge_macros(Design &, RigidMacroMap &,
                              std::span<const MacroMerge>);
// Expand against the host's current member netlist. All pin/net identities and
// original member orientations survive; grouped nodes cannot flip. The caller
// receives a new design, so failed audits cannot partially write back members.
Design expand_macros(const Design &grouped, const Design &source,
                     const RigidMacroMap &);
// Accept host changes only when every member still describes the same rigid
// rectangle. This does not re-run merge placement checks: GP may temporarily
// overlap another movable or fixed obstacle before legalization.
Design contract_macros(const Design &expanded, const Design &grouped,
                       const Design &source, const RigidMacroMap &);
} // namespace DreamPlace::native
