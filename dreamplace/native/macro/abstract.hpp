// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/topology.hpp"
namespace DreamPlace::native::macro {
struct Abstracted {
  Design design;
  std::vector<NodeInfo> info;
  std::vector<std::string> source_names;
  std::vector<int> source_to_node, source_of_node;
  std::map<std::string, int> group_node;
};
Abstracted abstract_design(const Design &, std::span<const NodeInfo>,
                           const Topology &, const Parameters &,
                           bool freeze_groups = false);
// Transactional writeback, with pin transforms and rigid member offsets.
void expand(const Abstracted &, const Topology &, Design &source);
void update_groups(const Abstracted &, Topology &);
void add_pseudonets(Abstracted &, const Topology &, const Parameters &);
} // namespace DreamPlace::native::macro
