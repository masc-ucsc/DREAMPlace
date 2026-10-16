// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/topology.hpp"
namespace DreamPlace::native::macro {
enum class FlowKind { combinational, sequential, macro, io };
struct FlowNode {
  std::string name, parent;
  FlowKind kind = FlowKind::combinational;
  std::vector<int> sinks;
};
struct FlowGraph {
  std::vector<FlowNode> nodes;
};
// Register hops, not combinational gate depth. Each endpoint pair receives its
// strongest shortest-hop connection; repeated paths do not amplify cycles.
// Sequential sources/destinations share hierarchy-cluster endpoints.
std::vector<Flow> dataflow(const FlowGraph &, int hops_max = 4,
                           int fanout_cap = 4096);
} // namespace DreamPlace::native::macro
