// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/dataflow.hpp"
#include "hhds/graph.hpp"
#include <functional>
namespace DreamPlace::native::macro {
struct HierarchyCallbacks {
  std::function<std::string(hhds::Node)> instance_name;
  std::function<FlowKind(hhds::Node)> kind;
  std::function<bool(hhds::Node)> opaque;
  std::function<std::string(hhds::Pin)> terminal_name;
};
struct Hierarchy {
  struct Instance {
    std::string name, parent, master;
    bool stable = true;
  };
  std::vector<Instance> instances;
  FlowGraph flow;
  std::vector<std::string> diagnostics;
};
// Reads the hierarchical graph before flattening. Interfaces are joined by port
// identity; module inputs are never shorted directly to module outputs. No
// LiveHD attributes, benchmark inputs, or OpenROAD runtime dependency.
Hierarchy read_hierarchy(hhds::Graph &, const HierarchyCallbacks &);
Topology seed(const Design &, std::span<const NodeInfo>, const Hierarchy &,
              const Parameters &, const Topology &input = {});
} // namespace DreamPlace::native::macro
