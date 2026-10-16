// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/hierarchy.hpp"
#include "dreamplace/native/macro/seed.hpp"
#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
namespace DreamPlace::native::macro {
Hierarchy read_hierarchy(hhds::Graph &root, const HierarchyCallbacks &cb) {
  if (!cb.instance_name || !cb.kind || !cb.opaque || !cb.terminal_name)
    throw std::invalid_argument("missing hierarchy callbacks");
  Hierarchy out;
  struct Context {
    hhds::Graph *graph;
    std::string prefix;
  };
  struct Leaf {
    hhds::Node node;
    int context, index;
  };
  std::vector<Context> contexts;
  std::vector<Leaf> leaves;
  std::map<std::pair<int, hhds::Pid>, int> pins;
  std::vector<int> parent;
  auto pin = [&](int context, hhds::Pin p) {
    auto [it, added] =
        pins.emplace(std::make_pair(context, p.get_debug_pid()), parent.size());
    if (added)
      parent.push_back(it->second);
    return it->second;
  };
  auto find = [&](int i) {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };
  auto join = [&](int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b)
      parent[b] = a;
  };
  std::set<hhds::Gid> stack;
  std::set<std::string> stable_names;
  std::function<int(hhds::Graph &, std::string, bool)> visit =
      [&](hhds::Graph &graph, std::string prefix, bool parent_stable) {
        if (!stack.insert(graph.get_gid()).second)
          throw std::invalid_argument("recursive H2D hierarchy");
        int context = contexts.size();
        contexts.push_back({&graph, prefix});
        for (auto p : graph.get_input_node().out_sorted_pins())
          pin(context, p);
        for (auto p : graph.get_output_node().inp_sorted_pins())
          pin(context, p);
        for (auto node : graph.body().nodes()) {
          auto local = cb.instance_name(node);
          bool stable = parent_stable && !node.get_name().empty();
          if (local.empty())
            throw std::invalid_argument("empty H2D instance name");
          std::string name = prefix + local;
          if (!stable_names.insert(name).second)
            throw std::invalid_argument("duplicate H2D hierarchical name: " +
                                        name);
          auto io = node.get_subnode_io();
          std::string master = io ? std::string(io->get_name()) : std::string{};
          out.instances.push_back({name,
                                   prefix.empty()
                                       ? std::string{}
                                       : prefix.substr(0, prefix.size() - 1),
                                   master, stable});
          if (!stable && cb.kind(node) == FlowKind::macro)
            out.diagnostics.push_back("unstable generated instance path: " +
                                      name);
          for (auto p : node.inp_sorted_pins())
            pin(context, p);
          for (auto p : node.out_sorted_pins())
            pin(context, p);
          auto child = node.get_subnode_graph();
          if (child && !cb.opaque(node)) {
            int c = visit(*child, name + ".", stable);
            for (auto p : child->get_input_node().out_sorted_pins()) {
              auto outer = node.get_sink_pin(p.get_port_id());
              if (outer.is_valid())
                join(pin(context, outer), pin(c, p));
            }
            for (auto p : child->get_output_node().inp_sorted_pins()) {
              auto outer = node.get_driver_pin(p.get_port_id());
              if (outer.is_valid())
                join(pin(context, outer), pin(c, p));
            }
          } else {
            int index = out.flow.nodes.size();
            out.flow.nodes.push_back({name,
                                      prefix.empty()
                                          ? std::string{}
                                          : prefix.substr(0, prefix.size() - 1),
                                      cb.kind(node),
                                      {}});
            leaves.push_back({node, context, index});
          }
        }
        for (auto driver : graph.get_input_node().out_sorted_pins())
          for (auto edge : driver.out_edges())
            join(pin(context, driver), pin(context, edge.sink));
        for (auto node : graph.body().nodes())
          for (auto driver : node.out_sorted_pins())
            for (auto edge : driver.out_edges())
              join(pin(context, driver), pin(context, edge.sink));
        stack.erase(graph.get_gid());
        return context;
      };
  visit(root, "", true);
  std::map<int, std::vector<int>> drivers, sinks;
  for (const auto &leaf : leaves) {
    for (auto p : leaf.node.out_sorted_pins())
      drivers[find(pin(leaf.context, p))].push_back(leaf.index);
    for (auto p : leaf.node.inp_sorted_pins())
      sinks[find(pin(leaf.context, p))].push_back(leaf.index);
  }
  auto terminal = [&](hhds::Pin p, bool input) {
    int i = out.flow.nodes.size();
    out.flow.nodes.push_back(
        {"@io:" + cb.terminal_name(p), {}, FlowKind::io, {}});
    (input ? drivers : sinks)[find(pin(0, p))].push_back(i);
  };
  for (auto p : root.get_input_node().out_sorted_pins())
    terminal(p, true);
  for (auto p : root.get_output_node().inp_sorted_pins())
    terminal(p, false);
  for (const auto &[net, from] : drivers) {
    auto it = sinks.find(net);
    if (it != sinks.end())
      for (int source : from)
        for (int sink : it->second)
          if (source != sink)
            out.flow.nodes[source].sinks.push_back(sink);
  }
  for (auto &n : out.flow.nodes) {
    std::sort(n.sinks.begin(), n.sinks.end());
    n.sinks.erase(std::unique(n.sinks.begin(), n.sinks.end()), n.sinks.end());
  }
  return out;
}
Topology seed(const Design &d, std::span<const NodeInfo> info,
              const Hierarchy &h, const Parameters &p, const Topology &input) {
  std::vector<NodeInfo> named(info.begin(), info.end());
  std::map<std::string, const Hierarchy::Instance *> instances;
  for (const auto &i : h.instances)
    instances[i.name] = &i;
  for (auto &n : named)
    if (auto i = instances.find(n.name); i != instances.end()) {
      n.parent = i->second->parent;
      n.stable = i->second->stable;
    }
  auto flows = dataflow(h.flow, p.hops_max, p.fanout_cap);
  auto out = seed_topology(d, named, flows, p, input);
  out.diagnostics.insert(out.diagnostics.end(), h.diagnostics.begin(),
                         h.diagnostics.end());
  return out;
}
} // namespace DreamPlace::native::macro
