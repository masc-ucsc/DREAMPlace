// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/hhds.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace DreamPlace::native {
HhdsDesign from_hhds(hhds::Graph& graph, const Layout& layout, const HhdsGeometry& geom) {
  if (!geom.cell || !geom.pin_offset || !geom.terminal_position)
    throw std::invalid_argument("HHDS adapter requires cell, pin-offset and terminal-position callbacks");
  HhdsDesign out;
  out.graph = &graph;
  out.design.layout = layout;
  auto& d = out.design;
  struct Cell {
    hhds::Node node;
    CellGeometry geom;
  };
  std::vector<Cell> cells;
  for (auto node : graph.body().nodes()) {
    auto io = node.get_subnode_io();
    if (io && io->has_graph()) throw std::invalid_argument("HHDS placement requires a flat mapped graph");
    cells.push_back({node, geom.cell(node)});
  }
  std::stable_sort(cells.begin(), cells.end(), [](const auto& a, const auto& b) {
    if ((a.geom.status == Status::movable) != (b.geom.status == Status::movable))
      return a.geom.status == Status::movable;
    return a.node.get_class_index().value < b.node.get_class_index().value;
  });
  std::unordered_map<hhds::Nid, int> node_index;
  for (const auto& c : cells) {
    if (c.geom.status == Status::movable) ++d.num_movable;
    node_index.emplace(c.node.get_debug_nid(), d.node_x.size());
    d.node_x.push_back(c.geom.x);
    d.node_y.push_back(c.geom.y);
    d.node_size_x.push_back(c.geom.width);
    d.node_size_y.push_back(c.geom.height);
    out.cells.push_back(c.node);
  }
  std::unordered_map<hhds::Pid, int> terminal_index;
  auto terminal = [&](hhds::Pin p) {
    auto [it, added] = terminal_index.emplace(p.get_debug_pid(), d.node_x.size());
    if (added) {
      auto [x, y] = geom.terminal_position(p);
      d.node_x.push_back(x);
      d.node_y.push_back(y);
      d.node_size_x.push_back(0);
      d.node_size_y.push_back(0);
      out.terminals.push_back(p);
    }
    return it->second;
  };
  auto pin = [&](hhds::Pin p) {
    if (p.get_graph() != &graph) throw std::invalid_argument("HHDS edge crosses unflattened hierarchy");
    auto node = p.get_master_node();
    int index;
    double x = 0, y = 0;
    if (node.is_input_node() || node.is_output_node())
      index = terminal(p);
    else {
      auto it = node_index.find(node.get_debug_nid());
      if (it == node_index.end()) throw std::invalid_argument("net endpoint is not a mapped physical cell");
      index = it->second;
      auto offset = geom.pin_offset(p);
      x = offset.first;
      y = offset.second;
    }
    if (d.pin2node.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
      throw std::invalid_argument("HHDS placement exceeds kernel pin indices");
    d.flat_net2pin.push_back(d.pin2node.size());
    out.pins.push_back(p);
    d.pin2node.push_back(index);
    d.pin_offset_x.push_back(x);
    d.pin_offset_y.push_back(y);
  };
  std::vector<hhds::Pin> drivers;
  for (auto node : out.cells)
    for (auto driver : node.out_sorted_pins()) drivers.push_back(driver);
  for (auto driver : graph.get_input_node().out_sorted_pins()) drivers.push_back(driver);
  std::sort(drivers.begin(), drivers.end(),
            [](auto a, auto b) { return a.get_class_index().value < b.get_class_index().value; });
  for (auto driver : drivers) {
    std::vector<hhds::Pin> sinks;
    for (auto edge : driver.out_edges()) sinks.push_back(edge.sink);
    if (sinks.empty()) continue;
    std::sort(sinks.begin(), sinks.end(),
              [](auto a, auto b) { return a.get_class_index().value < b.get_class_index().value; });
    pin(driver);
    for (auto sink : sinks) pin(sink);
    auto net = geom.net ? geom.net(driver) : NetProperties{};
    d.net_weights.push_back(net.weight);
    d.net_mask.push_back(net.enabled);
    d.flat_net2pin_start.push_back(d.flat_net2pin.size());
    out.net_drivers.push_back(driver);
  }
  d.build_maps();
  d.validate();
  out.epoch = graph.body_epoch();
  return out;
}

void apply_positions(const HhdsDesign& out, const std::function<void(hhds::Node, double, double)>& writer) {
  if (!out.graph || out.graph->body_epoch() != out.epoch)
    throw std::invalid_argument("HHDS graph changed after placement export; rebuild arrays");
  if (!writer) throw std::invalid_argument("missing HHDS position writer");
  out.design.validate();
  if (out.cells.size() + out.terminals.size() != out.design.node_x.size() ||
      out.pins.size() != out.design.pin2node.size() || out.cells.size() < static_cast<size_t>(out.design.num_movable))
    throw std::invalid_argument("HHDS export mapping does not match placement arrays");
  for (int i = 0; i < out.design.num_movable; ++i) {
    if (!out.cells[i].is_valid()) throw std::invalid_argument("HHDS cell deleted since export");
    writer(out.cells[i], out.design.node_x[i], out.design.node_y[i]);
    if (out.graph->body_epoch() != out.epoch)
      throw std::invalid_argument("HHDS position writer changed graph topology");
  }
}
}  // namespace DreamPlace::native
