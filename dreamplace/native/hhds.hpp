// BSD-3-Clause. See LICENSE.
#pragma once
#include <functional>
#include <utility>

#include "dreamplace/native/placement.hpp"
#include "hhds/graph.hpp"

namespace DreamPlace::native {
enum class Status { movable, fixed, locked };
struct CellGeometry {
  double x = 0, y = 0, width = 0, height = 0;
  Status status = Status::movable;
};
struct NetProperties {
  double weight = 1;
  bool enabled = true;
};
struct HhdsGeometry {
  // Every non-boundary node must be a mapped physical cell. Leaf black-box
  // subnodes (GraphIO without a body) are accepted; hierarchy must be flattened.
  std::function<CellGeometry(hhds::Node)> cell;
  // Offsets from the oriented cell's lower-left corner, in layout units.
  std::function<std::pair<double, double>(hhds::Pin)> pin_offset;
  // Each connected top-level pin becomes its own fixed, zero-area terminal.
  std::function<std::pair<double, double>(hhds::Pin)> terminal_position;
  std::function<NetProperties(hhds::Pin)> net;  // optional; default weight=1
};

struct HhdsDesign {
  Design design;
  // Cell handles follow design node order; terminals are recorded separately.
  std::vector<hhds::Node> cells;
  std::vector<hhds::Pin> terminals;
  std::vector<hhds::Pin> net_drivers;
  // Original pin identities, indexed exactly like Design::pin2node. Exporters
  // can retain physical pin names without reconstructing traversal order.
  std::vector<hhds::Pin> pins;
  hhds::Graph* graph = nullptr;  // non-owning; graph/library must outlive this object
  std::uint64_t epoch = 0;
};

// In-memory adapter. Read persisted databases with GraphLibrary::load first.
// Preserves one hypernet per driver pin, including fanout and multiple outputs.
// Constant nets are excluded (no physical source); materialize tie cells first
// when physical constant wiring is required.
HhdsDesign from_hhds(hhds::Graph& graph, const Layout& layout, const HhdsGeometry& geometry);

// Only movable cells are written. Callback owns rounding to integer DBU and
// updating LiveHD's place_t/status; no LiveHD attributes are built into this fork.
// Rejects structural edits made since export. Caller keeps the graph alive and
// does not mutate topology from inside the writer.
void apply_positions(const HhdsDesign& exported, const std::function<void(hhds::Node, double, double)>& writer);
}  // namespace DreamPlace::native
