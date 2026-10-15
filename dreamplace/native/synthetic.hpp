// BSD-3-Clause. See LICENSE.
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "dreamplace/native/hhds.hpp"

namespace DreamPlace::native::synthetic {
enum class Kind { mesh, torus, tree };
enum class Initialization { random, mesh_grid };
struct Spec {
  Kind kind = Kind::mesh;
  Initialization initialization = Initialization::random;
  bool side_io = false;  // One fixed input/output pair per mesh row, at left/right die edges.
  int width = 16, height = 16;
  int tree_depth = 3;  // 4^depth leaves, then pairwise reduction to one output
  unsigned seed = 42;
  double utilization = 0.6;
};

// Fixture-only physical/provenance attribute; never imposed on LiveHD graphs.
// Coordinates are microns here; physical files use 1000 DBU per micron.
struct CellState {
  double x = 0, y = 0;
  int group = 0, logical_x = 0, logical_y = 0, placed = 0;
};
struct CellTag {
  using value_type = CellState;
  using storage = hhds::flat_storage;
};
inline constexpr CellTag cell_state{};

struct Graph {
  Spec spec;
  std::string name;
  Layout layout;
  std::unique_ptr<hhds::GraphLibrary> library;
  std::shared_ptr<hhds::Graph> body;
};

Graph generate(const Spec& spec);
void save(const Graph& graph, const std::filesystem::path& directory);
Graph load(const std::filesystem::path& directory);
HhdsGeometry geometry(const Graph& graph);
std::string physical_pin_name(hhds::Pin pin);

// Full native flow: generate/checkpoint/reload -> export arrays -> place ->
// write HHDS positions -> checkpoint/reload -> LEF/DEF/PNG/CSV/JSON artifacts.
Result run(Graph graph, const std::filesystem::path& directory, const Options& options, bool animate = false);

struct AnimationFrame {
  PlacementStage stage;
  int iteration;
  std::vector<double> x, y;
};
void write_animation(const HhdsDesign& initial, const std::vector<AnimationFrame>& frames,
                     const std::filesystem::path& directory, const std::string& title, Initialization initialization);
void write_frame_png(const HhdsDesign& initial, const HhdsDesign& current, const std::filesystem::path& file,
                     const std::string& title, Initialization initialization, const AnimationFrame& frame,
                     size_t overlap_pairs);

// Synthetic format writers; these are not a general technology/LEF database.
void write_lef(const std::filesystem::path& file);
void write_def(const HhdsDesign& exported, const std::filesystem::path& file, const std::string& name);
void write_png(const HhdsDesign& initial, const HhdsDesign& placed, const std::filesystem::path& file,
               const std::string& title, const Result& result, Initialization initialization);
}  // namespace DreamPlace::native::synthetic

namespace hhds {
template <>
inline std::string attr_tag_name<DreamPlace::native::synthetic::CellTag>() {
  return "dreamplace.synthetic.cell.v1";
}
}  // namespace hhds
