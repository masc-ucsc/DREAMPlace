// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/synthetic.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <vector>

namespace DreamPlace::native::synthetic {
namespace {
void register_attributes() { hhds::register_attr_tag<CellTag>(hhds::attr_tag_name<CellTag>()); }
void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
void validate(const Spec& s) {
  check(s.kind == Kind::mesh || s.kind == Kind::torus || s.kind == Kind::tree, "unknown synthetic topology");
  check(s.initialization == Initialization::random ||
            (s.initialization == Initialization::mesh_grid && s.kind == Kind::mesh),
        "mesh-grid initialization requires a mesh");
  check(!s.side_io || (s.kind == Kind::mesh && s.initialization == Initialization::mesh_grid),
        "side IO requires a grid-start mesh");
  check(s.width >= 2 && s.width <= 64 && s.height >= 2 && s.height <= 64, "grid dimensions must be 2..64");
  check(s.tree_depth >= 1 && s.tree_depth <= 5, "tree depth must be 1..5");
  check(std::isfinite(s.utilization) && s.utilization >= 0.1 && s.utilization <= 0.8, "utilization must be 0.1..0.8");
}
}  // namespace

Graph generate(const Spec& spec) {
  validate(spec);
  register_attributes();
  Graph g;
  g.spec = spec;
  g.library = std::make_unique<hhds::GraphLibrary>();
  g.name = spec.kind == Kind::tree ? "tree_4to2_d" + std::to_string(spec.tree_depth)
                                   : (spec.kind == Kind::mesh ? "mesh_" : "torus_") + std::to_string(spec.width) + "x" +
                                         std::to_string(spec.height);
  if (spec.initialization == Initialization::mesh_grid) g.name += "_grid";
  if (spec.side_io) g.name += "_io";
  auto io = g.library->create_io(g.name);
  if (spec.kind == Kind::tree) {
    io->add_input("in", 1);
    io->add_output("out", 2);
  }
  if (spec.side_io)
    for (int y = 0; y < spec.height; ++y) {
      io->add_input("in_" + std::to_string(y), y + 1);
      io->add_output("out_" + std::to_string(y), spec.height + y + 1);
    }
  g.body = io->create_graph();
  std::vector<hhds::Node> nodes;
  std::vector<unsigned> inputs, outputs;
  auto node = [&](const std::string& name, int group, int lx, int ly) {
    auto n = g.body->create_node();
    n.set_name(name);
    n.attr(cell_state).set(CellState{0, 0, group, lx, ly, 0});
    nodes.push_back(n);
    inputs.push_back(0);
    outputs.push_back(0);
    return static_cast<int>(nodes.size()) - 1;
  };
  auto connect = [&](int a, int b) {
    check(outputs[a] < 4 && inputs[b] < 4, "synthetic cell exceeds four physical input/output pins");
    nodes[a].create_driver_pin(outputs[a]++).connect_sink(nodes[b].create_sink_pin(inputs[b]++));
  };
  if (spec.kind != Kind::tree) {
    for (int y = 0; y < spec.height; ++y)
      for (int x = 0; x < spec.width; ++x) node("r" + std::to_string(y) + "_c" + std::to_string(x), x, x, y);
    auto pair = [&](int a, int b) {
      connect(a, b);
      connect(b, a);
    };
    for (int y = 0; y < spec.height; ++y)
      for (int x = 0; x < spec.width; ++x) {
        int a = y * spec.width + x;
        if (x + 1 < spec.width || spec.kind == Kind::torus) pair(a, y * spec.width + (x + 1) % spec.width);
        if (y + 1 < spec.height || spec.kind == Kind::torus) pair(a, ((y + 1) % spec.height) * spec.width + x);
      }
    if (spec.side_io)
      for (int y = 0; y < spec.height; ++y) {
        int left = y * spec.width, right = left + spec.width - 1;
        check(inputs[left] < 4 && outputs[right] < 4, "mesh IO exceeds physical macro ports");
        g.body->get_input_pin("in_" + std::to_string(y)).connect_sink(nodes[left].create_sink_pin(inputs[left]++));
        nodes[right]
            .create_driver_pin(outputs[right]++)
            .connect_sink(g.body->get_output_pin("out_" + std::to_string(y)));
      }
  } else {
    int root = node("expand_0_0", 0, 0, 0);
    g.body->get_input_pin("in").connect_sink(nodes[root].create_sink_pin(inputs[root]++));
    std::vector<int> frontier{root};
    for (int depth = 1; depth <= spec.tree_depth; ++depth) {
      std::vector<int> next;
      for (int parent : frontier)
        for (int child = 0; child < 4; ++child) {
          int n =
              node("expand_" + std::to_string(depth) + "_" + std::to_string(next.size()), depth, depth, next.size());
          connect(parent, n);
          next.push_back(n);
        }
      frontier = std::move(next);
    }
    int stage = spec.tree_depth + 1;
    while (frontier.size() > 1) {
      std::vector<int> next;
      for (size_t i = 0; i < frontier.size(); i += 2) {
        int n = node("reduce_" + std::to_string(stage - spec.tree_depth - 1) + "_" + std::to_string(i / 2), stage,
                     stage, i / 2);
        connect(frontier[i], n);
        connect(frontier[i + 1], n);
        next.push_back(n);
      }
      frontier = std::move(next);
      ++stage;
    }
    auto out = frontier.front();
    nodes[out].create_driver_pin(outputs[out]++).connect_sink(g.body->get_output_pin("out"));
  }
  // Mesh variants use identical physical geometry; zero-area IO adds no cell area.
  double aspect = spec.kind == Kind::tree ? 1.6 : double(spec.width) / spec.height;
  double area = nodes.size() * 2.0 / spec.utilization;
  int rows = std::max(4, static_cast<int>(std::ceil(std::sqrt(area / aspect))));
  int sites = std::max(4, static_cast<int>(std::ceil(area / rows)));
  g.layout = {2, 2, double(sites + 2), double(rows + 2), 1, 1};
  std::mt19937 rng(spec.seed);
  std::uniform_real_distribution<double> ux(g.layout.xl, g.layout.xh - 2), uy(g.layout.yl, g.layout.yh - 1);
  for (auto n : nodes) {
    auto s = n.attr(cell_state).get();
    if (spec.initialization == Initialization::mesh_grid) {
      // Spread the logical grid over the core, allowing initial overlap when
      // its columns cannot fit 2um-wide cells. These are movable seed positions,
      // not fences, fixed cells, or persistent relative-placement constraints.
      s.x = g.layout.xl + s.logical_x * (g.layout.xh - g.layout.xl - 2) / (spec.width - 1);
      s.y = g.layout.yl + s.logical_y * (g.layout.yh - g.layout.yl - 1) / (spec.height - 1);
    } else {
      s.x = ux(rng);
      s.y = uy(rng);
    }
    n.attr(cell_state).set(s);
  }
  return g;
}

void save(const Graph& g, const std::filesystem::path& directory) {
  std::filesystem::create_directories(directory);
  g.library->save(directory.string());
  std::ofstream f(directory / "fixture.txt");
  f << std::setprecision(17) << "DREAMPLACE_SYNTHETIC_V3\n"
    << g.name << '\n'
    << static_cast<int>(g.spec.kind) << ' ' << g.spec.width << ' ' << g.spec.height << ' ' << g.spec.tree_depth << ' '
    << g.spec.seed << ' ' << g.spec.utilization << ' ' << static_cast<int>(g.spec.initialization) << ' '
    << g.spec.side_io << '\n'
    << g.layout.xl << ' ' << g.layout.yl << ' ' << g.layout.xh << ' ' << g.layout.yh << ' ' << g.layout.site_width
    << ' ' << g.layout.row_height << '\n';
  f.close();
  check(bool(f), "failed to write synthetic HHDS manifest");
}

Graph load(const std::filesystem::path& directory) {
  register_attributes();
  Graph g;
  std::ifstream f(directory / "fixture.txt");
  std::string magic;
  int kind;
  check(bool(f >> magic >> g.name >> kind >> g.spec.width >> g.spec.height >> g.spec.tree_depth >> g.spec.seed >>
             g.spec.utilization) &&
            (magic == "DREAMPLACE_SYNTHETIC_V1" || magic == "DREAMPLACE_SYNTHETIC_V2" ||
             magic == "DREAMPLACE_SYNTHETIC_V3"),
        "invalid synthetic HHDS manifest");
  if (magic != "DREAMPLACE_SYNTHETIC_V1") {
    int initialization = -1;
    check(bool(f >> initialization), "missing synthetic initialization");
    g.spec.initialization = static_cast<Initialization>(initialization);
  }
  if (magic == "DREAMPLACE_SYNTHETIC_V3") {
    int side_io = -1;
    check(bool(f >> side_io) && (side_io == 0 || side_io == 1), "invalid synthetic side IO flag");
    g.spec.side_io = side_io;
  }
  check(
      bool(f >> g.layout.xl >> g.layout.yl >> g.layout.xh >> g.layout.yh >> g.layout.site_width >> g.layout.row_height),
      "invalid synthetic HHDS layout");
  g.spec.kind = static_cast<Kind>(kind);
  validate(g.spec);
  g.library = std::make_unique<hhds::GraphLibrary>();
  g.library->load(directory.string());
  auto io = g.library->find_io(g.name);
  check(io && io->has_graph(), "synthetic HHDS top body is missing");
  g.body = io->get_graph();
  return g;
}

HhdsGeometry geometry(const Graph& g) {
  HhdsGeometry geom;
  geom.cell = [](hhds::Node node) {
    check(node.attr(cell_state).has(), "missing synthetic cell attribute");
    auto s = node.attr(cell_state).get();
    return CellGeometry{s.x, s.y, 2, 1, Status::movable};
  };
  geom.pin_offset = [](hhds::Pin pin) {
    check(pin.get_port_id() < 4, "synthetic pin port exceeds physical macro");
    return std::pair<double, double>{pin.is_driver() ? 1.9 : 0.1, 0.125 + 0.25 * pin.get_port_id()};
  };
  geom.terminal_position = [l = g.layout, side_io = g.spec.side_io, height = g.spec.height](hhds::Pin pin) {
    // Fixed at the die perimeter, outside the 2um margin around the core.
    bool input = pin.get_master_node().is_input_node();
    double y = (l.yl + l.yh) / 2;
    if (side_io) {
      int row = int(pin.get_port_id()) - 1 - (input ? 0 : height);
      check(row >= 0 && row < height, "invalid mesh IO port");
      y = l.yl + 0.5 + row * (l.yh - l.yl - 1) / (height - 1);
    }
    return std::pair<double, double>{input ? 0 : l.xh + 2, y};
  };
  return geom;
}

std::string physical_pin_name(hhds::Pin pin) {
  auto node = pin.get_master_node();
  if (node.is_input_node() || node.is_output_node()) {
    auto io = pin.get_graph()->get_io();
    const auto& declarations = node.is_input_node() ? io->get_input_pin_decls() : io->get_output_pin_decls();
    for (const auto& decl : declarations)
      if (decl.port_id == pin.get_port_id()) return decl.name;
    throw std::invalid_argument("synthetic top pin has no declared name");
  }
  return (pin.is_driver() ? "O" : "I") + std::to_string(pin.get_port_id());
}
}  // namespace DreamPlace::native::synthetic
