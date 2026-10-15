// BSD-3-Clause. See LICENSE. End-to-end fixtures, not snapshots of optimizer coordinates.
#include "dreamplace/native/synthetic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace syn = DreamPlace::native::synthetic;
using namespace DreamPlace::native;
namespace {
void check(bool c, const char* message) {
  if (!c) throw std::runtime_error(message);
}
std::string read(const std::filesystem::path& path) {
  std::ifstream f(path, std::ios::binary);
  check(bool(f), "missing exported file");
  return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}
void topology(const syn::Graph& g, int expected_nodes, int expected_nets) {
  int nodes = 0, nets = 0, links = 0;
  std::unordered_map<hhds::Nid, std::vector<hhds::Nid>> adjacency;
  hhds::Nid start = 0;
  for (auto n : g.body->body().nodes()) {
    if (!start) start = n.get_debug_nid();
    ++nodes;
    int degree = 0, inputs = 0;
    std::multiset<std::pair<int, int>> neighbors;
    for (auto pin : n.inp_sorted_pins()) {
      (void)pin;
      ++inputs;
    }
    for (auto driver : n.out_sorted_pins()) {
      ++nets;
      ++degree;
      for (auto edge : driver.out_edges()) {
        ++links;
        auto sink = edge.sink.get_master_node();
        if (sink.is_output_node() && g.spec.side_io) {
          auto state = n.attr(syn::cell_state).get();
          check(state.logical_x == g.spec.width - 1 &&
                    syn::physical_pin_name(edge.sink) == "out_" + std::to_string(state.logical_y),
                "mesh output must connect to the matching rightmost row cell");
        }
        if (!sink.is_output_node()) {
          adjacency[n.get_debug_nid()].push_back(sink.get_debug_nid());
          auto state = sink.attr(syn::cell_state).get();
          neighbors.emplace(state.logical_x, state.logical_y);
        }
      }
    }
    if (g.spec.kind == syn::Kind::torus)
      check(degree == 4 && inputs == 4, "torus must have four neighbors in each direction");
    if (g.spec.kind == syn::Kind::mesh) {
      auto s = n.attr(syn::cell_state).get();
      int want =
          (s.logical_x > 0) + (s.logical_y > 0) + (s.logical_x + 1 < g.spec.width) + (s.logical_y + 1 < g.spec.height);
      check(degree == want + (g.spec.side_io && s.logical_x == g.spec.width - 1) &&
                inputs == want + (g.spec.side_io && s.logical_x == 0),
            "mesh boundary/interior/IO degree mismatch");
    }
    if (g.spec.kind != syn::Kind::tree) {
      auto state = n.attr(syn::cell_state).get();
      std::multiset<std::pair<int, int>> expected;
      for (auto [dx, dy] : {std::pair{1, 0}, {-1, 0}, {0, 1}, {0, -1}}) {
        int x = state.logical_x + dx, y = state.logical_y + dy;
        if (g.spec.kind == syn::Kind::torus) {
          x = (x + g.spec.width) % g.spec.width;
          y = (y + g.spec.height) % g.spec.height;
        }
        if (x >= 0 && x < g.spec.width && y >= 0 && y < g.spec.height) expected.emplace(x, y);
      }
      check(neighbors == expected, "grid edges or torus wraparound do not match logical neighbors");
    }
    if (g.spec.kind == syn::Kind::tree) check(degree <= 4 && inputs <= 2, "invalid expansion/reduction degree");
  }
  for (auto pin : g.body->get_input_node().out_sorted_pins()) {
    ++nets;
    int fanout = 0;
    for (auto edge : pin.out_edges()) {
      ++fanout;
      ++links;
      start = edge.sink.get_master_node().get_debug_nid();
      if (g.spec.side_io) {
        auto state = edge.sink.get_master_node().attr(syn::cell_state).get();
        check(state.logical_x == 0 && syn::physical_pin_name(pin) == "in_" + std::to_string(state.logical_y),
              "mesh input must connect to the matching leftmost row cell");
      }
    }
    if (g.spec.side_io) check(fanout == 1, "mesh input must have one connected cell");
  }
  check(nodes == expected_nodes && nets == expected_nets && links == expected_nets, "fixture node/net/link count");
  std::set<hhds::Nid> visited;
  std::queue<hhds::Nid> todo;
  todo.push(start);
  while (!todo.empty()) {
    auto n = todo.front();
    todo.pop();
    if (!visited.insert(n).second) continue;
    for (auto next : adjacency[n]) todo.push(next);
  }
  check(visited.size() == static_cast<size_t>(expected_nodes), "fixture has unreachable/disconnected cells");
  if (g.spec.kind == syn::Kind::tree) {
    check(g.body->get_io()->get_input_pin_decls().size() == 1 && g.body->get_io()->get_output_pin_decls().size() == 1,
          "tree needs exactly one input and one output");
    int leaves = 1;
    for (int i = 0; i < g.spec.tree_depth; ++i) leaves *= 4;
    int expand = 0, reduce = 0;
    for (auto n : g.body->body().nodes()) {
      if (n.get_name().starts_with("expand_"))
        ++expand;
      else if (n.get_name().starts_with("reduce_"))
        ++reduce;
    }
    check(expand == (4 * leaves - 1) / 3 && reduce == leaves - 1, "4-way expansion / 2-way reduction counts");
  }
}
void grid_seed(const syn::Graph& g) {
  auto e = from_hhds(*g.body, g.layout, syn::geometry(g));
  check(e.design.num_movable == static_cast<int>(e.cells.size()) &&
            e.terminals.size() == (g.spec.side_io ? 2 * g.spec.height : 0),
        "mesh grid seeds must remain movable");
  std::map<std::pair<int, int>, syn::CellState> states;
  for (auto node : e.cells) {
    auto s = node.attr(syn::cell_state).get();
    states.emplace(std::pair{s.logical_x, s.logical_y}, s);
  }
  for (const auto& [logical, s] : states) {
    auto [x, y] = logical;
    check(s.placed == 0 && s.x >= g.layout.xl && s.x + 2 <= g.layout.xh + 1e-10 && s.y >= g.layout.yl &&
              s.y + 1 <= g.layout.yh + 1e-10,
          "invalid mesh grid seed bounds/status");
    check(s.x == states.at({x, 0}).x && s.y == states.at({0, y}).y,
          "mesh seed rows/columns must align with logical coordinates");
    if (x) check(s.x > states.at({x - 1, y}).x, "mesh seed must progress left to right");
    if (y) check(s.y > states.at({x, y - 1}).y, "mesh seed must progress bottom to top");
  }
}
void side_anchors(const syn::Graph& g) {
  auto e = from_hhds(*g.body, g.layout, syn::geometry(g));
  check(e.terminals.size() == 2 * g.spec.height && g.body->get_io()->get_input_pin_decls().size() == g.spec.height &&
            g.body->get_io()->get_output_pin_decls().size() == g.spec.height,
        "mesh needs exactly one input/output pair per row");
  std::set<std::pair<bool, int>> seen;
  for (size_t i = 0; i < e.terminals.size(); ++i) {
    auto pin = e.terminals[i];
    bool input = pin.get_master_node().is_input_node();
    auto name = syn::physical_pin_name(pin);
    check(name.starts_with(input ? "in_" : "out_"), "mesh IO names must retain row identity");
    int row = std::stoi(name.substr(input ? 3 : 4));
    check(row >= 0 && row < g.spec.height && seen.emplace(input, row).second, "missing/duplicate mesh IO row");
    int node = e.cells.size() + i;
    double y = g.layout.yl + 0.5 + row * (g.layout.yh - g.layout.yl - 1) / (g.spec.height - 1);
    check(node >= e.design.num_movable && e.design.node_size_x[node] == 0 && e.design.node_size_y[node] == 0 &&
              e.design.node_x[node] == (input ? 0 : g.layout.xh + 2) && std::abs(e.design.node_y[node] - y) < 1e-12 &&
              y > 0 && y < g.layout.yh + 2,
          "mesh IO must be fixed on left/right edges, row aligned and away from top/bottom");
  }
}
void physical_files(const HhdsDesign& e, const std::filesystem::path& dir, const syn::Graph& graph) {
  // Read the emitted subset back independently of the writers, and compare
  // each instance and each net endpoint against the actual HHDS export.
  auto contents = read(dir / "placed.def");
  check(contents.find("UNITS DISTANCE MICRONS 1000 ;") != std::string::npos, "DEF units");
  check(contents.ends_with("END NETS\nEND DESIGN\n"), "truncated DEF");
  std::istringstream f(contents);
  std::string line, net, terminal;
  enum { other, components, pins, nets } section = other;
  std::map<std::string, std::pair<long long, long long>> positions;
  std::map<std::string, std::pair<long long, long long>> terminal_positions;
  std::set<std::string> endpoints;
  std::regex component(R"(- (\S+) SYNTH4 \+ PLACED \( ([0-9]+) ([0-9]+) \) N ;)");
  std::regex connection(R"(\s*\( (\S+) (\S+) \)\s*;?)");
  std::regex fixed_pin(R"(\s*\+ FIXED \( ([0-9]+) ([0-9]+) \) N ;)");
  int count_components = -1, count_nets = -1, rows = 0;
  for (; std::getline(f, line);) {
    if (line.starts_with("ROW ")) ++rows;
    if (line.starts_with("COMPONENTS ")) {
      section = components;
      count_components = std::stoi(line.substr(11));
      continue;
    }
    if (line == "END COMPONENTS") {
      section = other;
      continue;
    }
    if (line.starts_with("PINS ")) {
      section = pins;
      continue;
    }
    if (line == "END PINS") {
      section = other;
      continue;
    }
    if (line.starts_with("NETS ")) {
      section = nets;
      count_nets = std::stoi(line.substr(5));
      continue;
    }
    if (line == "END NETS") {
      section = other;
      continue;
    }
    std::smatch match;
    if (section == components) {
      check(std::regex_match(line, match, component), "invalid DEF component record");
      check(positions.emplace(match[1].str(), std::pair{std::stoll(match[2].str()), std::stoll(match[3].str())}).second,
            "duplicate DEF instance");
    }
    if (section == pins) {
      if (line.starts_with("- "))
        terminal = line.substr(2, line.find(' ', 2) - 2);
      else if (std::regex_match(line, match, fixed_pin))
        check(terminal_positions.emplace(terminal, std::pair{std::stoll(match[1].str()), std::stoll(match[2].str())})
                  .second,
              "duplicate DEF fixed IO pin");
    }
    if (section == nets) {
      if (line.starts_with("- "))
        net = line.substr(2);
      else {
        check(std::regex_match(line, match, connection), "invalid DEF net endpoint");
        check(endpoints.insert(net + ":" + match[1].str() + ":" + match[2].str()).second, "duplicate DEF endpoint");
      }
    }
  }
  const auto& d = e.design;
  check(rows == std::llround((d.layout.yh - d.layout.yl) / d.layout.row_height), "DEF row count");
  check(count_components == static_cast<int>(e.cells.size()) && count_nets == static_cast<int>(e.net_drivers.size()),
        "DEF counts");
  check(positions.size() == e.cells.size() && endpoints.size() == e.pins.size(), "DEF omitted cells/pins");
  check(terminal_positions.size() == e.terminals.size(), "DEF omitted fixed IO positions");
  for (auto pin : e.terminals) {
    bool input = pin.get_master_node().is_input_node();
    double y = (d.layout.yl + d.layout.yh) / 2;
    if (graph.spec.side_io) {
      int row = std::stoi(syn::physical_pin_name(pin).substr(input ? 3 : 4));
      y = d.layout.yl + 0.5 + row * (d.layout.yh - d.layout.yl - 1) / (graph.spec.height - 1);
    }
    check(terminal_positions.at(syn::physical_pin_name(pin)) ==
              std::pair{input ? 0LL : std::llround((d.layout.xh + 2) * 1000), std::llround(y * 1000)},
          "DEF input/output pins must match their fixed left/right die-edge coordinates");
  }
  for (size_t i = 0; i < e.cells.size(); ++i)
    check(positions.at(std::string(e.cells[i].get_name())) ==
              std::pair{std::llround(d.node_x[i] * 1000), std::llround(d.node_y[i] * 1000)},
          "DEF position differs from placed HHDS");
  for (size_t i = 0; i < e.pins.size(); ++i) {
    int node = d.pin2node[i];
    auto instance = node < static_cast<int>(e.cells.size()) ? std::string(e.cells[node].get_name()) : "PIN";
    auto key = "n" + std::to_string(d.pin2net[i]) + ":" + instance + ":" + syn::physical_pin_name(e.pins[i]);
    check(endpoints.contains(key), "DEF connectivity differs from HHDS");
  }
  auto lef = read(dir / "cells.lef");
  check(lef.find("DATABASE MICRONS 1000 ;") != std::string::npos &&
            lef.find("SIZE 2.000 BY 1.000 ;") != std::string::npos,
        "LEF geometry/units");
  for (const auto& p : {"I0", "I1", "I2", "I3", "O0", "O1", "O2", "O3"})
    check(lef.find(std::string("  PIN ") + p + "\n") != std::string::npos, "missing LEF physical pin");
  auto png = read(dir / "placement.png");
  check(png.size() > 1600 * 900 * 3 && png.substr(0, 8) == std::string("\x89PNG\r\n\x1a\n", 8),
        "invalid PNG signature/data");
  auto be = [&](int off) {
    unsigned v = 0;
    for (int j = 0; j < 4; ++j) v = (v << 8) | static_cast<unsigned char>(png[off + j]);
    return v;
  };
  check(png.substr(12, 4) == "IHDR" && be(16) == 1600 && be(20) == 900, "PNG dimensions");
}
}  // namespace
int main(int argc, char** argv) {
  try {
    check(argc == 2, "expected mesh, mesh_grid, mesh_grid_io, torus or tree");
    std::string kind = argv[1];
    syn::Spec spec;
    int nodes, nets;
    if (kind == "mesh" || kind == "mesh_grid" || kind == "mesh_grid_io") {
      spec.kind = syn::Kind::mesh;
      if (kind != "mesh") spec.initialization = syn::Initialization::mesh_grid;
      spec.side_io = kind == "mesh_grid_io";
      nodes = 256;
      nets = spec.side_io ? 992 : 960;
    } else if (kind == "torus") {
      spec.kind = syn::Kind::torus;
      spec.width = 12;
      spec.height = 8;
      nodes = 96;
      nets = 384;
    } else if (kind == "tree") {
      spec.kind = syn::Kind::tree;
      nodes = 148;
      nets = 212;
    } else
      throw std::invalid_argument("unknown test topology");
    auto graph = syn::generate(spec);
    topology(graph, nodes, nets);
    if (spec.initialization == syn::Initialization::mesh_grid) {
      grid_seed(graph);
      auto small = spec;
      for (auto [w, h] : {std::pair{7, 5}, {2, 2}}) {
        small.width = w;
        small.height = h;
        auto rectangular = syn::generate(small);
        grid_seed(rectangular);
        if (spec.side_io) side_anchors(rectangular);
      }
    }
    if (spec.side_io) {
      side_anchors(graph);
      auto control_spec = spec;
      control_spec.side_io = false;
      auto control = syn::generate(control_spec);
      auto a = from_hhds(*graph.body, graph.layout, syn::geometry(graph));
      auto b = from_hhds(*control.body, control.layout, syn::geometry(control));
      check(a.cells.size() == b.cells.size() && graph.layout.xl == control.layout.xl &&
                graph.layout.xh == control.layout.xh && graph.layout.yl == control.layout.yl &&
                graph.layout.yh == control.layout.yh,
            "side IO must not change core or cell count");
      for (size_t i = 0; i < a.cells.size(); ++i)
        check(a.cells[i].get_name() == b.cells[i].get_name() && a.design.node_x[i] == b.design.node_x[i] &&
                  a.design.node_y[i] == b.design.node_y[i],
              "side IO must not change mesh starting coordinates");
    }
    if (spec.kind == syn::Kind::tree) {
      auto e = from_hhds(*graph.body, graph.layout, syn::geometry(graph));
      check(e.terminals.size() == 2, "tree must export two fixed IO anchors");
      for (size_t i = 0; i < e.terminals.size(); ++i) {
        int node = e.cells.size() + i;
        bool input = e.terminals[i].get_master_node().is_input_node();
        check(node >= e.design.num_movable && e.design.node_size_x[node] == 0 && e.design.node_size_y[node] == 0 &&
                  e.design.node_x[node] == (input ? 0 : graph.layout.xh + 2) &&
                  e.design.node_y[node] == (graph.layout.yl + graph.layout.yh) / 2,
              "tree input/output must constrain placement at fixed left/right die edges");
      }
    }
    auto dir = std::filesystem::path(std::getenv("TEST_TMPDIR") ? std::getenv("TEST_TMPDIR") : "/tmp") / graph.name;
    Options options;
    options.max_iterations = 1500;
    options.bins_x = options.bins_y = 16;
    bool animate = spec.initialization == syn::Initialization::mesh_grid;
    options.observer_interval = 1000;  // Sparse end-to-end export; core test covers regular sampling.
    auto result = syn::run(std::move(graph), dir, options, animate);
    check(result.legal && result.converged && result.overflow <= options.stop_overflow && std::isfinite(result.hpwl),
          "placement must converge and be legal");
    // A mesh seed may begin with overlaps; making it legal can increase HPWL.
    // Keep the improvement requirement on all original random-start fixtures.
    if (spec.initialization == syn::Initialization::random)
      check(result.hpwl < result.initial_hpwl, "placement must improve random seed HPWL");
    auto before = syn::load(dir / "initial.hhds"), after = syn::load(dir / "placed.hhds");
    topology(before, nodes, nets);
    topology(after, nodes, nets);
    check(before.spec.initialization == spec.initialization && after.spec.initialization == spec.initialization &&
              before.spec.side_io == spec.side_io && after.spec.side_io == spec.side_io,
          "checkpoint lost initialization provenance");
    if (spec.initialization == syn::Initialization::mesh_grid) grid_seed(before);
    if (spec.side_io) {
      side_anchors(before);
      side_anchors(after);
    }
    auto e = from_hhds(*after.body, after.layout, syn::geometry(after));
    for (auto n : e.cells) check(n.attr(syn::cell_state).get().placed == 1, "checkpoint lost placed status");
    check(is_legal(e.design), "persisted placement is illegal");
    physical_files(e, dir, after);
    if (animate) {
      auto html = read(dir / "animation.html");
      check(html.find("id=\"seek\"") != std::string::npos && html.find("id=\"play\"") != std::string::npos &&
                html.find("detailed placement") != std::string::npos,
            "animation needs playback, seeking and phase labels");
      auto csv = read(dir / "animation.csv");
      check(csv.starts_with("frame,stage,gp_iteration,hpwl_um,legal,duration_seconds,overlap_pairs,pairwise_overlap_"
                            "area_um2\n0,initial,0,") &&
                csv.find(",global placement,") != std::string::npos &&
                csv.find(",legalization,") != std::string::npos && csv.find(",finished,") != std::string::npos,
            "animation needs real iteration/phase metrics");
      check(csv.ends_with(",2,0,0\n"), "finished animation must report zero physical overlap");
      // Independently check the default seed's 16 rows of 15 horizontal overlaps.
      auto data_start = csv.find('\n') + 1;
      auto data_end = csv.find('\n', data_start);
      auto first_row = csv.substr(data_start, data_end - data_start);
      auto area_comma = first_row.rfind(',');
      auto pairs_comma = first_row.rfind(',', area_comma - 1);
      check(std::stoi(first_row.substr(pairs_comma + 1, area_comma - pairs_comma - 1)) == 240,
            "initial grid overlap count must match its 16 rows of 15 adjacent pairs");
      check(std::abs(std::stod(first_row.substr(area_comma + 1)) - 48) < 1e-8,
            "initial pairwise overlap area must match the 0.2um by 1um intersections");
      auto concat = read(dir / "animation.ffconcat");
      check(concat.starts_with("ffconcat version 1.0\n"), "video export manifest");
      size_t frame_count = 0;
      for (const auto& file : std::filesystem::directory_iterator(dir / "animation_frames")) {
        auto png = read(file.path());
        check(png.substr(0, 8) == std::string("\x89PNG\r\n\x1a\n", 8) &&
                  png.substr(16, 8) == std::string("\0\0\3\x20\0\0\1\xc2", 8),
              "animation frame must be 800x450 PNG");
        ++frame_count;
      }
      check(frame_count >= 5 && frame_count <= 6, "sparse animation needs all phases plus optional interval frame");
    }
    std::cout << kind << ": " << nodes << " cells, " << nets << " nets; HPWL " << result.initial_hpwl << " -> "
              << result.hpwl << "; legal=" << result.legal << ", GP converged=" << result.converged << '\n';
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
}
