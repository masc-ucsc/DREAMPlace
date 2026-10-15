// BSD-3-Clause. See LICENSE.
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>

#include "dreamplace/native/synthetic.hpp"

namespace DreamPlace::native::synthetic {
namespace {
constexpr int dbu = 1000;
long long units(double x) { return std::llround(x * dbu); }
std::ofstream output(const std::filesystem::path& path) {
  std::ofstream f(path);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  return f;
}
void mapping_check(const HhdsDesign& e) {
  e.design.validate();
  if (e.pins.size() != e.design.pin2node.size() || e.cells.size() + e.terminals.size() != e.design.node_x.size())
    throw std::invalid_argument("physical exporter needs complete HHDS node/pin identities");
  const auto& d = e.design;
  if (d.layout.site_width != 1 || d.layout.row_height != 1 || d.layout.xl < 2 || d.layout.yl < 2)
    throw std::invalid_argument("synthetic physical export requires its 1um site/row and die margin");
  for (size_t i = 0; i < e.cells.size(); ++i)
    if (d.node_size_x[i] != 2 || d.node_size_y[i] != 1)
      throw std::invalid_argument("synthetic LEF describes only 2um by 1um SYNTH4 cells");
  for (auto pin : e.pins) {
    auto node = pin.get_master_node();
    if (!node.is_input_node() && !node.is_output_node() && pin.get_port_id() >= 4)
      throw std::invalid_argument("synthetic physical export supports only I0..I3/O0..O3");
  }
}
}  // namespace

void write_lef(const std::filesystem::path& path) {
  auto f = output(path);
  f << "# Synthetic visualization/test technology, not a fabrication library.\n"
    << "VERSION 5.8 ;\nBUSBITCHARS \"[]\" ;\nDIVIDERCHAR \"/\" ;\n"
    << "UNITS\n  DATABASE MICRONS 1000 ;\nEND UNITS\nMANUFACTURINGGRID 0.001 ;\n";
  for (int layer = 1; layer <= 2; ++layer)
    f << "LAYER metal" << layer << "\n  TYPE ROUTING ;\n  DIRECTION " << (layer == 1 ? "HORIZONTAL" : "VERTICAL")
      << " ;\n  PITCH 0.160 ;\n  WIDTH 0.080 ;\n  SPACING 0.080 ;\nEND metal" << layer << '\n';
  f << "SITE fake_site\n  CLASS CORE ;\n  SYMMETRY Y ;\n  SIZE 1.000 BY 1.000 ;\nEND fake_site\n"
    << "MACRO SYNTH4\n  CLASS CORE ;\n  ORIGIN 0 0 ;\n  SIZE 2.000 BY 1.000 ;\n  SYMMETRY X Y ;\n  SITE fake_site ;\n";
  f << std::fixed << std::setprecision(3);
  for (int driver = 0; driver < 2; ++driver)
    for (int pin = 0; pin < 4; ++pin) {
      std::string name = (driver ? "O" : "I") + std::to_string(pin);
      double x = driver ? 1.9 : 0.1, y = 0.125 + 0.25 * pin;
      f << "  PIN " << name << "\n    DIRECTION " << (driver ? "OUTPUT" : "INPUT") << " ;\n    USE SIGNAL ;\n"
        << "    PORT\n      LAYER metal1 ;\n      RECT " << x - 0.04 << ' ' << y - 0.04 << ' ' << x + 0.04 << ' '
        << y + 0.04 << " ;\n    END\n  END " << name << '\n';
    }
  f << "END SYNTH4\nEND LIBRARY\n";
  f.close();
}

void write_def(const HhdsDesign& e, const std::filesystem::path& path, const std::string& name) {
  mapping_check(e);
  if (!is_legal(e.design)) throw std::invalid_argument("refusing to export an illegal placed DEF");
  const auto& d = e.design;
  const auto& l = d.layout;
  auto f = output(path);
  f << "VERSION 5.8 ;\nDIVIDERCHAR \"/\" ;\nBUSBITCHARS \"[]\" ;\nDESIGN " << name << " ;\n"
    << "UNITS DISTANCE MICRONS 1000 ;\nDIEAREA ( 0 0 ) ( " << units(l.xh + 2) << ' ' << units(l.yh + 2) << " ) ;\n";
  int rows = std::llround((l.yh - l.yl) / l.row_height), sites = std::llround((l.xh - l.xl) / l.site_width);
  for (int r = 0; r < rows; ++r)
    f << "ROW row_" << r << " fake_site " << units(l.xl) << ' ' << units(l.yl + r * l.row_height) << " N DO " << sites
      << " BY 1 STEP " << units(l.site_width) << " 0 ;\n";
  f << "COMPONENTS " << e.cells.size() << " ;\n";
  for (size_t i = 0; i < e.cells.size(); ++i)
    f << "- " << e.cells[i].get_name() << " SYNTH4 + " << (i < static_cast<size_t>(d.num_movable) ? "PLACED" : "FIXED")
      << " ( " << units(d.node_x[i]) << ' ' << units(d.node_y[i]) << " ) N ;\n";
  f << "END COMPONENTS\nPINS " << e.terminals.size() << " ;\n";
  for (auto terminal : e.terminals) {
    int p = -1;
    for (size_t i = 0; i < e.pins.size(); ++i)
      if (e.pins[i].get_class_index() == terminal.get_class_index()) p = i;
    if (p < 0) throw std::invalid_argument("terminal missing from exported pins");
    int n = d.pin2node[p];
    bool input = terminal.get_master_node().is_input_node();
    f << "- " << physical_pin_name(terminal) << " + NET n" << d.pin2net[p] << " + DIRECTION "
      << (input ? "INPUT" : "OUTPUT") << " + USE SIGNAL\n  + LAYER metal1 ( " << (input ? 0 : -80) << " -40 ) ( "
      << (input ? 80 : 0) << " 40 )\n  + FIXED ( " << units(d.node_x[n]) << ' ' << units(d.node_y[n]) << " ) N ;\n";
  }
  f << "END PINS\nNETS " << d.net_weights.size() << " ;\n";
  for (size_t n = 0; n < d.net_weights.size(); ++n) {
    f << "- n" << n;
    for (int j = d.flat_net2pin_start[n]; j < d.flat_net2pin_start[n + 1]; ++j) {
      int p = d.flat_net2pin[j], node = d.pin2node[p];
      f << "\n  ( " << (node < static_cast<int>(e.cells.size()) ? std::string(e.cells[node].get_name()) : "PIN") << ' '
        << physical_pin_name(e.pins[p]) << " )";
    }
    f << " ;\n";
  }
  f << "END NETS\nEND DESIGN\n";
  f.close();
}

Result run(Graph graph, const std::filesystem::path& directory, const Options& options, bool animate) {
  std::filesystem::create_directories(directory);
  save(graph, directory / "initial.hhds");
  auto restored = load(directory / "initial.hhds");
  auto initial = from_hhds(*restored.body, restored.layout, geometry(restored));
  auto placed = initial;
  std::vector<AnimationFrame> frames;
  Options traced = options;
  if (animate) {
    traced.observer = [&](const Design& d, PlacementStage stage, int iteration) {
      if (options.observer) options.observer(d, stage, iteration);
      frames.push_back({stage, iteration, d.node_x, d.node_y});
    };
  }
  auto begin = std::chrono::steady_clock::now();
  auto result = place(placed.design, traced);
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  if (!result.legal) throw std::runtime_error("synthetic flow requires a legal final placement");
  apply_positions(placed, [](hhds::Node node, double x, double y) {
    auto state = node.attr(cell_state).get();
    state.x = x;
    state.y = y;
    state.placed = 1;
    node.attr(cell_state).set(state);
  });
  save(restored, directory / "placed.hhds");
  // Verify the persisted result, not just the live array modified by place().
  auto disk = load(directory / "placed.hhds");
  auto exported = from_hhds(*disk.body, disk.layout, geometry(disk));
  if (exported.design.node_x != placed.design.node_x || exported.design.node_y != placed.design.node_y ||
      exported.design.pin2node != placed.design.pin2node || !is_legal(exported.design))
    throw std::runtime_error("placed HHDS checkpoint did not round-trip");
  write_lef(directory / "cells.lef");
  write_def(exported, directory / "placed.def", restored.name);
  write_png(initial, exported, directory / "placement.png", restored.name, result, restored.spec.initialization);
  if (animate) write_animation(initial, frames, directory, restored.name, restored.spec.initialization);
  auto csv = output(directory / "positions.csv");
  csv << "instance,group,initial_x_um,initial_y_um,placed_x_um,placed_y_um\n" << std::setprecision(17);
  for (size_t i = 0; i < exported.cells.size(); ++i)
    csv << exported.cells[i].get_name() << ',' << exported.cells[i].attr(cell_state).get().group << ','
        << initial.design.node_x[i] << ',' << initial.design.node_y[i] << ',' << exported.design.node_x[i] << ','
        << exported.design.node_y[i] << '\n';
  csv.close();
  auto internal = exported.design;
  for (size_t p = 0; p < internal.pin2node.size(); ++p)
    if (internal.pin2node[p] >= static_cast<int>(exported.cells.size())) internal.net_mask[internal.pin2net[p]] = 0;
  double internal_hpwl = hpwl(internal);
  auto json = output(directory / "qor.json");
  json << std::setprecision(12) << "{\n  \"design\": \"" << restored.name
       << "\",\n  \"cells\": " << exported.cells.size() << ",\n  \"nets\": " << exported.net_drivers.size()
       << ",\n  \"terminals\": " << exported.terminals.size() << ",\n  \"seed\": " << restored.spec.seed
       << ",\n  \"initialization\": \""
       << (restored.spec.initialization == Initialization::mesh_grid ? "mesh-grid" : "random") << '"'
       << ",\n  \"initial_legal\": " << (is_legal(initial.design) ? "true" : "false")
       << ",\n  \"side_io\": " << (restored.spec.side_io ? "true" : "false")
       << ",\n  \"requested_utilization\": " << restored.spec.utilization
       << ",\n  \"core_width_um\": " << restored.layout.xh - restored.layout.xl
       << ",\n  \"core_height_um\": " << restored.layout.yh - restored.layout.yl
       << ",\n  \"threads\": " << options.num_threads << ",\n  \"bins_x\": " << options.bins_x
       << ",\n  \"bins_y\": " << options.bins_y << ",\n  \"max_iterations\": " << options.max_iterations
       << ",\n  \"target_density\": " << options.target_density << ",\n  \"stop_overflow\": " << options.stop_overflow
       << ",\n  \"initial_density_weight\": " << options.initial_density_weight
       << ",\n  \"max_step\": " << options.max_step << ",\n  \"edge_weight\": 1,\n  \"timing_model\": \"none\",\n"
       << "  \"coordinate_units\": \"microns\",\n  \"def_dbu_per_micron\": 1000,\n"
       << "  \"initial_hpwl\": " << result.initial_hpwl << ",\n  \"hpwl\": " << result.hpwl
       << ",\n  \"internal_hpwl\": " << internal_hpwl << ",\n  \"io_hpwl\": " << result.hpwl - internal_hpwl
       << ",\n  \"gp_overflow\": " << result.overflow
       << ",\n  \"gp_converged\": " << (result.converged ? "true" : "false")
       << ",\n  \"legal\": " << (result.legal ? "true" : "false") << ",\n  \"iterations\": " << result.iterations
       << ",\n  \"placement_seconds\": " << seconds << "\n}\n";
  json.close();
  return result;
}
}  // namespace DreamPlace::native::synthetic
