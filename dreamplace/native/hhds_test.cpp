// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/hhds.hpp"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

using namespace DreamPlace::native;
namespace {
void check(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
template <typename F>
void rejects(F f) {
  bool rejected = false;
  try {
    f();
  } catch (const std::exception&) {
    rejected = true;
  }
  check(rejected, "expected HHDS rejection");
}
void test() {
  hhds::GraphLibrary lib;
  auto io = lib.create_io("mapped_top");
  io->add_input("a", 1);
  io->add_input("b", 2);
  io->add_output("z", 3);
  auto g = io->create_graph();
  auto fixed = g->create_node(), a = g->create_node(), b = g->create_node();
  fixed.set_name("fixed");
  a.set_name("a");
  b.set_name("b");
  // Distinct outputs on the same physical cell must remain distinct hypernets.
  auto a0 = a.create_driver_pin(0), a1 = a.create_driver_pin(1);
  g->get_input_pin("a").connect_sink(a.create_sink_pin(1));
  g->get_input_pin("b").connect_sink(b.create_sink_pin(2));
  a0.connect_sink(b.create_sink_pin(3));
  a0.connect_sink(fixed.create_sink_pin(4));
  a1.connect_sink(b.create_sink_pin(5));
  b.create_driver_pin(0).connect_sink(g->get_output_pin("z"));
  Layout layout{0, 0, 40, 10, 1, 1};
  HhdsGeometry geom;
  geom.cell = [](auto n) {
    return CellGeometry{20, 5, 2, 1, n.get_name() == "fixed" ? Status::locked : Status::movable};
  };
  geom.pin_offset = [](auto p) { return std::pair<double, double>{0.25 * p.get_port_id(), 0.5}; };
  geom.terminal_position = [](auto p) {
    return std::pair<double, double>{p.get_master_node().is_input_node() ? 0 : 40, 5};
  };
  geom.net = [](auto p) { return NetProperties{p.get_port_id() == 1 ? 2.0 : 1.0, true}; };
  auto e = from_hhds(*g, layout, geom);
  check(e.design.num_movable == 2 && e.cells.size() == 3 && e.terminals.size() == 3, "cell/terminal ordering");
  check(e.cells[2] == fixed, "locked cell must follow movable cells");
  check(e.net_drivers.size() == 5 && e.design.pin2node.size() == 11, "driver hypernets/fanout lost");
  check(e.design.flat_net2pin_start.back() == 11, "CSR end");
  check(e.pins.size() == e.design.pin2node.size(), "missing original pin identities");
  for (size_t i = 0; i < e.pins.size(); ++i) {
    auto master = e.pins[i].get_master_node();
    int node = e.design.pin2node[i];
    if (node < static_cast<int>(e.cells.size()))
      check(master == e.cells[node], "pin identity disagrees with flat node map");
    else
      check(master.is_input_node() || master.is_output_node(), "terminal pin identity lost");
  }
  auto e2 = from_hhds(*g, layout, geom);
  check(e2.design.flat_net2pin == e.design.flat_net2pin && e2.design.pin2node == e.design.pin2node,
        "stable adapter order");
  Options opts;
  opts.max_iterations = 20;
  opts.bins_x = 8;
  opts.bins_y = 8;
  Placer placer(opts);
  check(placer.place(e.design).legal, "place HHDS graph through configured instance");
  std::unordered_map<hhds::Nid, std::pair<double, double>> written;
  apply_positions(e, [&](auto n, double x, double y) { written[n.get_debug_nid()] = {x, y}; });
  check(written.size() == 2 && !written.contains(fixed.get_debug_nid()), "writeback moved locked cell");
  // Round-trip an actual persisted HHDS library, not just a synthetic array.
  auto path =
      std::filesystem::path(std::getenv("TEST_TMPDIR") ? std::getenv("TEST_TMPDIR") : "/tmp") / "dreamplace_hhds_test";
  std::filesystem::remove_all(path);
  lib.save(path.string());
  hhds::GraphLibrary restored;
  restored.load(path.string());
  auto restored_io = restored.find_io("mapped_top");
  check(restored_io != nullptr, "persisted graph missing");
  auto loaded = from_hhds(*restored_io->get_graph(), layout, geom);
  check(loaded.design.pin2node == e2.design.pin2node && loaded.design.pin_offset_x == e2.design.pin_offset_x,
        "persisted adapter mismatch");
  check(placer.place(loaded.design).legal, "reuse configured instance on persisted HHDS graph");
  std::filesystem::remove_all(path);
  auto extra = g->create_node();
  rejects([&] { apply_positions(e, [](auto, double, double) {}); });
  auto leaf = lib.create_io("unflattened");
  auto body = leaf->create_graph();
  (void)body;
  extra.set_subnode(leaf);
  rejects([&] { from_hhds(*g, layout, geom); });
  auto missing = geom;
  missing.cell = {};
  rejects([&] { from_hhds(*g, layout, missing); });
}
}  // namespace
int main() {
  try {
    test();
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
  std::cout << "HHDS adapter checks passed\n";
}
