// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/check.hpp"
#include "dreamplace/native/macro/dataflow.hpp"
#include "dreamplace/native/macro/h2d.hpp"
#include "dreamplace/native/macro/report.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void check(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
struct Fixture {
  Design d;
  std::vector<m::NodeInfo> info;
  m::FlowGraph flow;
};
Fixture pipeline() {
  Fixture f;
  auto &d = f.d;
  d.layout = {0, 0, 120, 40, 1, 1};
  d.num_movable = 6;
  d.node_x = {10, 50, 90, 24, 64, 104, 0, 120};
  d.node_y = {16, 16, 16, 20, 20, 20, 20, 20};
  d.node_size_x = {10, 10, 10, 2, 2, 2, 0, 0};
  d.node_size_y = {8, 8, 8, 1, 1, 1, 0, 0};
  d.node_halo_x = {1, 1, 1, 0, 0, 0, 0, 0};
  d.node_halo_y = d.node_halo_x;
  for (int i = 0; i < 3; ++i)
    f.info.push_back({"s" + std::to_string(i) + ".ram", "s" + std::to_string(i),
                      "ram", true});
  for (int i = 0; i < 3; ++i)
    f.info.push_back(
        {"s" + std::to_string(i) + ".logic", "s" + std::to_string(i), "buf"});
  f.info.push_back({"@io:in", {}, {}, false, false, true});
  f.info.push_back({"@io:out", {}, {}, false, false, true});
  const int chain[] = {6, 0, 3, 1, 4, 2, 5, 7};
  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k < 2; ++k) {
      int node = chain[n + k];
      d.flat_net2pin.push_back(d.pin2node.size());
      d.pin2node.push_back(node);
      d.pin_offset_x.push_back(k ? 0 : d.node_size_x[node]);
      d.pin_offset_y.push_back(d.node_size_y[node] / 2);
    }
    d.flat_net2pin_start.push_back(d.flat_net2pin.size());
  }
  d.build_maps();
  for (const auto &i : f.info)
    f.flow.nodes.push_back({i.name,
                            i.parent,
                            i.macro ? m::FlowKind::macro
                            : i.io  ? m::FlowKind::io
                                    : m::FlowKind::combinational,
                            {}});
  for (int i = 0; i < 7; ++i)
    f.flow.nodes[chain[i]].sinks.push_back(chain[i + 1]);
  return f;
}
int main() {
  try {
    auto f = pipeline();
    m::Parameters p;
    p.K = 1;
    p.iters_per_outer = 5;
    p.seed_sa_evals = 8;
    p.evals = 8;
    p.bins = 8;
    m::Topology input;
    input.design = "pipeline3";
    input.dataflow = m::dataflow(f.flow);
    Options cells;
    cells.max_iterations = 10;
    cells.bins_x = cells.bins_y = 8;
    auto single = f.d;
    auto r1 = m::place_h2d(single, f.info, p, input, cells);
    check(r1.legal && is_legal(single), "pipeline illegal");
    check(r1.topology.groups.size() == 3, "hierarchy grouping lost");
    auto again = f.d;
    auto repeat = m::place_h2d(again, f.info, p, input, cells);
    check(single.node_x == again.node_x && single.node_y == again.node_y &&
              single.node_orientation == again.node_orientation,
          "H2D nondeterministic");
    p.K = 2;
    auto multiple = f.d;
    auto r2 = m::place_h2d(multiple, f.info, p, input, cells);
    check(r2.merit <= r1.merit + 1e-8, "best iterate worse than K=1");
    check(r2.iterations.size() >= 1, "no H2D iterations");
    const auto report = m::report_json(r2, multiple, f.info);
    auto cached = m::cached_topology(report, f.d, f.info, r2.resolved);
    check(cached.groups.size() == 3, "topology cache lost groups");
    auto restored = m::abstract_design(f.d, f.info, cached, p, true);
    auto expanded = f.d;
    m::expand(restored, cached, expanded);
    for (int i = 0; i < 3; ++i)
      check(expanded.node_x[i] == multiple.node_x[i] &&
                expanded.node_y[i] == multiple.node_y[i],
            "macro cache moved coordinates");
    // The coarse-pitch recompaction runs its own cells pass. Like the main
    // rounds it must keep solver frames private (one snapshot per round) and
    // apply the same resolved controls.
    {
      auto coarse = m::Parameters(p);
      coarse.K = 1;
      coarse.pitch_x = coarse.pitch_y = 2;
      Options observed = cells;
      std::vector<PlacementStage> stages;
      observed.observer = [&](const Design &, PlacementStage stage, int) {
        stages.push_back(stage);
      };
      auto pitched = f.d;
      auto pr = m::place_h2d(pitched, f.info, coarse, input, observed);
      check(pr.legal && is_legal(pitched), "coarse-pitch pipeline illegal");
      check(stages.size() == 1 && stages[0] == PlacementStage::global_placement,
            "coarse-pitch pass leaked internal solver frames");
    }
    // Whole hierarchical bus ports resolve to actual fixed physical bit IOs.
    auto bus_info = f.info;
    bus_info[6].name = "@io:din[0]";
    bus_info[7].name = "@io:din[1]";
    m::Topology bus;
    bus.dataflow = {{bus_info[0].name, "@io:din", 2}};
    auto bt = m::seed_topology(f.d, bus_info, {}, p, bus);
    check(bt.dataflow.size() == 2 && bt.dataflow[0].to != "@io:din" &&
              bt.dataflow[1].to != "@io:din",
          "bus IO bundle lost");
    auto ba = m::abstract_design(f.d, bus_info, bt, p);
    auto before_nets = ba.design.net_weights.size();
    m::add_pseudonets(ba, bt, p);
    check(ba.design.net_weights.size() == before_nets + 2,
          "bus IO pseudo-nets lost");
    p.K = 1;
    p.search = "sa";
    p.chains = 2;
    p.evals = 6;
    auto annealed = f.d;
    auto rs = m::place_h2d(annealed, f.info, p, input, cells);
    check(rs.legal, "SA fixture illegal");
    p.toposearch = false;
    p.pseudonets = false;
    p.anchors = false;
    p.order = false;
    auto ablated = f.d;
    check(m::place_h2d(ablated, f.info, p, input, cells).legal,
          "ablations illegal");
    // Topology search may permute free groups but must keep a frozen group's
    // relation to every other group in both sequences.
    {
      Fixture row;
      row.d.layout = {0, 0, 160, 50, 1, 1};
      row.d.num_movable = 5;
      for (int i = 0; i < 5; ++i) {
        row.d.node_x.push_back(4 + 30 * i);
        row.d.node_y.push_back(20);
        row.d.node_size_x.push_back(10);
        row.d.node_size_y.push_back(8);
        row.d.node_halo_x.push_back(1);
        row.d.node_halo_y.push_back(1);
        row.info.push_back({"r" + std::to_string(i) + ".m",
                            "r" + std::to_string(i), "ram", true});
      }
      for (int i = 0; i < 4; ++i)
        for (int pin : {i, i + 1}) {
          row.d.flat_net2pin.push_back(row.d.pin2node.size());
          row.d.pin2node.push_back(pin);
          row.d.pin_offset_x.push_back(5);
          row.d.pin_offset_y.push_back(4);
          if (pin == i + 1)
            row.d.flat_net2pin_start.push_back(row.d.flat_net2pin.size());
        }
      row.d.build_maps();
      m::Parameters sp;
      sp.seed_sa_evals = 0;
      sp.anchors = sp.order = false;
      sp.search = "sa";
      sp.evals = 60;
      sp.bins = 8;
      auto topology = m::seed_topology(row.d, row.info, {}, sp);
      topology.groups[2].freeze_topology = true;
      const auto abstracted = m::abstract_design(row.d, row.info, topology, sp);
      auto incumbent = m::group_legalize(abstracted, topology, sp);
      check(incumbent.legal, "frozen search incumbent illegal");
      const auto initial = m::sequence_pair(incumbent.topology);
      const auto search = m::toposearch(row.d, row.info, incumbent, sp);
      check(search.evaluations > 0, "topology search did not run");
      const std::string frozen = topology.groups[2].id;
      auto side = [&](const std::vector<std::string> &order,
                      const std::string &id) {
        std::vector<bool> before;
        bool seen = false;
        for (const auto &other : order) {
          if (other == id)
            seen = true;
          else
            before.push_back(!seen);
        }
        return before;
      };
      auto others = [&](const std::vector<std::string> &order) {
        std::vector<std::string> out;
        for (const auto &id : order)
          if (id != frozen)
            out.push_back(id);
        return out;
      };
      for (auto [was, now] :
           {std::pair{&initial.plus, &search.sequence.plus},
            std::pair{&initial.minus, &search.sequence.minus}}) {
        const auto key = others(*was);
        std::map<std::string, bool> was_before, now_before;
        auto a = side(*was, frozen), b = side(*now, frozen);
        const auto now_others = others(*now);
        for (size_t k = 0; k < key.size(); ++k)
          was_before[key[k]] = a[k];
        for (size_t k = 0; k < now_others.size(); ++k)
          now_before[now_others[k]] = b[k];
        check(was_before == now_before, "search crossed a frozen group");
      }
    }
    // Factor-grid grouping is based on hierarchy/master, not flat node ids.
    Fixture array;
    array.d.layout = {0, 0, 100, 100, 1, 1};
    array.d.num_movable = 16;
    for (int i = 0; i < 16; ++i) {
      array.d.node_x.push_back(10 + (i % 4) * 12);
      array.d.node_y.push_back(10 + (i / 4) * 12);
      array.d.node_size_x.push_back(8);
      array.d.node_size_y.push_back(8);
      array.d.node_halo_x.push_back(1);
      array.d.node_halo_y.push_back(1);
      array.info.push_back(
          {"banks.m" + std::to_string(i), "banks", "ram", true});
    }
    array.d.build_maps();
    auto topology = m::seed_topology(array.d, array.info, {}, p);
    check(topology.groups.size() == 1 && topology.groups[0].cols == 4 &&
              topology.groups[0].rows == 4,
          "array shape curve failed");
    p.seed_sa_evals = 0;
    p.iters_per_outer = 0;
    cells.max_iterations = 0;
    auto ar = m::place_h2d(array.d, array.info, p, {}, cells);
    check(ar.legal, "rigid array illegal");
    // A user-fixed macro keeps its exact position and pins across all stages.
    auto fixed = f.d;
    fixed.num_movable = 2; // nodes 2..end are fixed, including std cells
    auto before = fixed;
    p.K = 1;
    p.iters_per_outer = 0;
    p.evals = 0;
    p.search = "none";
    auto fr = m::place_h2d(fixed, f.info, p, input, cells);
    check(fr.legal, "fixed fixture illegal");
    for (int i = 2; i < 8; ++i)
      check(fixed.node_x[i] == before.node_x[i] &&
                fixed.node_y[i] == before.node_y[i],
            "fixed input moved");
    // Fixing one leaf in an inferred array creates a fixed singleton instead
    // of accidentally fixing or rejecting the entire hierarchy/master bucket.
    auto guided = f.d;
    auto guided_info = f.info;
    guided_info[1].parent = guided_info[0].parent;
    m::Topology guidance;
    guidance.fixed = {{guided_info[0].name, 4, 4, Orientation::FN, {true, 1}}};
    auto gr = m::place_h2d(guided, guided_info, p, guidance, cells);
    check(gr.legal && guided.node_x[0] == 4 && guided.node_y[0] == 4 &&
              guided.node_orientation[0] == Orientation::FN,
          "typed fixed leaf ignored");
    // Impossible packing is an error and never publishes a partially placed
    // netlist.
    Fixture bad;
    bad.d.layout = {0, 0, 10, 10, 1, 1};
    bad.d.num_movable = 4;
    for (int i = 0; i < 4; ++i) {
      bad.d.node_x.push_back(1);
      bad.d.node_y.push_back(1);
      bad.d.node_size_x.push_back(8);
      bad.d.node_size_y.push_back(8);
      bad.info.push_back({"g" + std::to_string(i) + ".m",
                          "g" + std::to_string(i), "ram", true});
    }
    bad.d.build_maps();
    auto original = bad.d;
    bool failed = false;
    try {
      m::place_h2d(bad.d, bad.info, p, {}, cells);
    } catch (const std::exception &) {
      failed = true;
    }
    check(failed && bad.d.node_x == original.node_x &&
              bad.d.node_y == original.node_y,
          "infeasible placement published partial result");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
