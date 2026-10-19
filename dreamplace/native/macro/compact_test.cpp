// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/compact.hpp"
#include "dreamplace/native/macro/seed.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void require(bool v, const char *message) {
  if (!v)
    throw std::runtime_error(message);
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 40, 20, 1, 1};
    d.num_movable = 2;
    d.node_x = {24, 4};
    d.node_y = {4, 4};
    d.node_size_x = {8, 8};
    d.node_size_y = {4, 4};
    d.node_halo_x = {1, 1};
    d.node_halo_y = {1, 1};
    d.build_maps();
    std::vector<m::NodeInfo> info = {{"a.m", "a", "ram", true},
                                     {"b.m", "b", "ram", true}};
    m::Parameters p;
    p.seed_sa_evals = 0;
    p.anchors = p.order = false;
    auto t = m::seed_topology(d, info, {}, p);
    auto aid = t.groups[0].id, bid = t.groups[1].id;
    t.anchors = {{"left", aid, Side::left, 0, {true, 1}},
                 {"right", bid, Side::right, 0, {true, 1}}};
    t.order = {{"a-before-b", aid, bid, Axis::x, {true, 1}}};
    auto a = m::abstract_design(d, info, t, p);
    auto legal = m::group_legalize(a, t, p);
    require(legal.legal && legal.topology.groups[0].x == 1 &&
                legal.topology.groups[1].x == 31,
            "hard anchors/order lost");
    t.order.push_back({"b-before-a", bid, aid, Axis::x, {true, 1}});
    bool cycle = false;
    try {
      m::group_legalize(m::abstract_design(d, info, t, p), t, p);
    } catch (const m::Infeasible &e) {
      cycle = !e.cycle.empty();
    }
    require(cycle, "hard cycle not diagnosed");
    t.order[0].strength = {false, 1};
    t.order[1].strength = {false, 2};
    t.anchors.clear();
    p.order_in_lp = true;
    auto soft = m::compact(m::abstract_design(d, info, t, p), t, p, {}, true);
    require(soft.legal &&
                soft.relaxed == std::vector<std::string>{"a-before-b"},
            "wrong soft edge relaxed");
    // Exact HPWL oracle on a one-group/two-net instance. No l/u auxiliaries
    // are needed; enumerate all integer legal x positions independently.
    Design one;
    one.layout = {0, 0, 40, 10, 1, 1};
    one.num_movable = 1;
    one.node_x = {2, 20};
    one.node_y = {2, 2};
    one.node_size_x = {4, 0};
    one.node_size_y = {4, 0};
    one.pin2node = {0, 1, 0, 1};
    one.pin_offset_x = {1, 0, 3, 0};
    one.pin_offset_y = {2, 0, 2, 0};
    one.flat_net2pin = {0, 1, 2, 3};
    one.flat_net2pin_start = {0, 2, 4};
    one.build_maps();
    one.net_weights = {1, 3};
    std::vector<m::NodeInfo> oi = {{"g.m", "g", "ram", true},
                                   {"@io:out", {}, {}, false, false, true}};
    p.order_in_lp = false;
    p.min_channel = 1;
    auto ot = m::seed_topology(one, oi, {}, p);
    auto oa = m::abstract_design(one, oi, ot, p);
    auto opt = m::compact(oa, ot, p, {});
    auto expanded = one;
    m::expand(m::abstract_design(one, oi, opt.topology, p), opt.topology,
              expanded);
    double oracle = 1e100;
    for (int x = 1; x <= 35; ++x) {
      auto candidate = one;
      candidate.node_x[0] = x;
      candidate.node_y[0] = 1;
      oracle = std::min(oracle, hpwl(candidate));
    }
    require(opt.legal && std::abs(hpwl(expanded) - oracle) < 1e-8,
            "single-group LP differs from exhaustive HPWL");
    // Fixed objects retain residuals when the requested snap pitch is coarser.
    d.num_movable = 1;
    d.node_x = {3, 25};
    d.node_y = {3, 5};
    t = m::seed_topology(d, info, {}, p);
    p.pitch_x = 4;
    p.pitch_y = 2;
    auto snapped =
        m::compact(m::abstract_design(d, info, t, p), t, p,
                   m::position_graph(m::abstract_design(d, info, t, p), t));
    require(snapped.legal && snapped.topology.groups[1].x == 25 &&
                snapped.topology.groups[1].y == 5,
            "coarse snap moved a fixed object");
    // Cache guidance keys include strengths and constraints, and unnamed leaves
    // can be rematched by a size/pin signature and ordinal under the same
    // parent.
    auto key = m::topology_key(t);
    auto changed = t;
    changed.groups[0].freeze_topology = true;
    require(key != m::topology_key(changed), "guidance hash collision");
    info[0].stable = false;
    auto old = m::remap_unstable(t, d, info);
    info[0].name = "a.generated_99";
    auto remapped = m::merge({}, old, d, info);
    require(remapped.groups[0].members[0].name == info[0].name &&
                !remapped.diagnostics.empty(),
            "unstable fallback failed");
    // Many edge-crowded rigid groups exercise the recovery path used when
    // a prototype has collapsed its soft anchors onto one side.
    Design crowded;
    crowded.layout = {0, 0, 1000, 800, 1, 1};
    crowded.num_movable = 48;
    std::vector<m::NodeInfo> ci;
    for (int i = 0; i < 48; ++i) {
      const double widths[] = {50, 72, 40, 220, 12, 57},
                   heights[] = {120, 80, 35, 24, 22, 42};
      crowded.node_x.push_back(2);
      crowded.node_y.push_back(2);
      crowded.node_size_x.push_back(widths[i % 6] + .3);
      crowded.node_size_y.push_back(heights[i % 6] + .7);
      crowded.node_halo_x.push_back(2);
      crowded.node_halo_y.push_back(2);
      ci.push_back({"cluster" + std::to_string(i) + ".m",
                    "cluster" + std::to_string(i), "ram", true});
    }
    crowded.build_maps();
    p.pitch_x = p.pitch_y = 0;
    p.seed_sa_evals = 0;
    auto ct = m::seed_topology(crowded, ci, {}, p);
    for (const auto &g : ct.groups)
      ct.anchors.push_back({"edge:" + g.id, g.id, Side::left, 0, {false, 1}});
    auto cr = m::group_legalize(m::abstract_design(crowded, ci, ct, p), ct, p);
    require(cr.legal, "edge-crowded group recovery failed");
    // A frozen group keeps the realized relation even when the geometric
    // prototype would pick the opposite side; obstacle relations to fixed
    // objects are always retained, and a hard order owns its pair.
    {
      Design row;
      row.layout = {0, 0, 120, 30, 1, 1};
      row.num_movable = 3;
      row.node_x = {10, 40, 70, 100};
      row.node_y = {10, 10, 10, 10};
      row.node_size_x = {10, 10, 10, 10};
      row.node_size_y = {8, 8, 8, 8};
      row.build_maps();
      std::vector<m::NodeInfo> ri = {{"a.m", "a", "ram", true},
                                     {"b.m", "b", "ram", true},
                                     {"c.m", "c", "ram", true},
                                     {"f.m", "f", "ram", true}};
      m::Parameters rp;
      rp.seed_sa_evals = 0;
      rp.anchors = rp.order = false;
      rp.min_channel = 1;
      auto rt = m::seed_topology(row, ri, {}, rp);
      const auto ida = rt.groups[0].id, idb = rt.groups[1].id,
                 idc = rt.groups[2].id;
      auto position = [&](const m::Topology &topology, const std::string &id) {
        for (const auto &g : topology.groups)
          if (g.id == id)
            return g.x;
        return -1.;
      };
      // Geometry says a < b < c; freeze b with the opposite relation to a.
      rt.groups[1].freeze_topology = true;
      rt.realized_graph = {{idb, ida, "sp:" + idb + ":" + ida, Axis::x}};
      auto ra = m::abstract_design(row, ri, rt, rp);
      auto realized = m::compact(ra, rt, rp, m::position_graph(ra, rt));
      require(realized.legal && position(realized.topology, idb) <
                                    position(realized.topology, ida),
              "frozen realized relation ignored");
      require(position(realized.topology, idb) + 10 <=
                  position(realized.topology, ida),
              "frozen relation left an overlap");
      // The geometric prototype still separates the unfrozen pairs, and the
      // fixed group (index 3) remains an obstacle for every movable group.
      for (const auto &id : {ida, idb, idc}) {
        const double x = position(realized.topology, id);
        require(x + 10 <= 100 || x >= 110, "fixed obstacle not retained");
      }
      bool pinned = false;
      for (const auto &e : realized.graph)
        pinned |= e.from == idb && e.to == ida;
      require(pinned, "realized edge missing from the returned graph");
      for (const auto &e : realized.graph)
        require(!(e.from == ida && e.to == idb),
                "opposite prototype edge survived for a frozen pair");
    }
    // LPbench: 50 groups and 5,000 physical nets exercise multigroup l/u rows.
    Design large;
    large.layout = {0, 0, 1000, 100, 1, 1};
    large.num_movable = 50;
    std::vector<m::NodeInfo> li;
    for (int i = 0; i < 50; ++i) {
      large.node_x.push_back(i * 15 + 2);
      large.node_y.push_back(10);
      large.node_size_x.push_back(8);
      large.node_size_y.push_back(8);
      li.push_back({"g" + std::to_string(i) + ".m", "g" + std::to_string(i),
                    "ram", true});
    }
    for (int n = 0; n < 5000; ++n) {
      for (int v : {n % 50, (n * 7 + 3) % 50}) {
        large.flat_net2pin.push_back(large.pin2node.size());
        large.pin2node.push_back(v);
        large.pin_offset_x.push_back(4);
        large.pin_offset_y.push_back(4);
      }
      large.flat_net2pin_start.push_back(large.flat_net2pin.size());
    }
    large.build_maps();
    p.pitch_x = p.pitch_y = 0;
    auto lt = m::seed_topology(large, li, {}, p);
    auto la = m::abstract_design(large, li, lt, p);
    auto lr = m::compact(la, lt, p, m::position_graph(la, lt));
    require(lr.legal && std::isfinite(lr.hpwl), "LPbench illegal");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
