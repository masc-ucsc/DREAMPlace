// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/check.hpp"
#include "dreamplace/native/macro/h2d.hpp"
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void require(bool b, const char *s) {
  if (!b)
    throw std::runtime_error(s);
}
struct Fixture {
  Design d;
  std::vector<m::NodeInfo> info;
};
Fixture macros(int n, double w, double h, double cw, double ch,
               bool same_parent = false) {
  Fixture f;
  f.d.layout = {0, 0, cw, ch, 1, 1};
  f.d.num_movable = n;
  for (int i = 0; i < n; ++i) {
    auto parent = same_parent ? std::string("array") : "g" + std::to_string(i);
    f.info.push_back({parent + ".m" + std::to_string(i), parent, "ram", true});
    f.d.node_x.push_back(1 + (i % 4) * (w + 2));
    f.d.node_y.push_back(1 + (i / 4) * (h + 2));
    f.d.node_size_x.push_back(w);
    f.d.node_size_y.push_back(h);
    f.d.node_halo_x.push_back(1);
    f.d.node_halo_y.push_back(1);
  }
  f.d.build_maps();
  return f;
}
void net(Design &d, int a, int b) {
  for (int i : {a, b}) {
    d.flat_net2pin.push_back(d.pin2node.size());
    d.pin2node.push_back(i);
    d.pin_offset_x.push_back(d.node_size_x[i]);
    d.pin_offset_y.push_back(d.node_size_y[i] / 2);
  }
  d.flat_net2pin_start.push_back(d.flat_net2pin.size());
}
int main() {
  try {
    m::Parameters p;
    p.K = 2;
    p.bins = 8;
    p.iters_per_outer = 5;
    p.seed_sa_evals = 8;
    p.evals = 24;
    p.anchors = p.order = false;
    p.min_channel = 1;
    Options cells;
    cells.max_iterations = 5;
    cells.bins_x = cells.bins_y = 8;
    auto sides = macros(4, 10, 8, 80, 80);
    for (auto [x, y] :
         {std::pair{0., 40.}, {80., 40.}, {40., 0.}, {40., 80.}}) {
      int id = sides.info.size();
      sides.info.push_back(
          {"@io:" + std::to_string(id), {}, {}, false, false, true});
      sides.d.node_x.push_back(x);
      sides.d.node_y.push_back(y);
      sides.d.node_size_x.push_back(0);
      sides.d.node_size_y.push_back(0);
      sides.d.node_halo_x.push_back(0);
      sides.d.node_halo_y.push_back(0);
    }
    for (int i = 0; i < 4; ++i)
      net(sides.d, i, i + 4);
    sides.d.build_maps();
    auto st = m::seed_topology(sides.d, sides.info, {}, p);
    for (int i = 0; i < 4; ++i)
      st.anchors.push_back({"edge" + std::to_string(i),
                            st.groups[i].id,
                            static_cast<Side>(i),
                            0,
                            {true, 1}});
    auto sr = m::place_h2d(sides.d, sides.info, p, st, cells);
    require(sr.legal && m::anchor_deviation(sr.topology, sides.d.layout) == 0,
            "edge IO hard anchor lost");
    // Ring and scrambled ring are robustness fixtures, not an unmeasured QoR
    // discrimination claim. Repeated seeds must be exact in the serial build.
    for (bool scrambled : {false, true}) {
      auto ring = macros(8, 8, 8, 80, 80);
      for (int i = 0; i < 8; ++i)
        net(ring.d, i, (i + 1) % 8);
      ring.d.build_maps();
      if (scrambled)
        for (int i = 0; i < 4; ++i)
          std::swap(ring.d.node_x[i], ring.d.node_x[7 - i]);
      auto original = ring.d;
      auto r = m::place_h2d(ring.d, ring.info, p, {}, cells);
      auto repeat = original;
      m::place_h2d(repeat, ring.info, p, {}, cells);
      require(r.legal && repeat.node_x == ring.d.node_x &&
                  repeat.node_y == ring.d.node_y &&
                  repeat.node_orientation == ring.d.node_orientation,
              "ring legality/determinism");
    }
    auto wall = macros(16, 16, 16, 80, 80, true);
    auto wr = m::place_h2d(wall.d, wall.info, p, {}, cells);
    require(wr.legal && wr.topology.groups.size() == 1 &&
                wr.topology.groups[0].cols == 4 &&
                wr.topology.groups[0].rows == 4,
            "64 percent SRAM array failed");
    // Three edge anchors cannot share a side this short. Soft mode reports the
    // compromise; hard mode must reject the unchanged input with named causes.
    auto tight = macros(3, 8, 8, 40, 20);
    auto original = tight.d;
    p.iters_per_outer = 0;
    p.seed_sa_evals = 0;
    p.evals = 0;
    p.search = "none";
    cells.max_iterations = 0;
    auto tt = m::seed_topology(tight.d, tight.info, {}, p);
    for (int i = 0; i < 3; ++i)
      tt.anchors.push_back({"tight" + std::to_string(i),
                            tt.groups[i].id,
                            Side::left,
                            0,
                            {false, 1}});
    auto tr = m::place_h2d(tight.d, tight.info, p, tt, cells);
    require(tr.legal && !tr.violations.empty(),
            "tight soft anchors not reported");
    for (auto &a : tt.anchors)
      a.strength.hard = true;
    bool rejected = false;
    try {
      m::place_h2d(original, tight.info, p, tt, cells);
    } catch (const std::runtime_error &e) {
      rejected = std::string(e.what()).find("tight") != std::string::npos;
    }
    require(rejected && original.node_x == std::vector<double>({1, 11, 21}),
            "tight hard anchors not diagnosed/transactional");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
