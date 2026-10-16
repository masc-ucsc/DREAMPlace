// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/abstract.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("abstract assertion failed");
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 80, 40, 1, 1};
    d.num_movable = 3;
    d.node_x = {10, 24, 50, 0};
    d.node_y = {10, 10, 10, 20};
    d.node_size_x = {10, 10, 2, 0};
    d.node_size_y = {8, 8, 1, 0};
    d.node_halo_x = {2, 2, 0, 0};
    d.node_halo_y = {1, 1, 0, 0};
    d.pin2node = {0, 1, 1, 2, 0, 3};
    d.pin_offset_x = {1, 2, 3, 1, 9, 0};
    d.pin_offset_y = {1, 2, 3, .5, 7, 0};
    d.flat_net2pin = {0, 1, 2, 3, 4, 5};
    d.flat_net2pin_start = {0, 2, 4, 6};
    d.build_maps();
    std::vector<m::NodeInfo> info = {{"a.m0", "a", "ram", true},
                                     {"a.m1", "a", "ram", true},
                                     {"b.cell", "b", "inv"},
                                     {"@io:in", {}, {}, false, false, true}};
    m::Group g;
    g.id = "a";
    g.x = 10;
    g.y = 10;
    g.width = 24;
    g.height = 8;
    g.halo_x = 2;
    g.halo_y = 1;
    g.cols = 2;
    g.members = {{"a.m0", 0, 0}, {"a.m1", 14, 0}};
    m::Topology t;
    t.groups = {g};
    m::Parameters p;
    auto a = m::abstract_design(d, info, t, p);
    check(a.design.node_x.size() == 3);
    check(a.design.net_mask[0] == 0);
    Design back = d;
    m::expand(a, t, back);
    check(back.node_x == d.node_x && back.node_y == d.node_y &&
          back.pin_offset_x == d.pin_offset_x);
    auto original_mask = d.net_mask;
    d.net_mask[0] = 0;
    check(std::abs(hpwl(a.design) - hpwl(d)) < 1e-9);
    d.net_mask = original_mask;
    for (auto orientation :
         {Orientation::N, Orientation::FN, Orientation::FS, Orientation::S}) {
      t.groups[0].orientation = orientation;
      auto oriented = m::abstract_design(d, info, t, p);
      auto expanded = d;
      m::expand(oriented, t, expanded);
      expanded.net_mask[0] = 0;
      check(std::abs(hpwl(oriented.design) - hpwl(expanded)) < 1e-9);
      check(expanded.node_orientation[0] == orientation);
    }
    t.groups[0].orientation = Orientation::N;
    t.dataflow = {{"a.m0", "@io:in", 4}, {"a.m1", "b.cell", 2}};
    p.theta_q = 0;
    m::add_pseudonets(a, t, p);
    check(a.design.net_weights.size() == 5);
    check(a.design.node_x.size() == 4 && a.design.num_movable == 2);
    auto merged = m::merge(t, {}, d, info);
    check(!merged.macro_fingerprint.empty());
    for (const auto &s : m::parameter_specs()) {
      m::Parameters parsed;
      m::set_parameter(parsed, s.name, s.default_value);
      check(m::parameter_value(parsed, s.name) ==
                m::parameter_value(p, s.name) ||
            s.name == "theta_q");
    }
    p.validate();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
