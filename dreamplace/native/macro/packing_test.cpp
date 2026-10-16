// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/check.hpp"
#include "dreamplace/native/macro/packing.hpp"
#include <iostream>
#include <numeric>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void require(bool b, const char *s) {
  if (!b)
    throw std::runtime_error(s);
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 1000, 800, 1, 1};
    d.num_movable = 96;
    for (int i = 0; i < 96; ++i) {
      const double widths[] = {50, 72, 40, 220, 12, 57},
                   heights[] = {120, 80, 35, 24, 22, 42};
      d.node_x.push_back(2);
      d.node_y.push_back(2);
      d.node_size_x.push_back(widths[i % 6] + .3);
      d.node_size_y.push_back(heights[i % 6] + .7);
      d.node_halo_x.push_back(2);
      d.node_halo_y.push_back(2);
    }
    d.node_x.push_back(320);
    d.node_y.push_back(260);
    d.node_size_x.push_back(160);
    d.node_size_y.push_back(160);
    d.node_halo_x.push_back(2);
    d.node_halo_y.push_back(2);
    d.build_maps();
    auto original = d;
    require(m::pack_rectangles(d), "mixed groups did not pack");
    std::vector<int> nodes(d.node_x.size());
    std::iota(nodes.begin(), nodes.end(), 0);
    require(m::check(d, nodes).legal(), "packing not independently legal");
    require(d.node_x.back() == 320 && d.node_y.back() == 260 &&
                d.node_size_x == original.node_size_x &&
                d.node_size_y == original.node_size_y,
            "fixed input or physical sizes changed");
    auto repeat = original;
    require(m::pack_rectangles(repeat) && repeat.node_x == d.node_x &&
                repeat.node_y == d.node_y,
            "packing nondeterministic");
    Design bad;
    bad.layout = {0, 0, 10, 10, 1, 1};
    bad.num_movable = 2;
    bad.node_x = {0, 0};
    bad.node_y = {0, 0};
    bad.node_size_x = {8, 8};
    bad.node_size_y = {8, 8};
    bad.build_maps();
    auto before = bad;
    require(!m::pack_rectangles(bad) && bad.node_x == before.node_x &&
                bad.node_y == before.node_y,
            "failed packing modified input");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
