// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/fillers.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
void require(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
int main() {
  try {
    Design d;
    d.layout = {0, 0, 1500, 1500, .05, .25};
    d.num_movable = 1;
    d.node_x = {0, 1000};
    d.node_y = {0, 1000};
    d.node_size_x = {.5, 600};
    d.node_size_y = {.25, 600};
    d.node_halo_x = {0, .11};
    d.node_halo_y = {0, .26};
    d.build_maps();
    Options o;
    o.target_density = .5;
    const double w = m::bounded_filler_width(d, o, .5);
    // Fixed box extends beyond the core: its clipped padded area is
    // 500.15 * 500.5, not its full 600.3 * 601 footprint.
    const double target = .5 * (1500. * 1500. - 500.15 * 500.5) - .125;
    const double count = std::floor(target / (w * .25));
    require(w > .5 && count <= 1000000 && count >= 999999,
            "sparse design exceeds filler budget or loses target density");
    require(target - count * w * .25 < w * .25 + 1e-6,
            "filler area does not represent target density");
    require(d.node_size_x[0] == .5 && d.node_size_x[1] == 600,
            "physical sizes changed");
    d.layout.xh = d.layout.yh = 100;
    require(m::bounded_filler_width(d, o, .5) == .5,
            "ordinary preferred width changed");
    o.add_fillers = false;
    require(m::bounded_filler_width(d, o, .5) == .5,
            "disabled fillers changed");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
