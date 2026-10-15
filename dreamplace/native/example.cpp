// BSD-3-Clause. See LICENSE.
#include <iostream>

#include "dreamplace/native/placement.hpp"

int main() {
  using namespace DreamPlace::native;
  Design d;
  d.layout = {0, 0, 40, 10, 1, 1};
  d.num_movable = 4;
  d.node_x = {20, 20, 20, 20, 0, 40};
  d.node_y = {5, 5, 5, 5, 5, 5};
  d.node_size_x = {2, 2, 2, 2, 0, 0};
  d.node_size_y = {1, 1, 1, 1, 0, 0};
  d.pin2node = {4, 0, 0, 1, 1, 2, 2, 3, 3, 5};
  d.pin_offset_x.assign(10, 0);
  d.pin_offset_y.assign(10, 0);
  d.flat_net2pin = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  d.flat_net2pin_start = {0, 2, 4, 6, 8, 10};
  d.build_maps();
  Placer placer;
  placer.options().max_iterations = 500;
  placer.options().target_density = 0.8;
  placer.options().initial_density_weight = 0.01;
  auto result = placer.place(d);
  std::cout << "{\"hpwl\":" << result.hpwl << ",\"overflow\":" << result.overflow
            << ",\"legal\":" << (result.legal ? "true" : "false") << ",\"iterations\":" << result.iterations << "}\n";
  return result.legal ? 0 : 1;
}
