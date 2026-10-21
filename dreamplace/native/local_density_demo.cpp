// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
using namespace DreamPlace::native;
int main(int argc, char **) {
  Design source;
  source.layout = {0, 0, 100, 100, .1, .25};
  source.num_movable = 1600;
  std::mt19937 rng(1000);
  std::uniform_real_distribution<double> u(-2, 2);
  for (int i = 0; i < 1600; ++i) {
    source.node_x.push_back(25 + 50 * (i / 400 % 2) + u(rng));
    source.node_y.push_back(25 + 50 * (i / 800) + u(rng));
    source.node_size_x.push_back(.5);
    source.node_size_y.push_back(.25);
  }
  for (int i = 0; i < 4; ++i) {
    source.node_x.push_back(25 + 50 * (i % 2));
    source.node_y.push_back(25 + 50 * (i / 2));
    source.node_size_x.push_back(0);
    source.node_size_y.push_back(0);
  }
  for (int i = 0; i < 1600; ++i) {
    source.pin2node.push_back(i);
    source.pin2node.push_back(1600 + i / 400);
    source.pin_offset_x.insert(source.pin_offset_x.end(), {.25, 0});
    source.pin_offset_y.insert(source.pin_offset_y.end(), {.125, 0});
    source.flat_net2pin.push_back(2 * i);
    source.flat_net2pin.push_back(2 * i + 1);
    source.flat_net2pin_start.push_back(2 * i + 2);
  }
  source.node_macro.assign(source.node_x.size(), 0);
  source.build_maps();
  std::cout << "weight,bins,seconds,hpwl,local_seconds,median_local_"
               "utilization,max_local_utilization,force,virtual_fraction,"
               "virtual_count\n";
  for (bool dynamic :
       (argc > 1 ? std::vector<bool>{true} : std::vector<bool>{false, true}))
    for (double weight : (argc > 1 ? std::vector<double>{4.}
                                   : std::vector<double>{0., 2., 4., 8., 12.}))
      for (double fraction : (argc > 1 ? std::vector<double>{0., .5, 1.}
                                       : std::vector<double>{0.})) {
        auto d = source;
        Options o;
        o.bins_x = o.bins_y = 16;
        o.dynamic_bins = dynamic;
        o.bin_max = 128;
        o.auxiliary_bins = 16;
        o.legalize = o.detailed_place = false;
        o.add_fillers = fraction > 0;
        o.filler_area_fraction = fraction > 0 ? fraction : 1;
        o.filler_max_count = 4096;
        o.min_iterations = o.max_iterations = 200;
        o.local_density_weight = weight;
        o.local_density_bins = 128;
        o.target_density = .5;
        o.adaptive_density = true;
        const auto start = std::chrono::steady_clock::now();
        const auto r = place(d, o);
        std::vector<double> density(50 * 50);
        for (int i = 0; i < 1600; ++i) {
          const double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
                       h = d.node_size_y[i];
          for (int bx = std::max(0, int(x / 2));
               bx < std::min(50, int(std::ceil((x + w) / 2))); ++bx)
            for (int by = std::max(0, int(y / 2));
                 by < std::min(50, int(std::ceil((y + h) / 2))); ++by)
              density[bx * 50 + by] +=
                  std::max(0., std::min(x + w, 2. * bx + 2) -
                                   std::max(x, 2. * bx)) *
                  std::max(0., std::min(y + h, 2. * by + 2) -
                                   std::max(y, 2. * by)) /
                  4;
        }
        std::vector<double> local;
        for (int i = 0; i < 1600; ++i)
          local.push_back(
              density[std::clamp(int((d.node_x[i] + .25) / 2), 0, 49) * 50 +
                      std::clamp(int((d.node_y[i] + .125) / 2), 0, 49)]);
        std::sort(local.begin(), local.end());
        std::cout << weight << ',' << r.bins_x << ','
                  << std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - start)
                         .count()
                  << ',' << r.hpwl << ',' << r.local_density_seconds << ','
                  << local[local.size() / 2] << ','
                  << *std::max_element(density.begin(), density.end()) << ','
                  << r.local_density_force << ',' << fraction << ','
                  << r.filler_count << std::endl;
      }
}
