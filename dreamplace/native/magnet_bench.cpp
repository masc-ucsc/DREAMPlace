// BSD-3-Clause. See LICENSE.
// Bounded geometry-only microbenchmark for candidate/grid/rigid-node costs.
#include "dreamplace/native/magnet.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
int main(int argc, char **argv) {
  try {
    const int cells = argc > 1 ? std::stoi(argv[1]) : 200000;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 20;
    if (argc > 3 || cells < 0 || cells > 2000000 || repeats < 1 ||
        repeats > 1000)
      throw std::invalid_argument(
          "magnet_bench [cells 0..2000000] [repeats 1..1000]");
    Design source;
    source.layout = {0, 0, 1000, 3000, 1, 1};
    source.num_movable = 1040 + cells;
    for (int i = 0; i < source.num_movable; ++i) {
      const bool macro = i < 1040;
      const int k = i - 1040;
      source.node_x.push_back(macro ? 20 + (i % 40) * 20
                                    : 10 + (k * 37LL % 980));
      source.node_y.push_back(macro ? 20 + (i / 40) * 60
                                    : 1800 + (k * 71LL % 1100));
      source.node_size_x.push_back(macro ? 10 : .5);
      source.node_size_y.push_back(macro ? 10 : 1);
      source.node_halo_x.push_back(macro ? 2 : 0);
      source.node_halo_y.push_back(macro ? 2 : 0);
      source.node_macro.push_back(macro);
    }
    source.build_maps();
    double grid_s = 0, search_s = 0, merge_s = 0, expand_s = 0;
    for (int repeat = 0; repeat < repeats; ++repeat) {
      const auto start = std::chrono::steady_clock::now();
      MagnetDensity density(source, 128, 128);
      const auto gridded = std::chrono::steady_clock::now();
      const auto found = magnet_candidates(source, {}, density);
      const auto searched = std::chrono::steady_clock::now();
      auto d = source;
      RigidMacroMap map;
      const auto result = merge_macros(d, map, found.merges);
      if (!result.applied || found.merges.size() != 520)
        throw std::runtime_error("microbenchmark failed rigid batch: " +
                                 result.reason);
      const auto merged = std::chrono::steady_clock::now();
      const auto expanded = expand_macros(d, source, map);
      if (expanded.node_x.back() != source.node_x.back() ||
          expanded.node_y.back() != source.node_y.back())
        throw std::runtime_error("microbenchmark moved an untouched cell");
      const auto end = std::chrono::steady_clock::now();
      grid_s += std::chrono::duration<double>(gridded - start).count();
      search_s += std::chrono::duration<double>(searched - gridded).count();
      merge_s += std::chrono::duration<double>(merged - searched).count();
      expand_s += std::chrono::duration<double>(end - merged).count();
    }
    std::cout << "{\"macros\":1040,\"cells\":" << cells
              << ",\"repeats\":" << repeats
              << ",\"pairs_per_round\":520,\"grid_s\":" << grid_s / repeats
              << ",\"search_s\":" << search_s / repeats
              << ",\"merge_s\":" << merge_s / repeats
              << ",\"expand_s\":" << expand_s / repeats << "}\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
