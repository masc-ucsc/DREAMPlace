// This file is distributed under the BSD 3-Clause License.
// Standalone exact-cost regression for the global-swap high-degree bounds
// cache.
#include "global_swap/src/net_bounds.h"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>

namespace DreamPlace {
template <typename T> int globalSwapCPULauncher(DetailedPlaceDB<T>, int, int);
}

namespace {
using DB = DreamPlace::DetailedPlaceDB<double>;
struct Fixture {
  std::vector<double> x, y, ox, oy, weights;
  std::vector<int> n2p, n2p_start{0}, p2n, p2net, net2p, net_start{0};
  std::vector<unsigned char> mask;
  DB db{};
  Fixture() {
    std::mt19937 random(270131);
    for (int n = 0; n < 512; ++n) {
      x.push_back(double(int(random() % 560) - 80) / 4);
      y.push_back(double(int(random() % 560) - 80) / 4);
    }
    for (int net = 0; net < 72; ++net) {
      const int count = net == 0   ? 1024
                        : net == 1 ? 128
                        : net == 2 ? 15
                        : net == 3 ? 16
                        : net == 5 ? 1
                        : net == 6 ? 0
                                   : 16 + random() % 240;
      mask.push_back(net != 4);
      weights.push_back(double(1 + random() % 1900) / 791);
      for (int i = 0; i < count; ++i) {
        // A high-degree net containing only two owners exercises the empty
        // unaffected bounds, and repeated owner pins have different offsets.
        p2n.push_back(net == 1 ? i % 2 : random() % x.size());
        p2net.push_back(net);
        ox.push_back(double(int(random() % 17) - 8) / 4);
        oy.push_back(double(int(random() % 17) - 8) / 4);
        net2p.push_back(p2n.size() - 1);
      }
      net_start.push_back(p2n.size());
    }
    for (int n = 0; n < int(x.size()); ++n) {
      for (int p = 0; p < int(p2n.size()); ++p)
        if (p2n[p] == n)
          n2p.push_back(p);
      n2p_start.push_back(n2p.size());
    }
    db.x = x.data();
    db.y = y.data();
    db.flat_net2pin_map = net2p.data();
    db.flat_net2pin_start_map = net_start.data();
    db.flat_node2pin_map = n2p.data();
    db.flat_node2pin_start_map = n2p_start.data();
    db.pin2node_map = p2n.data();
    db.pin2net_map = p2net.data();
    db.pin_offset_x = ox.data();
    db.pin_offset_y = oy.data();
    db.net_mask = mask.data();
    db.net_weights = weights.data();
    db.num_nets = mask.size();
    db.xl = db.yl = 0;
    db.xh = db.yh = 100;
  }
  double full_scan(int net, int a, double ax, double ay, int b, double bx,
                   double by) const {
    double xl = db.xh, yl = db.yh, xh = db.xl, yh = db.yl;
    for (int i = net_start[net]; i < net_start[net + 1]; ++i) {
      const int p = net2p[i], node = p2n[p];
      const double px = (node == a ? ax : node == b ? bx : x[node]) + ox[p];
      const double py = (node == a ? ay : node == b ? by : y[node]) + oy[p];
      xl = std::min(xl, px);
      xh = std::max(xh, px);
      yl = std::min(yl, py);
      yh = std::max(yh, py);
    }
    return xh - xl + yh - yl;
  }
  double original_pair_scan(int a, int b) const {
    double cost = 0;
    for (int node : {a, b}) {
      for (int i = n2p_start[node]; i < n2p_start[node + 1]; ++i) {
        const int net = p2net[n2p[i]];
        if (!mask[net])
          continue;
        if (node == b && node != a) {
          bool shared = false;
          for (int j = net_start[net]; j < net_start[net + 1]; ++j)
            shared |= p2n[net2p[j]] == a;
          if (shared)
            continue;
        }
        double hpwl = full_scan(net, a, x[a], y[a], b, x[b], y[b]);
#ifdef DREAMPLACE_TORCH_FREE
        if (db.net_weights)
          hpwl *= db.net_weights[net];
#endif
        cost += hpwl;
      }
      if (a == b)
        break;
    }
    return cost;
  }
};
void check(bool yes, const char *message) {
  if (!yes)
    throw std::runtime_error(message);
}

void check_parallel_swaps(int count) {
  constexpr int cols = 64;
  check(count >= 512 && count % cols == 0,
        "cell count must be a multiple of 64, at least 512");
  const int rows = count / cols;
  std::vector<double> x, y, width, height(count, 1), offset, weights;
  std::vector<int> p2node, p2net, net2pin, net_start{0}, node2pin,
      node_start{0};
  std::vector<int> regions(count, 0);
  std::vector<unsigned char> mask;
  for (int n = 0; n < count; ++n) {
    x.push_back(4 * (n % cols));
    y.push_back(n / cols);
    width.push_back(1 + n % 2);
  }
  auto net = [&](const std::vector<int> &nodes) {
    for (int n : nodes) {
      p2node.push_back(n);
      p2net.push_back(weights.size());
      net2pin.push_back(net2pin.size());
      offset.push_back(.25);
    }
    weights.push_back(1 + (weights.size() % 7) * .125);
    mask.push_back(1);
    net_start.push_back(p2node.size());
  };
  for (int n = 0; n < count; ++n)
    net({n, (37 * n + 137) % count});
  std::vector<int> high_degree;
  for (int n = 0; n < count; n += 4)
    high_degree.push_back(n);
  net(high_degree);
  for (int n = 0; n < count; ++n) {
    for (int p = 0; p < int(p2node.size()); ++p)
      if (p2node[p] == n)
        node2pin.push_back(p);
    node_start.push_back(node2pin.size());
  }
  auto initial_x = x, initial_y = y;
  DB db{};
  db.init_x = initial_x.data();
  db.init_y = initial_y.data();
  db.node_size_x = width.data();
  db.node_size_y = height.data();
  db.node2fence_region_map = regions.data();
  db.flat_net2pin_map = net2pin.data();
  db.flat_net2pin_start_map = net_start.data();
  db.flat_node2pin_map = node2pin.data();
  db.flat_node2pin_start_map = node_start.data();
  db.pin2node_map = p2node.data();
  db.pin2net_map = p2net.data();
  db.pin_offset_x = db.pin_offset_y = offset.data();
  db.net_weights = weights.data();
  db.net_mask = mask.data();
  db.xh = 4 * cols;
  db.yh = rows;
  db.site_width = db.row_height = 1;
  db.num_bins_x = db.num_bins_y =
      1; // enough candidates to exercise parallel scoring
  db.bin_size_x = db.xh;
  db.bin_size_y = db.yh;
  db.num_sites_x = db.xh;
  db.num_sites_y = db.yh;
  db.num_nodes = count;
  db.num_movable_nodes = count - 2;
  db.num_pins = p2node.size();
  db.num_nets = weights.size();
  std::vector<double> serial_x, serial_y;
  for (int threads : {1, 2, 4, 8}) {
    x = initial_x;
    y = initial_y;
    db.x = x.data();
    db.y = y.data();
    const double before = db.compute_total_hpwl();
    const auto start = std::chrono::steady_clock::now();
    DreamPlace::globalSwapCPULauncher(db, 2, threads);
    std::cout << "global swap " << count << " cells, " << threads
              << " threads: "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
                     .count()
              << " s\n";
    check(db.compute_total_hpwl() < before, "swap fixture did not improve");
    for (int n = db.num_movable_nodes; n < count; ++n)
      check(x[n] == initial_x[n] && y[n] == initial_y[n], "fixed cell moved");
    for (int a = 0; a < count; ++a) {
      check(x[a] >= 0 && x[a] + width[a] <= db.xh && y[a] >= 0 &&
                y[a] + 1 <= db.yh && x[a] == std::floor(x[a]),
            "swap left the legal site grid");
      for (int b = a + 1; b < count; ++b)
        check(y[a] != y[b] || x[a] + width[a] <= x[b] ||
                  x[b] + width[b] <= x[a],
              "swap created overlap");
    }
    if (threads == 1) {
      serial_x = x;
      serial_y = y;
    } else {
      check(x == serial_x && y == serial_y,
            "parallel scoring changed serial move selection");
    }
  }
}
} // namespace

int main(int argc, char **argv) {
  try {
    check_parallel_swaps(argc > 1 ? std::stoi(argv[1]) : 512);
    Fixture f;
    DreamPlace::global_swap_detail::NetBounds<double> cache(f.db);
    check(cache.cached(0) && cache.cached(1) && cache.cached(3),
          "high-degree nets not cached");
    check(!cache.cached(2) && !cache.cached(4),
          "low-degree or masked net cached");
    std::mt19937 random(271732);
    size_t tested = 0;
    size_t originals = 0;
    for (int trial = 0; trial < 3000; ++trial) {
      const int a = trial % 7 == 0 ? 0 : random() % f.x.size();
      const int b = trial % 7 == 0 ? 1 : random() % f.x.size();
      check(cache.compute_original_pair_hpwl(a, b) ==
                f.original_pair_scan(a, b),
            "cached original pair differs from weighted full scan");
      check(cache.compute_original_pair_hpwl(a, a) ==
                f.original_pair_scan(a, a),
            "identical-node original pair differs from full scan");
      originals += 2;
      if (a == b)
        continue;
      const double ax = double(int(random() % 640) - 120) / 4;
      const double ay = double(int(random() % 640) - 120) / 4;
      const double bx = double(int(random() % 640) - 120) / 4;
      const double by = double(int(random() % 640) - 120) / 4;
      for (int net = 0; net < f.db.num_nets; ++net) {
        if (f.mask[net]) {
          check(cache.original_hpwl(net) ==
                    f.full_scan(net, a, f.x[a], f.y[a], b, f.x[b], f.y[b]),
                "cached original net differs from full pin scan");
          ++originals;
        }
        if (!cache.cached(net))
          continue;
        check(cache.hpwl(net, a, ax, ay, b, bx, by) ==
                  f.full_scan(net, a, ax, ay, b, bx, by),
              "cached proposal differs from full pin scan");
        bool found = false;
        for (int i = f.net_start[net]; i < f.net_start[net + 1]; ++i)
          found |= f.p2n[f.net2p[i]] == a;
        check(cache.contains(a, net) == found, "shared-net membership differs");
        ++tested;
      }
      if (trial % 11 == 0) {
        f.x[a] = ax;
        f.y[a] = ay;
        f.x[b] = bx;
        f.y[b] = by;
        cache.refresh(a, b);
      }
    }
    // Include ties among more than three distinct owners and proposals that
    // remove all extrema owners at the core boundaries.
    std::fill(f.x.begin(), f.x.end(), 0);
    std::fill(f.y.begin(), f.y.end(), 100);
    for (int n = 0; n < int(f.x.size()); ++n)
      cache.refresh(n, n);
    for (int net = 0; net < f.db.num_nets; ++net)
      if (f.mask[net]) {
        check(cache.original_hpwl(net) ==
                  f.full_scan(net, 0, f.x[0], f.y[0], 1, f.x[1], f.y[1]),
              "tied original extrema differ");
        if (cache.cached(net))
          check(cache.hpwl(net, 0, -50, 150, 1, 150, -50) ==
                    f.full_scan(net, 0, -50, 150, 1, 150, -50),
                "tied extrema differ");
      }
    f.db.net_weights = nullptr;
    check(cache.compute_original_pair_hpwl(0, 1) == f.original_pair_scan(0, 1),
          "unweighted original pair differs from full scan");
    std::cout << "global swap bounds: " << tested
              << " exact proposal comparisons, " << originals
              << " exact original-cost comparisons passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
