// BSD-3-Clause. See LICENSE.
#pragma once
#include "utility/src/detailed_place_db.h"
#include <algorithm>
#include <array>
#include <vector>

DREAMPLACE_BEGIN_NAMESPACE
// Exact proposal bounds for high-degree nets. Excluding two moved nodes needs
// only the first three extrema with distinct owners; repeated pins from an
// owner cannot evict the third unaffected owner. Low-degree nets keep the
// direct scan.
namespace global_swap_detail {
template <typename T> class NetBounds {
  struct Entry {
    T value{};
    int node = -1;
  };
  struct Extreme {
    std::array<Entry, 3> entries{};
    int count = 0;
    void add(T value, int node, bool minimum) {
      auto better = [&](T a, T b) { return minimum ? a < b : a > b; };
      int at = 0;
      while (at < count && entries[at].node != node)
        ++at;
      if (at < count) {
        if (!better(value, entries[at].value))
          return;
        entries[at].value = value;
      } else {
        if (count == 3 && !better(value, entries[2].value))
          return;
        at = std::min(count, 2);
        entries[at] = {value, node};
        count = std::min(count + 1, 3);
      }
      while (at > 0 && better(entries[at].value, entries[at - 1].value)) {
        std::swap(entries[at], entries[at - 1]);
        --at;
      }
    }
    T excluding(int a, int b, T sentinel) const {
      for (int i = 0; i < count; ++i)
        if (entries[i].node != a && entries[i].node != b)
          return entries[i].value;
      return sentinel;
    }
  };
  struct Bounds {
    Extreme xmin, xmax, ymin, ymax;
  };
  const DetailedPlaceDB<T> &db;
  std::vector<int> index;
  std::vector<Bounds> bounds;
  std::vector<T> current_hpwl;
  std::vector<size_t> refreshed;
  size_t generation = 0;

  void rebuild_current(int net) {
    DreamPlace::Utility::Box<T> box(db.xh, db.yh, db.xl, db.yl);
    for (int i = db.flat_net2pin_start_map[net];
         i < db.flat_net2pin_start_map[net + 1]; ++i) {
      const int pin = db.flat_net2pin_map[i];
      const int node = db.pin2node_map[pin];
      const T x = db.x[node] + db.pin_offset_x[pin];
      const T y = db.y[node] + db.pin_offset_y[pin];
      box.xl = std::min(box.xl, x);
      box.xh = std::max(box.xh, x);
      box.yl = std::min(box.yl, y);
      box.yh = std::max(box.yh, y);
    }
    current_hpwl[net] = box.xh - box.xl + box.yh - box.yl;
  }

  void rebuild(int net) {
    auto &b = bounds[index[net]];
    b = {};
    for (int i = db.flat_net2pin_start_map[net];
         i < db.flat_net2pin_start_map[net + 1]; ++i) {
      const int pin = db.flat_net2pin_map[i];
      const int node = db.pin2node_map[pin];
      const T x = db.x[node] + db.pin_offset_x[pin];
      const T y = db.y[node] + db.pin_offset_y[pin];
      b.xmin.add(x, node, true);
      b.xmax.add(x, node, false);
      b.ymin.add(y, node, true);
      b.ymax.add(y, node, false);
    }
  }

public:
  explicit NetBounds(const DetailedPlaceDB<T> &database)
      : db(database), index(db.num_nets, -1), current_hpwl(db.num_nets),
        refreshed(db.num_nets) {
    for (int net = 0; net < db.num_nets; ++net) {
      if (!db.net_mask[net])
        continue;
      rebuild_current(net);
      if (db.flat_net2pin_start_map[net + 1] - db.flat_net2pin_start_map[net] >=
          16) {
        index[net] = bounds.size();
        bounds.emplace_back();
        rebuild(net);
      }
    }
  }
  bool cached(int net) const { return index[net] >= 0; }
  T original_hpwl(int net) const { return current_hpwl[net]; }
  bool contains(int node, int net) const {
    for (int i = db.flat_node2pin_start_map[node];
         i < db.flat_node2pin_start_map[node + 1]; ++i)
      if (db.pin2net_map[db.flat_node2pin_map[i]] == net)
        return true;
    return false;
  }
  T compute_original_pair_hpwl(int a, int b) const {
    T cost = 0;
    auto add = [&](int net) {
      T hpwl = current_hpwl[net];
#ifdef DREAMPLACE_TORCH_FREE
      if (db.net_weights)
        hpwl *= db.net_weights[net];
#endif
      cost += hpwl;
    };
    // Preserve the full-scan cost's incident-pin order and repeated-net
    // contributions, including the target's exclusion of shared nets.
    for (int i = db.flat_node2pin_start_map[a];
         i < db.flat_node2pin_start_map[a + 1]; ++i) {
      const int net = db.pin2net_map[db.flat_node2pin_map[i]];
      if (db.net_mask[net])
        add(net);
    }
    for (int i = db.flat_node2pin_start_map[b];
         i < db.flat_node2pin_start_map[b + 1]; ++i) {
      const int net = db.pin2net_map[db.flat_node2pin_map[i]];
      if (db.net_mask[net] && !contains(a, net))
        add(net);
    }
    return cost;
  }
  T hpwl(int net, int a, T ax, T ay, int b, T bx, T by) const {
    const auto &cached = bounds[index[net]];
    DreamPlace::Utility::Box<T> box(db.xh, db.yh, db.xl, db.yl);
    // Keep the original core-bound sentinels even if every remaining pin lies
    // outside the core or there are no unaffected owners.
    box.xl = std::min(box.xl, cached.xmin.excluding(a, b, db.xh));
    box.xh = std::max(box.xh, cached.xmax.excluding(a, b, db.xl));
    box.yl = std::min(box.yl, cached.ymin.excluding(a, b, db.yh));
    box.yh = std::max(box.yh, cached.ymax.excluding(a, b, db.yl));
    auto add = [&](int node, T x, T y) {
      for (int i = db.flat_node2pin_start_map[node];
           i < db.flat_node2pin_start_map[node + 1]; ++i) {
        const int pin = db.flat_node2pin_map[i];
        if (db.pin2net_map[pin] != net)
          continue;
        box.xl = std::min(box.xl, x + db.pin_offset_x[pin]);
        box.xh = std::max(box.xh, x + db.pin_offset_x[pin]);
        box.yl = std::min(box.yl, y + db.pin_offset_y[pin]);
        box.yh = std::max(box.yh, y + db.pin_offset_y[pin]);
      }
    };
    add(a, ax, ay);
    if (b != a)
      add(b, bx, by);
    return box.xh - box.xl + box.yh - box.yl;
  }
  void refresh(int a, int b) {
    ++generation;
    for (int node : {a, b})
      for (int i = db.flat_node2pin_start_map[node];
           i < db.flat_node2pin_start_map[node + 1]; ++i) {
        const int net = db.pin2net_map[db.flat_node2pin_map[i]];
        if (db.net_mask[net] && refreshed[net] != generation) {
          rebuild_current(net);
          if (cached(net))
            rebuild(net);
          refreshed[net] = generation;
        }
      }
  }
};
} // namespace global_swap_detail

DREAMPLACE_END_NAMESPACE
