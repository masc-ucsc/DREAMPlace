// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/dataflow.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <stdexcept>
namespace DreamPlace::native::macro {
std::vector<Flow> dataflow(const FlowGraph &g, int max_hops, int cap) {
  if (max_hops < 0 || max_hops > 16 || cap < 1)
    throw std::invalid_argument("invalid dataflow BFS controls");
  std::map<std::string, std::vector<int>> sources;
  std::map<std::string, int> names;
  auto endpoint = [&](const FlowNode &n) {
    return n.kind == FlowKind::sequential ? "@cluster:" + n.parent : n.name;
  };
  for (int i = 0; i < int(g.nodes.size()); ++i) {
    const auto &n = g.nodes[i];
    if (n.name.empty() || !names.emplace(n.name, i).second)
      throw std::invalid_argument("duplicate/empty dataflow node name");
    for (int j : n.sinks)
      if (j < 0 || j >= int(g.nodes.size()))
        throw std::invalid_argument("invalid dataflow edge");
    if (n.kind != FlowKind::combinational)
      sources[endpoint(n)].push_back(i);
  }
  std::map<std::pair<std::string, std::string>, double> weights;
  // One distance table for every source: allocating and filling a table per
  // endpoint made the pass O(sources * nodes) on large flat netlists. Each
  // search only touches the nodes within its hop budget, so reset just those.
  std::vector<int> best(g.nodes.size(), max_hops + 1);
  std::vector<int> touched;
  for (const auto &[name, starts] : sources) {
    for (int i : touched)
      best[i] = max_hops + 1;
    touched.clear();
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                        std::greater<>>
        queue;
    auto push = [&](int i, int h) {
      if (h <= max_hops && h < best[i]) {
        best[i] = h;
        touched.push_back(i);
        queue.push({h, i});
      }
    };
    for (int start : starts) {
      best[start] = 0;
      touched.push_back(start);
      const auto &out = g.nodes[start].sinks;
      if (out.size() <= size_t(cap))
        for (int sink : out)
          push(sink, 0);
    }
    while (!queue.empty()) {
      auto [h, i] = queue.top();
      queue.pop();
      if (h != best[i])
        continue;
      const auto &n = g.nodes[i];
      bool sequential = n.kind == FlowKind::sequential;
      int next = h + (sequential ? 1 : 0);
      if (next > max_hops)
        continue;
      if (n.kind != FlowKind::combinational) {
        std::string target = endpoint(n);
        if (target != name) {
          auto &weight = weights[{name, target}];
          weight = std::max(weight, std::ldexp(1., -next));
        }
        if (!sequential)
          continue;
      }
      if (n.sinks.size() <= size_t(cap))
        for (int sink : n.sinks)
          push(sink, next);
    }
  }
  std::vector<Flow> out;
  for (const auto &[ends, weight] : weights)
    out.push_back({ends.first, ends.second, weight});
  return out;
}
} // namespace DreamPlace::native::macro
