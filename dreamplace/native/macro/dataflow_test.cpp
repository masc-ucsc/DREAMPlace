// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/dataflow.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <stdexcept>
using namespace DreamPlace::native::macro;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("dataflow oracle mismatch");
}
// Straightforward per-source reference of the same hop-distance semantics.
std::vector<Flow> reference(const FlowGraph &g, int max_hops, int cap) {
  std::map<std::string, std::vector<int>> sources;
  auto endpoint = [&](const FlowNode &n) {
    return n.kind == FlowKind::sequential ? "@cluster:" + n.parent : n.name;
  };
  for (int i = 0; i < int(g.nodes.size()); ++i)
    if (g.nodes[i].kind != FlowKind::combinational)
      sources[endpoint(g.nodes[i])].push_back(i);
  std::map<std::pair<std::string, std::string>, double> weights;
  for (const auto &[name, starts] : sources) {
    std::vector<int> best(g.nodes.size(), max_hops + 1);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                        std::greater<>>
        queue;
    auto push = [&](int i, int h) {
      if (h <= max_hops && h < best[i]) {
        best[i] = h;
        queue.push({h, i});
      }
    };
    for (int start : starts) {
      best[start] = 0;
      if (g.nodes[start].sinks.size() <= size_t(cap))
        for (int sink : g.nodes[start].sinks)
          push(sink, 0);
    }
    while (!queue.empty()) {
      auto [h, i] = queue.top();
      queue.pop();
      if (h != best[i])
        continue;
      const auto &n = g.nodes[i];
      const bool sequential = n.kind == FlowKind::sequential;
      const int next = h + (sequential ? 1 : 0);
      if (next > max_hops)
        continue;
      if (n.kind != FlowKind::combinational) {
        if (endpoint(n) != name) {
          auto &weight = weights[{name, endpoint(n)}];
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
void randomized() {
  std::mt19937 rng(7);
  for (int round = 0; round < 40; ++round) {
    FlowGraph g;
    const int n = 20 + rng() % 200;
    for (int i = 0; i < n; ++i) {
      const auto kind = static_cast<FlowKind>(rng() % 4);
      g.nodes.push_back({"n" + std::to_string(i), "p" + std::to_string(rng() % 6),
                         kind, {}});
    }
    for (int e = 0; e < 3 * n; ++e)
      g.nodes[rng() % n].sinks.push_back(rng() % n);
    for (int hops : {0, 1, 3, 4})
      for (int cap : {1, 3, 4096}) {
        const auto got = dataflow(g, hops, cap), want = reference(g, hops, cap);
        check(got.size() == want.size());
        for (size_t i = 0; i < got.size(); ++i)
          check(got[i].from == want[i].from && got[i].to == want[i].to &&
                got[i].weight == want[i].weight);
      }
  }
}
int main() {
  try {
    randomized();
    FlowGraph g;
    g.nodes = {{"m0", "a", FlowKind::macro, {1}},
               {"comb", "a", FlowKind::combinational, {2, 2}},
               {"reg0", "regs0", FlowKind::sequential, {3}},
               {"reg1", "regs1", FlowKind::sequential, {4}},
               {"m1", "b", FlowKind::macro, {}}};
    auto flows = dataflow(g);
    bool found = false;
    for (const auto &f : flows)
      if (f.from == "m0" && f.to == "m1") {
        check(f.weight == .25);
        found = true;
      }
    check(found);
    for (const auto &f : dataflow(g, 1))
      check(f.from != "m0" || f.to != "m1");
    g.nodes[1].sinks.push_back(1);
    flows = dataflow(g);
    for (const auto &f : flows)
      if (f.from == "m0" && f.to == "m1")
        check(f.weight == .25);
    g.nodes[0].sinks = {1, 2, 3, 4};
    for (const auto &f : dataflow(g, 4, 2))
      check(f.from != "m0");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
