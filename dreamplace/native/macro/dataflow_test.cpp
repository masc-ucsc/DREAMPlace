// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/dataflow.hpp"
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native::macro;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("dataflow oracle mismatch");
}
int main() {
  try {
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
