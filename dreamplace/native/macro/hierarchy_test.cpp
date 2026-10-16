// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/hierarchy.hpp"
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native::macro;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("hierarchical connectivity mismatch");
}
int main() {
  try {
    hhds::GraphLibrary lib;
    auto leaf = lib.create_io("leaf");
    leaf->add_input("a", 1);
    leaf->add_output("z", 2);
    auto child = leaf->create_graph();
    auto memory = child->create_node();
    memory.set_name("mem");
    child->get_input_pin("a").connect_sink(memory.create_sink_pin(1));
    memory.create_driver_pin(2).connect_sink(child->get_output_pin("z"));
    auto module = lib.create_io("module");
    module->add_input("a", 1);
    module->add_output("z", 2);
    auto body = module->create_graph();
    auto nested = body->create_node();
    nested.set_name("inner");
    nested.set_subnode(leaf);
    body->get_input_pin("a").connect_sink(nested.create_sink_pin(1));
    nested.create_driver_pin(2).connect_sink(body->get_output_pin("z"));
    auto topio = lib.create_io("top");
    topio->add_input("a", 1);
    topio->add_output("z", 2);
    auto top = topio->create_graph();
    auto first = top->create_node(), second = top->create_node();
    first.set_name("u0");
    first.set_subnode(module);
    second.set_name("u1");
    second.set_subnode(module);
    top->get_input_pin("a").connect_sink(first.create_sink_pin(1));
    first.create_driver_pin(2).connect_sink(second.create_sink_pin(1));
    second.create_driver_pin(2).connect_sink(top->get_output_pin("z"));
    HierarchyCallbacks cb;
    cb.instance_name = [](auto n) { return std::string(n.get_name()); };
    cb.kind = [](auto n) {
      return n.get_name() == "mem" ? FlowKind::macro : FlowKind::combinational;
    };
    cb.opaque = [](auto n) { return n.get_name() == "mem"; };
    cb.terminal_name = [](auto p) {
      return p.get_master_node().is_input_node() ? "a" : "z";
    };
    auto h = read_hierarchy(*top, cb);
    check(h.flow.nodes.size() == 4);
    auto flow = dataflow(h.flow);
    bool input = false, between = false, output = false;
    for (const auto &f : flow) {
      input |= f.from == "@io:a" && f.to == "u0.inner.mem";
      between |= f.from == "u0.inner.mem" && f.to == "u1.inner.mem";
      output |= f.from == "u1.inner.mem" && f.to == "@io:z";
    }
    check(input && between && output); // No spurious input->output shortcut
                                       // across the two module interfaces.
    for (const auto &f : flow)
      check(f.from != "@io:a" || f.to != "@io:z");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
