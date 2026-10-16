// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/abstract.hpp"
#include "dreamplace/native/macro/diff_lp.hpp"
namespace DreamPlace::native::macro {

struct Realization {
  Topology topology;
  std::vector<Separation> graph;
  std::vector<std::string> relaxed;
  std::vector<DifferenceLP::Slack> slack;
  double hpwl = 0, slack_amount = 0, anchor_deviation = 0, order_violation = 0;
  int lp_solves = 0;
  bool legal = false;
};
SequencePair sequence_pair(const Topology &);
std::vector<Separation> sequence_graph(const Topology &, const SequencePair &);
std::vector<Separation> position_graph(const Abstracted &, const Topology &);
// Fixed obstacles are added from the current placement, never permuted/moved.
Realization compact(const Abstracted &, const Topology &, const Parameters &,
                    std::vector<Separation>, bool displacement_only = false);
Realization group_legalize(const Abstracted &, const Topology &,
                           const Parameters &);
double anchor_deviation(const Topology &, const Layout &);
double order_violation(const Topology &);
} // namespace DreamPlace::native::macro
