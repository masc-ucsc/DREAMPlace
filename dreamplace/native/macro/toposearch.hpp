// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/compact.hpp"
namespace DreamPlace::native::macro {
struct SearchResult {
  Realization best;
  SequencePair sequence;
  double cost = 0, seed_cost = 0, rudy_overflow = 0;
  int evaluations = 0, accepted = 0, lp_solves = 0, screened = 0;
};
SearchResult toposearch(const Design &, std::span<const NodeInfo>,
                        const Realization &, const Parameters &,
                        const SequencePair &previous = {});
} // namespace DreamPlace::native::macro
