// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"
#include <span>
#include <vector>

namespace DreamPlace::native::macro {
enum class ViolationKind { bounds, pitch, overlap };
struct Violation {
  ViolationKind kind;
  int node, other = -1;
};
struct CheckResult {
  std::vector<Violation> violations;
  bool legal() const { return violations.empty(); }
};
// Explicit macro subset. Fixed nodes are checked identically to movable ones.
// Physical envelopes include per-side halos rounded up to site/row pitches.
// The caller maps transient indices to stable names for persisted diagnostics.
CheckResult check(const Design &design, std::span<const int> nodes,
                  double tolerance = 1e-6);
} // namespace DreamPlace::native::macro
