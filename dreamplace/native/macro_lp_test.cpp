// BSD-3-Clause. See LICENSE.
#include <cstdlib>
#include <iostream>
#include <limits>

#include "macro_legalize/src/native_lp.h"

int main() {
  using LP = DreamPlace::MacroDifferenceLP;
  // Independent exhaustive oracle for small bounded integer LPs. Covers
  // positive/negative anchors, active bounds, separation and origin offsets.
  for (int a = -3; a <= 5; ++a)
    for (int b = -3; b <= 5; ++b) {
      LP lp;
      lp.lower = {-2, -2, -2};
      lp.upper = {4, 4, 4};
      lp.anchor = {a, b, 1};
      lp.constraints = {{1, 0, -1}, {2, 1, -2}, {0, 2, 5}};
      long long best = std::numeric_limits<long long>::max();
      for (int x = -2; x <= 4; ++x)
        for (int y = -2; y <= 4; ++y)
          for (int z = -2; z <= 4; ++z)
            if (x + 1 <= y && y + 2 <= z && z - x <= 5)
              best = std::min(best, static_cast<long long>(std::abs(x - a) + std::abs(y - b) + std::abs(z - 1)));
      auto v = lp.solve();
      if (std::abs(v[0] - a) + std::abs(v[1] - b) + std::abs(v[2] - 1) != best) {
        std::cerr << "macro LP disagrees with exhaustive optimum\n";
        return 1;
      }
    }
  LP bad;
  bad.lower = {0, 0};
  bad.upper = {1, 1};
  bad.anchor = {0, 0};
  bad.constraints = {{0, 1, -2}};
  try {
    bad.solve();
  } catch (const std::runtime_error&) {
    return 0;
  }
  std::cerr << "infeasible macro LP accepted\n";
  return 1;
}
