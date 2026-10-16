// BSD-3-Clause. See LICENSE.
#include <cstdlib>
#include <iostream>
#include <limits>

#include "utility/src/legalization_db.h"
#include "macro_legalize/src/lp_legalize.h"

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
  // X has a valid displacement solution, but Y has infeasible bounds.
  // A failed heuristic LP must leave BOTH coordinates intact so the macro
  // launcher's grid fallback starts from the same placement.
  double x = -2, y = 1, anchor_x = 3, anchor_y = 1, width = 4, height = 12;
  DreamPlace::LegalizationDB<double> db{};
  db.x = &x;
  db.y = &y;
  db.init_x = &anchor_x;
  db.init_y = &anchor_y;
  db.node_size_x = &width;
  db.node_size_y = &height;
  db.xl = db.yl = 0;
  db.xh = db.yh = 10;
  DreamPlace::lpLegalizeGraphLauncher(db, {0}, {});
  if (x != -2 || y != 1) {
    std::cerr << "infeasible macro LP partially changed coordinates\n";
    return 1;
  }
  LP bad;
  bad.lower = {0, 0};
  bad.upper = {1, 1};
  bad.anchor = {0, 0};
  bad.constraints = {{0, 1, -2}};
  try {
    bad.solve();
  } catch (const DreamPlace::MacroLPInfeasible&) {
    return 0;
  }
  std::cerr << "infeasible macro LP accepted\n";
  return 1;
}
