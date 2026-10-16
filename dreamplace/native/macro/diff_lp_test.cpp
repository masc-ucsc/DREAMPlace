// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/diff_lp.hpp"
#include "macro_legalize/src/native_lp.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace DreamPlace::native::macro;
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("difference LP oracle mismatch");
}
int main() {
  try {
    for (int a = -2; a <= 4; ++a)
      for (int b = -2; b <= 4; ++b)
        for (int weight = 1; weight <= 3; ++weight) {
          DifferenceLP p;
          p.variable(0, 4, 0, "x");
          p.variable(0, 4, 0, "y");
          p.anchor(0, a, weight);
          p.anchor(1, b, 2);
          p.constrain(1, 0, -2, 0, "separate");
          auto s = p.solve();
          double best = 1e100;
          for (int x = 0; x <= 4; ++x)
            for (int y = 0; y <= 4; ++y)
              if (x + 2 <= y)
                best = std::min(best, double(weight * std::abs(x - a) +
                                             2 * std::abs(y - b)));
          check(std::abs(s.objective - best) < 1e-9);
        }
    // Linear HPWL endpoint costs, soft separation and L1 displacement.
    // Exhaustive enumeration independently checks dual supply signs and
    // soft-arc capacities.
    for (int weight = 1; weight <= 3; ++weight) {
      DifferenceLP p;
      p.variable(0, 4, -weight, "l");
      p.variable(0, 4, weight, "u");
      p.anchor(0, 3, 2);
      p.anchor(1, 1, 1);
      p.constrain(1, 0, -2, 3, "soft");
      auto s = p.solve();
      double best = 1e100;
      for (int x = 0; x <= 4; ++x)
        for (int y = 0; y <= 4; ++y)
          best = std::min(best,
                          double(weight * (y - x) + 2 * std::abs(x - 3) +
                                 std::abs(y - 1) + 3 * std::max(0, x - y + 2)));
      check(std::abs(s.objective - best) < 1e-9);
    }
    for (int a = -2; a <= 5; ++a) {
      DifferenceLP p;
      p.variable(0, 6);
      p.variable(0, 6);
      p.anchor(0, a, 1);
      p.anchor(1, 4, 1);
      p.constrain(1, 0, -2);
      auto s = p.solve();
      DreamPlace::MacroDifferenceLP old;
      old.lower = {0, 0};
      old.upper = {6, 6};
      old.anchor = {a, 4};
      old.constraints = {{1, 0, -2}};
      auto v = old.solve();
      check(s.objective == std::abs(v[0] - a) + std::abs(v[1] - 4));
    }
    DifferenceLP bad;
    bad.variable(0, 5, 0, "a");
    bad.variable(0, 5, 0, "b");
    bad.constrain(0, 1, -1, 0, "ab");
    bad.constrain(1, 0, -1, 0, "ba");
    bool caught = false;
    try {
      bad.solve();
    } catch (const Infeasible &e) {
      caught = true;
      check(e.cycle.size() == 2);
    }
    check(caught);
    bad.arcs[0].penalty = 10;
    check(!bad.solve().slack.empty());
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
