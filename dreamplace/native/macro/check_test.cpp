// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/check.hpp"
#include <array>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
namespace {
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("macro checker assertion failed");
}
template <typename F> void rejects(F f) {
  bool bad = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    bad = true;
  }
  check(bad);
}
bool has(const m::CheckResult &r, m::ViolationKind kind) {
  for (const auto &v : r.violations)
    if (v.kind == kind)
      return true;
  return false;
}
Design design() {
  Design d;
  d.layout = {10, 20, 110, 80, 2, 3};
  d.num_movable = 2;
  d.node_x = {12, 36, 70, 90, 10};
  d.node_y = {23, 47, 23, 50, 20};
  d.node_size_x = {12, 12, 8, 8, 0};
  d.node_size_y = {12, 12, 8, 8, 0};
  d.node_halo_x = {1.1, 1.1, 1.1, 1.1, 0};
  d.node_halo_y = {2.2, 2.2, 2.2, 2.2, 0};
  d.build_maps();
  return d;
}
} // namespace
int main() {
  try {
    const std::array ids{0, 1, 2, 3};
    const std::array fixed{2, 3};
    auto d = design();
    check(m::check(d, ids).legal());
    d.node_x[3] = 74;
    d.node_y[3] = 23;
    check(has(m::check(d, fixed), m::ViolationKind::overlap));
    d = design();
    d.node_x[1] = 24;
    d.node_y[1] = 23;
    check(has(m::check(d, ids),
              m::ViolationKind::overlap)); // physical edges only touch
    d = design();
    d.node_x[3] = 102;
    check(has(m::check(d, ids),
              m::ViolationKind::bounds)); // only halo outside core
    d = design();
    d.node_x[3] += .5;
    check(has(m::check(d, ids), m::ViolationKind::pitch));
    d = design();
    const auto before = d;
    rejects([&] { m::check(d, std::array{0, 0}); });
    rejects([&] { m::check(d, std::array{4}); });
    rejects([&] { m::check(d, std::array{-1}); });
    rejects([&] { m::check(d, ids, -1); });
    check(d.node_x == before.node_x && d.node_y == before.node_y);
    // Independent pair oracle over a moving fixed macro on the technology grid.
    for (int x = 12; x <= 98; x += 2)
      for (int y = 23; y <= 68; y += 3) {
        d = design();
        d.node_x[3] = x;
        d.node_y[3] = y;
        bool overlap = false;
        for (int i = 0; i < 4; ++i)
          for (int j = i + 1; j < 4; ++j)
            overlap |= d.node_x[i] - 2 < d.node_x[j] + d.node_size_x[j] + 2 &&
                       d.node_x[j] - 2 < d.node_x[i] + d.node_size_x[i] + 2 &&
                       d.node_y[i] - 3 < d.node_y[j] + d.node_size_y[j] + 3 &&
                       d.node_y[j] - 3 < d.node_y[i] + d.node_size_y[i] + 3;
        check(has(m::check(d, ids), m::ViolationKind::overlap) == overlap);
      }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
