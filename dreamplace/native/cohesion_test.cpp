// BSD-3-Clause. See LICENSE. Checks remain active under -c opt.
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool value, const char *message) {
  if (!value)
    throw std::runtime_error(message);
}
void near(double a, double b, double tolerance, const char *message) {
  check(std::abs(a - b) <= tolerance, message);
}
template <class F> void rejects(F f) {
  bool caught = false;
  try {
    f();
  } catch (const std::invalid_argument &) {
    caught = true;
  }
  check(caught, "invalid cohesion accepted");
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 1000, 1000, 1, 1};
  d.num_movable = 2;
  d.node_x = {100, 800, 700};
  d.node_y = {200, 700, 800};
  d.node_size_x = {4, 8, 5};
  d.node_size_y = {2, 3, 4};
  d.cohesion_groups = {{{0, 1}, 2, {}}, {{2}, 1, {}}};
  d.build_maps();
  d.validate();
  return d;
}
int main() {
  try {
    auto d = fixture();
    const double lambda = .3;
    auto grad = detail::cohesion(d, lambda);
    near(grad.x[0] + grad.x[1], 0, 1e-10, "cohesion translates group");
    near(grad.y[0] + grad.y[1], 0, 1e-10, "cohesion translates group");
    near(grad.x[2], 0, 0, "singleton force");
    for (int axis = 0; axis < 2; ++axis)
      for (int i = 0; i < 2; ++i) {
        auto hi = d, lo = d;
        auto &hp = axis ? hi.node_y : hi.node_x;
        auto &lp = axis ? lo.node_y : lo.node_x;
        const double eps = 1e-3;
        hp[i] += eps;
        lp[i] -= eps;
        double fd = (detail::cohesion(hi, lambda).value -
                     detail::cohesion(lo, lambda).value) /
                    (2 * eps);
        near(fd, axis ? grad.y[i] : grad.x[i], 1e-5,
             "centroid derivative disagrees with finite difference");
      }
    auto translated = d;
    for (auto &x : translated.node_x)
      x += 120;
    for (auto &y : translated.node_y)
      y -= 50;
    near(detail::cohesion(translated, lambda).value, grad.value, 1e-9,
         "translation changes cohesion energy");
    auto padded = d;
    for (auto &g : padded.cohesion_groups)
      for (int i : g.nodes)
        g.areas.push_back(d.node_size_x[i] * d.node_size_y[i]);
    for (int i = 0; i < 3; ++i) {
      padded.node_x[i] -= 10;
      padded.node_y[i] -= 10;
      padded.node_size_x[i] += 20;
      padded.node_size_y[i] += 20;
    }
    near(detail::cohesion(padded, lambda).value, grad.value, 1e-9,
         "halo changes physical-area cohesion");
    auto bad = d;
    bad.cohesion_groups[0].nodes.push_back(0);
    rejects([&] { bad.validate(); });
    bad = d;
    bad.cohesion_groups[0].nodes[0] = 100;
    rejects([&] { bad.validate(); });
    bad = d;
    bad.cohesion_groups[0].weight = -1;
    rejects([&] { bad.validate(); });
    Options o;
    o.max_iterations = 20;
    o.bins_x = o.bins_y = 8;
    o.add_fillers = false;
    o.legalize = o.detailed_place = false;
    o.initial_density_weight = 1e-8;
    auto plain = d;
    plain.cohesion_groups.clear();
    place(plain, o);
    auto disabled = d;
    place(disabled, o);
    check(plain.node_x == disabled.node_x && plain.node_y == disabled.node_y,
          "disabled cohesion changes GP");
    auto enabled = d;
    o.cohesion_weight = 1;
    place(enabled, o);
    plain.cohesion_groups = d.cohesion_groups;
    check(detail::cohesion(enabled, 1).value < detail::cohesion(plain, 1).value,
          "cohesion fails to reduce scatter");
    check(enabled.node_x[2] == d.node_x[2] && enabled.node_y[2] == d.node_y[2],
          "cohesion moves fixed node");
    o.cohesion_weight = -1;
    rejects([&] { place(enabled, o); });
    std::cout << "cohesion gradient, translation, halo weights, opt-in GP and "
                 "fixed-node checks passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
