// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <source_location>
#include <stdexcept>
using namespace DreamPlace::native;
namespace {
void check(bool ok,
           std::source_location loc = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("macro routing check at line " +
                             std::to_string(loc.line()));
}
Design fixture() {
  Design d;
  d.layout = {0, 0, 1000, 1000, 1, 1};
  d.num_movable = 2;
  d.node_x = {400, 490, 485, 485};
  d.node_y = {200, 200, 50, 950};
  d.node_size_x = {80, 80, 0, 0};
  d.node_size_y = {600, 600, 0, 0};
  d.node_macro = {1, 1, 0, 0};
  d.pin2node = {2, 3};
  d.pin_offset_x = d.pin_offset_y = {0, 0};
  d.flat_net2pin = {0, 1};
  d.flat_net2pin_start = {0, 2};
  d.net_weights = {1};
  d.build_maps();
  return d;
}
} // namespace
int main() {
  try {
    auto d = fixture();
    std::vector<double> demand(64 * 64, .8);
    auto f = detail::macro_routing_pressure(d, 64, 64, demand, 2);
    check(f.channels == 1 && f.x[0] > 0 && f.x[1] == -f.x[0]);
    check(f.y[0] == 0 && f.y[1] == 0 && f.x[2] == 0);
    const double expected = .8 * .8 * std::pow(1 - 10. / 31.25, 2);
    check(std::abs(f.x[0] - expected) < 1e-12);
    std::fill(demand.begin(), demand.end(), .4);
    const auto half = detail::macro_routing_pressure(d, 64, 64, demand, 2);
    check(std::abs(half.x[0] * 4 - f.x[0]) < 1e-12);
    std::fill(demand.begin(), demand.end(), 0);
    check(detail::macro_routing_pressure(d, 64, 64, demand, 2).channels == 0);
    std::fill(demand.begin(), demand.end(), .8);
    auto fixed = d;
    fixed.num_movable = 1;
    const auto one = detail::macro_routing_pressure(fixed, 64, 64, demand, 2);
    check(one.x[0] == f.x[0] && one.x[1] == 0);
    fixed.num_movable = 0;
    check(detail::macro_routing_pressure(fixed, 64, 64, demand, 2).channels ==
          0);
    auto far = d;
    far.node_x[1] += 40;
    check(detail::macro_routing_pressure(far, 64, 64, demand, 2).channels == 0);
    auto vertical = d;
    std::swap(vertical.node_x, vertical.node_y);
    std::swap(vertical.node_size_x, vertical.node_size_y);
    const auto vf = detail::macro_routing_pressure(vertical, 64, 64, demand, 2);
    check(vf.y[0] == f.x[0] && vf.x[0] == 0);
    auto padded = d;
    padded.node_pad_x = padded.node_pad_y = {5, 5, 0, 0};
    for (int i = 0; i < 2; ++i) {
      padded.node_x[i] -= 5;
      padded.node_y[i] -= 5;
      padded.node_size_x[i] += 10;
      padded.node_size_y[i] += 10;
    }
    const auto pf = detail::macro_routing_pressure(padded, 64, 64, demand, 2);
    check(pf.x == f.x && pf.y == f.y);
    auto occluded = d;
    occluded.node_x.push_back(482);
    occluded.node_y.push_back(200);
    occluded.node_size_x.push_back(6);
    occluded.node_size_y.push_back(600);
    occluded.node_macro.push_back(1);
    occluded.build_maps();
    const auto blocked =
        detail::macro_routing_pressure(occluded, 64, 64, demand, 2);
    check(blocked.channels == 2); // outer macros see only the intervening macro
    // A real GP run with no movable standard cells proves this is direct macro
    // feedback, not a side effect of standard-cell inflation. RUDY is below 1.
    Options o;
    o.bins_x = o.bins_y = 64;
    o.max_iterations = o.min_iterations = 100;
    o.routability = o.improved_routability = true;
    o.route_capacity_h = o.route_capacity_v = .25;
    o.density_upper_bound = o.scale_macro_density = o.exact_macro_density =
        true;
    o.add_fillers = o.legalize = o.detailed_place = false;
    o.target_density = .8;
    auto baseline = d, pushed = d;
    const auto before = rudy(d, {64, 64, .25, .25, true, true});
    check(before.max_utilization > 0 && before.max_utilization < 1);
    place(baseline, o);
    o.routability_macro_weight = 1;
    const auto history = place(pushed, o);
    const double gap0 = baseline.node_x[1] - baseline.node_x[0] - 80;
    const double gap1 = pushed.node_x[1] - pushed.node_x[0] - 80;
    check(gap1 > gap0 + 1);
    check(pushed.node_x[2] == d.node_x[2] && pushed.node_y[3] == d.node_y[3]);
    check(std::any_of(
        history.history.begin(), history.history.end(), [](const auto &s) {
          return s.macro_routing_force > 0 && s.macro_routing_channels > 0;
        }));
    const auto after = rudy(pushed, {64, 64, .25, .25, true, true});
    check(after.max_utilization < before.max_utilization);
    std::cout << "gap " << gap0 << " -> " << gap1 << "; RUDY peak "
              << before.max_utilization << " -> " << after.max_utilization
              << '\n';
    for (int i = 0; i < 2; ++i)
      std::cout << "macro " << i << " control " << baseline.node_x[i] << ' '
                << baseline.node_y[i] << " pushed " << pushed.node_x[i] << ' '
                << pushed.node_y[i] << '\n';
    auto halo_gp = d;
    halo_gp.node_halo_x = halo_gp.node_halo_y = {5, 5, 0, 0};
    auto legal_options = o;
    legal_options.legalize = true;
    const auto legal = place(halo_gp, legal_options);
    check(legal.legal && is_legal(halo_gp));
    check(halo_gp.node_x[1] - halo_gp.node_x[0] - 80 >= 10);
    auto locked = d;
    locked.num_movable = 1;
    place(locked, o);
    check(locked.node_x[1] == d.node_x[1] && locked.node_y[1] == d.node_y[1]);
    // Connected macros have a competing wirelength pull toward the channel.
    // This is a trade-off, not an instruction to scatter every SRAM.
    auto connected = d;
    connected.node_x.push_back(485);
    connected.node_y.push_back(500);
    connected.node_size_x.push_back(0);
    connected.node_size_y.push_back(0);
    connected.node_macro.push_back(0);
    connected.pin2node = {2, 3, 0, 4, 1, 4};
    connected.pin_offset_x = {0, 0, 80, 0, 0, 0};
    connected.pin_offset_y = {0, 0, 300, 0, 300, 0};
    connected.flat_net2pin = {0, 1, 2, 3, 4, 5};
    connected.flat_net2pin_start = {0, 2, 4, 6};
    connected.net_weights = {1, .1, .1};
    connected.net_mask.clear();
    connected.build_maps();
    connected.node_halo_x = connected.node_halo_y = {5, 5, 0, 0, 0};
    auto linked_control = connected, linked_push = connected;
    auto linked_options = o;
    linked_options.macro_overlap_weight = 50;
    linked_options.legalize = true;
    linked_options.routability_macro_weight = 0;
    place(linked_control, linked_options);
    linked_options.routability_macro_weight = 1;
    place(linked_push, linked_options);
    const double linked_gap0 =
        linked_control.node_x[1] - linked_control.node_x[0] - 80;
    const double linked_gap1 =
        linked_push.node_x[1] - linked_push.node_x[0] - 80;
    const auto route0 = rudy(linked_control, {64, 64, .25, .25, true, true});
    const auto route1 = rudy(linked_push, {64, 64, .25, .25, true, true});
    std::cout << "connected gap " << linked_gap0 << " -> " << linked_gap1
              << "; RUDY " << route0.max_utilization << " -> "
              << route1.max_utilization << "; HPWL " << hpwl(linked_control)
              << " -> " << hpwl(linked_push) << '\n';
    check(is_legal(linked_control) && is_legal(linked_push));
    check(linked_gap0 >= 10 && linked_gap1 > linked_gap0);
    for (int i = 0; i < 2; ++i)
      std::cout << "connected macro " << i << " control "
                << linked_control.node_x[i] << ' ' << linked_control.node_y[i]
                << " pushed " << linked_push.node_x[i] << ' '
                << linked_push.node_y[i] << '\n';
    check(route1.max_utilization < route0.max_utilization);
    // Full continuation retains the RUDY map between refreshes.
    auto full = d, split = d;
    GpState state;
    o.max_iterations = 60;
    place(full, o);
    o.max_iterations = 35;
    place(split, o, &state);
    o.max_iterations = 25;
    place(split, o, &state);
    check(full.node_x == split.node_x && full.node_y == split.node_y);
    o.routability_macro_weight = 2;
    bool rejected = false;
    try {
      place(split, o, &state);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected);
    auto invalid = o;
    invalid.routability_macro_range = 0;
    rejected = false;
    try {
      place(d, invalid);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    check(rejected);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
