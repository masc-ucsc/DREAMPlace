// BSD-3-Clause. See LICENSE. Assertions remain active under -c opt.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace DreamPlace::native;
namespace {
void check(bool ok, const char *why) {
  if (!ok)
    throw std::runtime_error(why);
}
void near(double a, double b, double eps, const char *why) {
  check(std::abs(a - b) <= eps, why);
}
template <class F> void rejects(F fn) {
  bool failed = false;
  try {
    fn();
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  check(failed, "invalid reservation was accepted");
}
Design network() {
  Design d;
  d.layout = {0, 0, 100, 8, 1, 1};
  d.num_movable = 2;
  d.node_x = {20, 70, 0, 98};
  d.node_y = {2, 2, 2, 2};
  d.node_size_x = {2, 2, 0, 0};
  d.node_size_y = {1, 1, 0, 0};
  d.node_macro = {0, 0, 0, 0};
  d.pin2node = {0, 1, 2, 3};
  d.pin_offset_x = {1, 1, 0, 0};
  d.pin_offset_y = {.5, .5, 0, 0};
  d.flat_net2pin = {0, 1, 2, 3};
  d.flat_net2pin_start = {0, 4};
  d.build_maps();
  d.validate();
  return d;
}
Options controls(int iterations = 20) {
  Options o;
  o.max_iterations = o.min_iterations = iterations;
  o.bins_x = 16;
  o.bins_y = 4;
  o.target_density = .5;
  o.legalize = false;
  o.detailed_place = false;
  o.num_threads = 1;
  o.filler_width = 1;
  o.wirelength_gamma_start = o.wirelength_gamma_end = 1;
  return o;
}
DensityReservation footprint() { return {"net0/buffer0", 40, 2, 2, 2, -1}; }
void rectangle_and_identity_validation() {
  auto base = network();
  base.density_reservations = {footprint()};
  base.validate();
  for (int i = 0; i < 10; ++i) {
    auto d = base;
    auto &r = d.density_reservations[0];
    if (i == 0)
      r.owner.clear();
    if (i == 1)
      d.density_reservations.push_back(r);
    if (i == 2)
      r.width = 0;
    if (i == 3)
      r.height = -1;
    if (i == 4)
      r.x = std::numeric_limits<double>::quiet_NaN();
    if (i == 5)
      r.x = 99;
    if (i == 6)
      r.materialized_node = -2;
    if (i == 7)
      r.materialized_node = 100;
    if (i == 8)
      r.materialized_node = 2;
    if (i == 9) {
      r.materialized_node = 0;
      d.node_macro[0] = 1;
    }
    rejects([&] { d.validate(); });
  }
  auto d = base;
  d.density_reservations[0].materialized_node = 0;
  auto r = footprint();
  r.owner = "other";
  r.materialized_node = 0;
  d.density_reservations.push_back(r);
  rejects([&] { d.validate(); });
}
void direct_materialized_binding_validation() {
  for (int i = 0; i < 7; ++i) {
    auto d = network();
    d.density_reservations = {footprint()};
    d.density_reservations[0].materialized_node = 1;
    if (i == 0)
      d.node_x.clear();
    if (i == 1)
      d.node_y.clear();
    if (i == 2)
      d.node_size_x.clear();
    if (i == 3)
      d.node_size_y.clear();
    if (i == 4)
      d.node_macro.resize(1);
    if (i == 5)
      d.num_movable = -1;
    if (i == 6) {
      d.num_movable = 100;
      d.density_reservations[0].materialized_node = 99;
    }
    rejects([&] { reservation_charge(d); });
  }
}
void overlap_and_touching() {
  auto d = network();
  d.density_reservations = {footprint()};
  auto other = footprint();
  other.owner = "second";
  d.density_reservations.push_back(other);
  rejects([&] { d.validate(); });
  d.density_reservations[1].x += 2;
  d.validate();
  d.density_reservations[1].x -= 1;
  d.density_reservations[1].y += 2;
  d.validate();
}
void native_map_oracle() {
  auto d = network();
  d.density_reservations = {footprint()};
  for (bool exact : {true, false}) {
    std::vector<double> actual;
    detail::reservation_density_map(d, 16, 4, exact, 1, actual);
    Design oracle;
    oracle.layout = d.layout;
    oracle.node_x = {40};
    oracle.node_y = {2};
    const double sx = exact ? 2 : std::sqrt(2.) * 100 / 16;
    const double sy = exact ? 2 : std::sqrt(2.) * 8 / 4;
    std::vector<double> expected;
    detail::density_map(oracle, 16, 4, 1, {sx}, {sy}, {(2 - sx) / 2},
                        {(2 - sy) / 2}, {4 / (sx * sy)}, expected);
    check(actual == expected,
          "reservation is not ordinary cell density charge");
    near(std::accumulate(actual.begin(), actual.end(), 0.), 4, 1e-12,
         "reservation area lost");
  }
  std::vector<double> actual;
  detail::reservation_density_map(d, 16, 4, true, 1, actual);
  // Independent rectangle-bin intersection, with unchanged target capacity.
  for (int x = 0; x < 16; ++x)
    for (int y = 0; y < 4; ++y) {
      double w =
          std::max(0., std::min(42., (x + 1) * 6.25) - std::max(40., x * 6.25));
      double h =
          std::max(0., std::min(4., (y + 1) * 2.) - std::max(2., y * 2.));
      near(actual[x * 4 + y], w * h, 1e-12, "rectangle oracle disagrees");
    }
  d.density_reservations[0].materialized_node = 0;
  detail::reservation_density_map(d, 16, 4, true, 1, actual);
  near(std::accumulate(actual.begin(), actual.end(), 0.), 0, 0,
       "physical buffer counted twice");
}
void charge_transfer_and_budget() {
  for (double fraction : {1., .25}) {
    const double budget = 100 * fraction, rounding = .25;
    auto d = network();
    d.density_reservations = {footprint()};
    const auto v = reservation_charge(d);
    near(v.virtual_area, 4, 0, "wrong virtual area");
    const double vf = detail::reservation_filler_area(budget, rounding, v);
    near(vf + v.virtual_area, budget - rounding, 1e-12,
         "reserve charge not conserved");
    d.density_reservations[0].materialized_node = 0;
    d.node_size_x[0] = 2;
    d.node_size_y[0] = 2;
    const auto p = reservation_charge(d);
    near(p.virtual_area, 0, 0, "activation retained virtual charge");
    near(p.materialized_area, 4, 0, "activation lost physical charge");
    const double pf = detail::reservation_filler_area(budget, rounding, p);
    near(pf, vf, 0, "reserve/activate changed filler charge");
    near(p.materialized_area + pf, budget - rounding, 1e-12,
         "activation double counted charge");
    d.density_reservations.clear();
    near(detail::reservation_filler_area(budget, rounding,
                                         reservation_charge(d)),
         budget - rounding, 0, "remove did not release filler charge");
    rejects([&] { detail::reservation_filler_area(1, 0, v); });
  }
}
void actual_gp_density_and_clear() {
  auto a = network(), b = a;
  auto o = controls();
  b.density_reservations = {footprint()};
  auto control = place(a, o), reserved = place(b, o);
  check(a.node_x != b.node_x || a.node_y != b.node_y,
        "reservation did not enter actual GP forces");
  near(reserved.filler_area - reserved.effective_filler_area, 4, 1e-10,
       "actual GP did not displace equal filler charge");
  near(control.filler_area, reserved.filler_area, 0,
       "reservation resized filler particles/count");
  near(reserved.reserved_area, 4, 0, "final reservation metric lost");
  check(a.pin2node == b.pin2node && a.node_x.size() == b.node_x.size(),
        "reservation entered physical maps");
  auto c = network();
  o.checkpoint_interval = 5;
  o.checkpoint = [](const Design &, int i) {
    GpUpdate u;
    if (i == 5)
      u.density_reservations = std::vector<DensityReservation>{footprint()};
    if (i == 10)
      u.density_reservations = std::vector<DensityReservation>{};
    return u;
  };
  auto result = place(c, o);
  check(c.density_reservations.empty(), "clear left obsolete reservation");
  near(result.effective_filler_area, result.filler_area, 1e-10,
       "clear did not release filler charge");
}
void continuation_and_atomicity() {
  auto d = network();
  d.density_reservations = {footprint()};
  auto o = controls(10);
  GpState state;
  place(d, o, &state);
  auto changed = d;
  changed.density_reservations.clear();
  auto old = state;
  rejects([&] { place(changed, o, &state); });
  check(state.x == old.x && changed.node_x == d.node_x,
        "failed continuation mutated state/design");
  auto uninterrupted = d;
  auto copied = state;
  auto a = place(d, o, &state), b = place(uninterrupted, o, &copied);
  check(d.node_x == uninterrupted.node_x && state.x == copied.x,
        "resumed reservation trajectory is not deterministic");
  near(a.effective_filler_area, b.effective_filler_area, 0,
       "continuation filler charge drift");
  auto bad = network(), before = bad;
  o.checkpoint_interval = 5;
  o.checkpoint = [](const Design &, int) {
    GpUpdate u;
    auto r = footprint();
    r.width = 99;
    u.density_reservations = std::vector<DensityReservation>{r};
    return u;
  };
  rejects([&] { place(bad, o); });
  check(bad.node_x == before.node_x && bad.density_reservations.empty(),
        "failed checkpoint partially committed");
}
void rollback_trajectory() {
  auto base = network(), trial = base, control = base;
  auto o = controls(30);
  o.checkpoint_interval = 5;
  o.checkpoint_rollback = true;
  o.checkpoint = [](const Design &, int i) {
    GpUpdate u;
    if (i == 5)
      u.action = GpCheckpointAction::retain;
    if (i == 10)
      u.density_reservations = std::vector<DensityReservation>{footprint()};
    if (i == 20)
      u.action = GpCheckpointAction::restore;
    return u;
  };
  auto result = place(trial, o);
  auto reference = controls(15);
  reference.checkpoint_interval = 5;
  reference.checkpoint_rollback = true;
  reference.checkpoint = [](const Design &, int i) {
    GpUpdate u;
    if (i == 5)
      u.action = GpCheckpointAction::retain;
    return u;
  };
  place(control, reference);
  check(trial.node_x == control.node_x && trial.node_y == control.node_y,
        "rollback did not restore exact reservation/filler trajectory");
  check(trial.density_reservations.empty() && result.iterations == 30,
        "rollback rewound budget or retained reservation");
  near(result.effective_filler_area, result.filler_area, 1e-10,
       "rollback retained spent filler charge");
}
void final_checkpoint_and_overflow() {
  auto metric = network();
  metric.density_reservations = {{"wide", 40, 0, 10, 6, -1}};
  near(density_overflow(metric, 16, 4, .5), 0, 0, "physical overflow changed");
  check(reservation_density_overflow(metric, 16, 4, .5) > 0,
        "solver overflow omitted virtual charge");

  auto d = network();
  auto o = controls(5);
  o.checkpoint_interval = 5;
  o.checkpoint = [](const Design &, int) {
    GpUpdate u;
    u.density_reservations = std::vector<DensityReservation>{footprint()};
    return u;
  };
  auto r = place(d, o);
  near(r.filler_area - r.effective_filler_area, 4, 1e-10,
       "final update reports stale filler charge");
  near(r.reserved_area, 4, 0, "final update lost reservation state");
  auto physical = d;
  physical.density_reservations.clear();
  near(density_overflow(d, 16, 4, .5), density_overflow(physical, 16, 4, .5), 0,
       "physical overflow silently includes virtual area");
  auto before = d;
  auto bad = controls(5);
  bad.filler_area_fraction = .001;
  rejects([&] { place(d, bad); });
  check(d.node_x == before.node_x &&
            d.density_reservations == before.density_reservations,
        "charge exhaustion partially committed");
}
void topology_activation_removal() {
  for (double fraction : {1., .25}) {
    auto d = network();
    d.density_reservations = {footprint()};
    auto o = controls(60);
    o.checkpoint_interval = 20;
    o.filler_area_fraction = fraction;
    o.checkpoint = [](const Design &snapshot, int i) {
      GpUpdate u;
      if (i != 20 && i != 40)
        return u;
      auto c = snapshot;
      if (i == 20) {
        // A physical occupancy fixture. Logical buffer pool/net splitting is a
        // separate host transaction; this verifies native density transfer.
        for (auto *v : {&c.node_x, &c.node_y, &c.node_size_x, &c.node_size_y})
          v->insert(v->begin() + 2, 0);
        c.node_x[2] = 40;
        c.node_y[2] = 2;
        c.node_size_x[2] = 2;
        c.node_size_y[2] = 2;
        c.node_macro.insert(c.node_macro.begin() + 2, 0);
        c.num_movable = 3;
        for (int &n : c.pin2node)
          if (n >= 2)
            ++n;
        c.density_reservations[0].materialized_node = 2;
        u.previous_node = {0, 1, -1, 2, 3};
      } else {
        for (auto *v : {&c.node_x, &c.node_y, &c.node_size_x, &c.node_size_y})
          v->erase(v->begin() + 2);
        c.node_macro.erase(c.node_macro.begin() + 2);
        c.num_movable = 2;
        for (int &n : c.pin2node)
          if (n > 2)
            --n;
        c.density_reservations.clear();
        u.previous_node = {0, 1, 3, 4};
      }
      c.build_maps();
      u.replacement = c;
      return u;
    };
    auto result = place(d, o);
    check(d.node_x.size() == 4 && d.density_reservations.empty(),
          "topology cycle did not remove buffer/reservation");
    near(result.effective_filler_area, result.filler_area, 1e-9,
         "topology removal did not restore filler charge");
    check(result.iterations == 60, "topology changes reset spent budget");
  }
}
} // namespace
int main() {
  try {
    rectangle_and_identity_validation();
    direct_materialized_binding_validation();
    overlap_and_touching();
    native_map_oracle();
    charge_transfer_and_budget();
    actual_gp_density_and_clear();
    continuation_and_atomicity();
    rollback_trajectory();
    topology_activation_removal();
    final_checkpoint_and_overflow();
    std::cout << "10 reservation density/charge/GP/checkpoint checks passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
