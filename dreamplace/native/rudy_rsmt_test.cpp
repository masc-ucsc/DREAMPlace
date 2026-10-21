// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy_tree.hpp"
#include "thirdparty/flute/flute.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <source_location>
#include <stdexcept>
using namespace DreamPlace::native;
void check(bool ok, std::source_location at = std::source_location::current()) {
  if (!ok)
    throw std::runtime_error("RSMT route check failed at line " +
                             std::to_string(at.line()));
}
void near(double a, double b,
          std::source_location at = std::source_location::current()) {
  check(std::abs(a - b) < 1e-9 * std::max({1., std::abs(a), std::abs(b)}), at);
}
template <class F> void rejects(F f) {
  bool caught = false;
  try {
    f();
  } catch (const std::exception &) {
    caught = true;
  }
  check(caught);
}
Design make(const std::vector<std::pair<double, double>> &points) {
  Design d;
  d.layout = {0, 0, 100, 100, 1, 1};
  for (auto [x, y] : points) {
    d.node_x.push_back(x);
    d.node_y.push_back(y);
    d.pin2node.push_back(d.pin2node.size());
  }
  d.node_size_x.resize(points.size());
  d.node_size_y.resize(points.size());
  d.pin_offset_x.resize(points.size());
  d.pin_offset_y.resize(points.size());
  d.flat_net2pin = d.pin2node;
  d.flat_net2pin_start = {0, int(points.size())};
  d.net_weights = {1};
  d.build_maps();
  return d;
}
int main(int argc, char **argv) {
  try {
    check(argc == 3);
    RudyOptions o;
    o.model = RudyModel::rsmt;
    o.bins_x = o.bins_y = 10;
    o.initialize = [&] {
      static std::once_flag once;
      std::call_once(once, [&] { flute::readLUT(argv[1], argv[2]); });
    };
    auto d = make({{10, 10}, {20, 10}, {10, 20}, {20, 20}});
    auto a = rudy(d, o), b = rudy(d, o);
    check(a.horizontal == b.horizontal && a.vertical == b.vertical);
    // A four-corner Steiner tree has length 3L, versus HPWL 2L: a 1.12
    // universal HPWL bound would be mathematically incorrect.
    const double scale = 100.;
    near((std::accumulate(a.horizontal.begin(), a.horizontal.end(), 0.) +
          std::accumulate(a.vertical.begin(), a.vertical.end(), 0.)) *
             scale,
         30.);
    RudyOptions missing = o;
    missing.initialize = {};
    rejects([&] { rudy(d, missing); });
    o.degree = 3;
    auto fallback = rudy(d, o);
    auto bbox = o;
    bbox.model = RudyModel::bbox;
    auto legacy = rudy(d, bbox);
    for (size_t k = 0; k < legacy.horizontal.size(); ++k) {
      near(fallback.horizontal[k], legacy.horizontal[k]);
      near(fallback.vertical[k], legacy.vertical[k]);
    }
    o.degree = 40;
    for (bool obstacles : {false, true}) {
      o.macro_blockages = obstacles;
      auto two = make({{15, 15}, {85, 85}});
      auto first = rudy(two, o);
      bbox = o;
      bbox.model = RudyModel::bbox;
      auto second = rudy(two, bbox);
      for (size_t k = 0; k < first.horizontal.size(); ++k) {
        near(first.horizontal[k], second.horizontal[k]);
        near(first.vertical[k], second.vertical[k]);
      }
      auto straight = rudy(make({{15, 15}, {85, 15}}), o);
      check(straight.measured_nets == 1);
      check(std::accumulate(straight.horizontal.begin(),
                            straight.horizontal.end(), 0.) > 0);
      auto zero = rudy(make({{15, 15}, {15, 15}}), o);
      check(zero.degenerate_nets == 1);
      near(zero.overflow, 0);
      auto dup = rudy(make({{15, 15}, {15, 15}, {85, 85}, {85, 85}}), o);
      check(dup.measured_nets == 1);
    }
    o.macro_blockages = false;
    rejects([&] { detail::route_tree({{1e20, 0}, {0, 0}}, o); });
    std::mt19937 random(1);
    for (int degree = 3; degree <= 15; ++degree)
      for (int trial = 0; trial < 20; ++trial) {
        std::vector<std::pair<double, double>> points;
        double xl = 100, yl = 100, xh = 0, yh = 0;
        for (int i = 0; i < degree; ++i) {
          double x = 1 + random() % 98, y = 1 + random() % 98;
          points.emplace_back(x, y);
          xl = std::min(xl, x);
          xh = std::max(xh, x);
          yl = std::min(yl, y);
          yh = std::max(yh, y);
        }
        double length = 0;
        for (auto e : detail::route_tree(points, o))
          length += e.xh - e.xl + e.yh - e.yl;
        check(length + 1e-9 >= xh - xl + yh - yl);
        // HPWL is a lower bound, not an upper bound within 12%. Check
        // against an independently computed Manhattan spanning tree instead.
        std::vector<double> distance(points.size(), 1e100);
        std::vector<bool> used(points.size());
        distance[0] = 0;
        double mst = 0;
        for (size_t step = 0; step < points.size(); ++step) {
          size_t best = points.size();
          for (size_t i = 0; i < points.size(); ++i)
            if (!used[i] &&
                (best == points.size() || distance[i] < distance[best]))
              best = i;
          used[best] = true;
          mst += distance[best];
          for (size_t i = 0; i < points.size(); ++i)
            distance[i] =
                std::min(distance[i],
                         std::abs(points[i].first - points[best].first) +
                             std::abs(points[i].second - points[best].second));
        }
        check(length <= 1.12 * mst + 1e-9);
      }
    auto gp = make({{10, 10}, {80, 80}});
    gp.num_movable = 2;
    gp.node_size_x = {1, 1};
    gp.node_size_y = {1, 1};
    Options opts;
    opts.routability = true;
    opts.legalize = false;
    opts.detailed_place = false;
    opts.max_iterations = 60;
    opts.bins_x = opts.bins_y = 8;
    opts.route_capacity_h = opts.route_capacity_v = 100;
    opts.route.model = RouteModel::automatic;
    opts.route.hold = 2;
    auto placed = place(gp, opts);
    check(placed.history.front().route_peak_utilization > 0);
    check(placed.route_model == RouteModel::rsmt);
    check(placed.route_switch_iteration == 21);
    check(placed.route_switch_rudy == 0 && placed.route_switch_rsmt == 0);
    bool seen = false;
    for (auto sample : placed.history) {
      if (sample.route_model == RouteModel::rsmt)
        seen = true;
      if (seen)
        check(sample.route_model == RouteModel::rsmt);
    }
    // A continuation retains the one-way latch and the hold count.
    GpState continuation;
    auto continued = make({{10, 10}, {80, 80}});
    continued.num_movable = 2;
    continued.node_size_x = {1, 1};
    continued.node_size_y = {1, 1};
    auto continuation_opts = opts;
    continuation_opts.max_iterations = 20;
    auto first_stage = place(continued, continuation_opts, &continuation);
    check(first_stage.route_model == RouteModel::rudy);
    continuation_opts.max_iterations = 40;
    auto resumed = place(continued, continuation_opts, &continuation);
    check(resumed.route_model == RouteModel::rsmt);
    check(resumed.route_switch_iteration == 21);
    continuation_opts.route.hold = 3;
    rejects([&] { place(continued, continuation_opts, &continuation); });
    // A diagnostic probe also runs with the legacy density-only stop mode.
    auto probe_design = make({{10, 10}, {90, 10}, {10, 90}, {90, 90}});
    probe_design.num_movable = 1;
    probe_design.node_size_x[0] = probe_design.node_size_y[0] = 1;
    auto probe_opts = opts;
    probe_opts.route_initialize = o.initialize;
    probe_opts.route.model = RouteModel::rudy;
    probe_opts.route.probe = true;
    probe_opts.route.probe_iterations = 20;
    probe_opts.route.target_overflow = 0;
    probe_opts.route_capacity_h = probe_opts.route_capacity_v = .1;
    probe_opts.post_density_factor = 0;
    probe_opts.add_fillers = false;
    probe_opts.max_iterations = 100;
    probe_opts.min_iterations = 20;
    auto probed = place(probe_design, probe_opts);
    check(probed.route_probed);
    check(probed.route_switch_rsmt >= 0);
    if (probed.route_switch_rsmt > 0) {
      check(probed.route_model == RouteModel::rsmt);
      check(probed.iterations <=
            probed.route_switch_iteration + probe_opts.route.probe_iterations);
    }
    // A passing probe accepts the candidate even when WDM-weighted RUDY
    // is stricter than the tree model on duplicated terminal locations.
    std::vector<std::pair<double, double>> repeated;
    for (int i = 0; i < 10; ++i)
      repeated.emplace_back(i % 2 ? 90 : 10, i % 2 ? 90 : 10);
    auto passed = make(repeated);
    passed.node_x.insert(passed.node_x.begin(), 50);
    passed.node_y.insert(passed.node_y.begin(), 50);
    passed.node_size_x.insert(passed.node_size_x.begin(), 1);
    passed.node_size_y.insert(passed.node_size_y.begin(), 1);
    for (auto &node : passed.pin2node)
      ++node;
    passed.num_movable = 1;
    passed.build_maps();
    probe_opts.post_density_factor = 4;
    probe_opts.route_capacity_h = probe_opts.route_capacity_v = .016;
    auto accepted = place(passed, probe_opts);
    check(accepted.route_probed);
    check(accepted.route_switch_rsmt == 0);
    check(accepted.route_switch_rudy > 0);
    check(accepted.converged);
    check(accepted.route_model == RouteModel::rudy);
    opts.route.model = RouteModel::off;
    opts.route_capacity_h = opts.route_capacity_v = 0;
    placed = place(gp, opts);
    check(placed.rudy_calls == 0);
    opts.route.target_overflow = .1;
    opts.route.switch_overflow = .01;
    rejects([&] { place(gp, opts); });
    std::cout << "RSMT and route policy checks passed\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
