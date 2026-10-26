// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace DreamPlace::native {
namespace {
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
} // namespace
ReservationCharge reservation_charge(const Design &d,
                                     const std::vector<double> &inflation) {
  ReservationCharge result;
  std::set<std::string> owners;
  std::set<int> nodes;
  for (const auto &r : d.density_reservations) {
    require(!r.owner.empty() && owners.insert(r.owner).second,
            "duplicate/empty density reservation owner");
    require(std::isfinite(r.x) && std::isfinite(r.y) &&
                std::isfinite(r.width) && std::isfinite(r.height) &&
                r.width > 0 && r.height > 0 &&
                std::isfinite(r.width * r.height) && r.x >= d.layout.xl &&
                r.y >= d.layout.yl && r.x + r.width <= d.layout.xh &&
                r.y + r.height <= d.layout.yh,
            "invalid density reservation rectangle");
    require(r.materialized_node >= -1, "invalid materialized reservation node");
    if (r.materialized_node < 0) {
      result.virtual_area += r.width * r.height;
      ++result.virtual_count;
    } else {
      const size_t n = r.materialized_node;
      require(d.num_movable >= 0 && n < size_t(d.num_movable) &&
                  n < d.node_x.size() && n < d.node_y.size() &&
                  n < d.node_size_x.size() && n < d.node_size_y.size() &&
                  (d.node_macro.empty() || n < d.node_macro.size()) &&
                  nodes.insert(n).second &&
                  (d.node_macro.empty() || !d.node_macro[n]) &&
                  d.node_size_x[n] > 0 && d.node_size_y[n] > 0,
              "stale/duplicate/non-cell materialized reservation binding");
      require(inflation.empty() ||
                  (n < inflation.size() && std::isfinite(inflation[n]) &&
                   inflation[n] >= 1),
              "invalid reservation inflation");
      result.materialized_area += d.node_size_x[n] * d.node_size_y[n] *
                                  (inflation.empty() ? 1 : inflation[n]);
      result.predicted_materialized_area += r.width * r.height;
      ++result.materialized_count;
    }
  }
  // Accepted active x intervals must be disjoint in y. Check nearest y
  // neighbors while expiring rectangles by x, without a quadratic all-pairs
  // scan on long columns of reserved row sites. Touching edges stay legal.
  std::vector<const DensityReservation *> rectangles;
  for (const auto &r : d.density_reservations)
    if (r.materialized_node < 0)
      rectangles.push_back(&r);
  std::sort(rectangles.begin(), rectangles.end(), [](auto a, auto b) {
    return std::tie(a->x, a->y, a->owner) < std::tie(b->x, b->y, b->owner);
  });
  std::map<double, std::pair<double, size_t>> active;
  using Expiry = std::tuple<double, size_t, double>;
  std::priority_queue<Expiry, std::vector<Expiry>, std::greater<Expiry>> expiry;
  for (size_t i = 0; i < rectangles.size(); ++i) {
    const auto &r = *rectangles[i];
    while (!expiry.empty() && std::get<0>(expiry.top()) <= r.x) {
      auto [end, id, y] = expiry.top();
      expiry.pop();
      (void)end;
      auto old = active.find(y);
      if (old != active.end() && old->second.second == id)
        active.erase(old);
    }
    auto next = active.lower_bound(r.y);
    require(next == active.end() || next->first >= r.y + r.height,
            "overlapping virtual density reservations");
    require(next == active.begin() || std::prev(next)->second.first <= r.y,
            "overlapping virtual density reservations");
    active.emplace(r.y, std::make_pair(r.y + r.height, i));
    expiry.emplace(r.x + r.width, i, r.y);
  }
  require(std::isfinite(result.virtual_area) &&
              std::isfinite(result.materialized_area) &&
              std::isfinite(result.predicted_materialized_area),
          "density reservation area overflow");
  return result;
}
namespace detail {
void reservation_density_map(const Design &d, int bx, int by, bool exact,
                             int threads, std::vector<double> &map) {
  require(bx > 0 && by > 0 && bx <= 4096 && by <= 4096 && threads > 0,
          "invalid reservation density grid");
  reservation_charge(d);
  Design charge;
  charge.layout = d.layout;
  const double dx = (d.layout.xh - d.layout.xl) / bx;
  const double dy = (d.layout.yh - d.layout.yl) / by;
  std::vector<double> sx, sy, ox, oy, ratio;
  for (const auto &r : d.density_reservations) {
    if (r.materialized_node >= 0)
      continue;
    const double w = exact ? r.width : std::max(r.width, std::sqrt(2.) * dx);
    const double h = exact ? r.height : std::max(r.height, std::sqrt(2.) * dy);
    charge.node_x.push_back(r.x);
    charge.node_y.push_back(r.y);
    sx.push_back(w);
    sy.push_back(h);
    ox.push_back((r.width - w) / 2);
    oy.push_back((r.height - h) / 2);
    ratio.push_back(r.width * r.height / (w * h));
  }
  density_map(charge, bx, by, threads, sx, sy, ox, oy, ratio, map);
}
double reservation_filler_area(double budget, double rounding,
                               const ReservationCharge &charge) {
  require(std::isfinite(budget) && budget >= 0 && std::isfinite(rounding) &&
              rounding >= 0,
          "invalid unreserved filler budget");
  const double available = std::max(0., budget - rounding);
  const double required = charge.virtual_area + charge.materialized_area;
  require(std::isfinite(required) && required >= 0 &&
              required <= available + 1e-12 * std::max(1., available),
          "insufficient filler charge for density reservations");
  return std::max(0., available - required);
}
} // namespace detail
} // namespace DreamPlace::native
