// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/magnet.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>

namespace DreamPlace::native {
namespace {
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
double pad(const Design &d, int i, bool x) {
  const auto &internal = x ? d.node_pad_x : d.node_pad_y;
  if (!internal.empty())
    return internal[i];
  const auto &halo = x ? d.node_halo_x : d.node_halo_y;
  const double pitch = x ? d.layout.site_width : d.layout.row_height;
  return halo.empty() ? 0 : std::ceil(halo[i] / pitch - 1e-9) * pitch;
}
MagnetBox body(const Design &d, int i) {
  const double hx = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
  const double hy = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
  return {d.node_x[i] + hx, d.node_y[i] + hy,
          d.node_x[i] + d.node_size_x[i] - hx,
          d.node_y[i] + d.node_size_y[i] - hy};
}
MagnetBox outer(const Design &d, int i) {
  const auto b = body(d, i);
  return {b.xl - pad(d, i, true), b.yl - pad(d, i, false),
          b.xh + pad(d, i, true), b.yh + pad(d, i, false)};
}
double overlap_area(const MagnetBox &a, const MagnetBox &b) {
  return std::max(0., std::min(a.xh, b.xh) - std::max(a.xl, b.xl)) *
         std::max(0., std::min(a.yh, b.yh) - std::max(a.yl, b.yl));
}
bool cell(const Design &d, int i) {
  return i < d.num_movable && (d.node_macro.empty() || !d.node_macro[i]);
}
MagnetBox bounds(const MagnetRegion &region) {
  MagnetBox box{region[0].first, region[0].second, region[0].first,
                region[0].second};
  for (const auto &[x, y] : region) {
    box.xl = std::min(box.xl, x);
    box.xh = std::max(box.xh, x);
    box.yl = std::min(box.yl, y);
    box.yh = std::max(box.yh, y);
  }
  return box;
}
// Sutherland-Hodgman clipping, with stack storage for a clipped quadrilateral.
// Translate the shoelace sum to a vertex to avoid cancellation at large
// origins.
double overlap_area(const MagnetRegion &region, const MagnetBox &box) {
  if (box.xh <= box.xl || box.yh <= box.yl ||
      overlap_area(bounds(region), box) <= 0)
    return 0;
  std::array<std::pair<double, double>, 12> polygon{}, clipped{};
  std::copy(region.begin(), region.end(), polygon.begin());
  int count = 4;
  for (int side = 0; side < 4 && count; ++side) {
    const bool x = side < 2, lower = side % 2 == 0;
    const double edge =
        x ? (lower ? box.xl : box.xh) : (lower ? box.yl : box.yh);
    const auto coordinate = [&](const auto &point) {
      return x ? point.first : point.second;
    };
    const auto inside = [&](const auto &point) {
      return lower ? coordinate(point) >= edge : coordinate(point) <= edge;
    };
    int size = 0;
    auto previous = polygon[count - 1];
    bool was_inside = inside(previous);
    for (int i = 0; i < count; ++i) {
      const auto current = polygon[i];
      const bool is_inside = inside(current);
      if (is_inside != was_inside) {
        const double t = (edge - coordinate(previous)) /
                         (coordinate(current) - coordinate(previous));
        clipped[size++] = {
            previous.first + t * (current.first - previous.first),
            previous.second + t * (current.second - previous.second)};
      }
      if (is_inside)
        clipped[size++] = current;
      previous = current;
      was_inside = is_inside;
    }
    polygon.swap(clipped);
    count = size;
  }
  double twice_area = 0;
  for (int i = 1; i + 1 < count; ++i)
    twice_area += (polygon[i].first - polygon[0].first) *
                      (polygon[i + 1].second - polygon[0].second) -
                  (polygon[i].second - polygon[0].second) *
                      (polygon[i + 1].first - polygon[0].first);
  return std::abs(twice_area) / 2;
}
// Exact union of clipped physical macro bodies. Overlapping macros during GP
// must not subtract the same part of the core twice. Halos are not cell area.
double macro_union_area(const Design &d) {
  struct Edge {
    double x, yl, yh;
    int delta;
  };
  std::vector<Edge> edges;
  std::vector<double> ys;
  for (size_t i = 0; i < d.node_x.size(); ++i) {
    if (d.node_macro.empty() || !d.node_macro[i])
      continue;
    auto b = body(d, i);
    b.xl = std::max(b.xl, d.layout.xl);
    b.xh = std::min(b.xh, d.layout.xh);
    b.yl = std::max(b.yl, d.layout.yl);
    b.yh = std::min(b.yh, d.layout.yh);
    if (b.xl >= b.xh || b.yl >= b.yh)
      continue;
    edges.push_back({b.xl, b.yl, b.yh, 1});
    edges.push_back({b.xh, b.yl, b.yh, -1});
    ys.push_back(b.yl);
    ys.push_back(b.yh);
  }
  if (edges.empty())
    return 0;
  std::sort(ys.begin(), ys.end());
  ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
  std::sort(edges.begin(), edges.end(),
            [](const auto &a, const auto &b) { return a.x < b.x; });
  std::vector<int> cover(4 * ys.size());
  std::vector<double> length(cover.size());
  const auto update = [&](auto &&self, int node, int lo, int hi, int from,
                          int to, int delta) -> void {
    if (from <= lo && hi <= to)
      cover[node] += delta;
    else {
      const int mid = (lo + hi) / 2;
      if (from < mid)
        self(self, 2 * node, lo, mid, from, to, delta);
      if (to > mid)
        self(self, 2 * node + 1, mid, hi, from, to, delta);
    }
    length[node] =
        cover[node] > 0
            ? ys[hi] - ys[lo]
            : (hi - lo == 1 ? 0 : length[2 * node] + length[2 * node + 1]);
  };
  double area = 0, previous = edges.front().x;
  for (const auto &edge : edges) {
    area += (edge.x - previous) * length[1];
    const int lo = std::lower_bound(ys.begin(), ys.end(), edge.yl) - ys.begin();
    const int hi = std::lower_bound(ys.begin(), ys.end(), edge.yh) - ys.begin();
    update(update, 1, 0, int(ys.size()) - 1, lo, hi, edge.delta);
    previous = edge.x;
  }
  return area;
}
double exact_density(const Design &d, const MagnetRegion &region) {
  const double area = overlap_area(region, bounds(region));
  if (area <= 0)
    return 0;
  double occupied = 0;
  for (int i = 0; i < d.num_movable; ++i)
    if (cell(d, i))
      occupied += overlap_area(region, body(d, i));
  return occupied / area;
}
} // namespace
void MagnetOptions::validate() const {
  require(std::isfinite(start_displacement) && start_displacement >= 0 &&
              start_displacement <= 1,
          "invalid magnet start displacement");
  require(settle_window >= 1 && settle_window <= 1000 &&
              std::isfinite(settle_displacement) && settle_displacement >= 0 &&
              settle_displacement <= 1,
          "invalid magnet settling controls");
  require(std::isfinite(max_gap) && max_gap >= 0 && max_gap <= 100 &&
              std::isfinite(density_fraction) && density_fraction >= 0 &&
              density_fraction <= 1 && max_per_round >= 0 && max_rounds >= 0,
          "invalid magnet candidate controls");
}
const char *magnet_event_name(MagnetEvent event) {
  switch (event) {
  case MagnetEvent::timing:
    return "timing";
  case MagnetEvent::pull:
    return "pull";
  case MagnetEvent::terminal:
    return "terminal";
  case MagnetEvent::snap:
    return "snap";
  case MagnetEvent::round:
    return "magnet";
  case MagnetEvent::ready:
    return "ready";
  }
  throw std::invalid_argument("invalid magnet instability event");
}
MagnetStability::MagnetStability(MagnetOptions options) : options_(options) {
  options_.validate();
}
void MagnetStability::event(int iteration, MagnetEvent kind) {
  (void)magnet_event_name(kind);
  require(iteration >= last_event_, "magnet event iteration moved backwards");
  last_event_ = iteration;
  frames_.clear();
  movement_ = {};
}
MacroMovement MagnetStability::observe(const Design &d,
                                       const std::vector<double> &x,
                                       const std::vector<double> &y,
                                       int iteration) {
  require(iteration >= last_event_ && x.size() >= size_t(d.num_movable) &&
              y.size() >= size_t(d.num_movable),
          "invalid macro movement observation");
  Frame frame{iteration, {}};
  std::vector<double> extents;
  std::vector<int> nodes;
  for (int i = 0; i < d.num_movable; ++i) {
    if (d.node_macro.empty() || !d.node_macro[i])
      continue;
    const double hx = d.node_pad_x.empty() ? 0 : d.node_pad_x[i];
    const double hy = d.node_pad_y.empty() ? 0 : d.node_pad_y[i];
    const double extent =
        std::min(d.node_size_x[i] - 2 * hx, d.node_size_y[i] - 2 * hy);
    require(extent > 0 && std::isfinite(extent),
            "invalid macro movement extent");
    extents.push_back(extent);
    nodes.push_back(i);
    const auto cx = x[i] + d.node_size_x[i] / 2,
               cy = y[i] + d.node_size_y[i] / 2;
    require(std::isfinite(cx) && std::isfinite(cy), "nonfinite macro movement");
    frame.centers.emplace_back(cx, cy);
  }
  require(!frame.centers.empty(), "magnet stability requires movable macros");
  if (!frames_.empty() &&
      (frames_.back().iteration + 1 != iteration || nodes_ != nodes))
    frames_.clear();
  nodes_ = std::move(nodes);
  displacements_.assign(nodes_.size(), 0);
  frames_.push_back(std::move(frame));
  while (frames_.front().iteration < iteration - options_.settle_window)
    frames_.pop_front();
  movement_ = {};
  movement_.window_complete =
      frames_.front().iteration == iteration - options_.settle_window;
  double sum = 0;
  for (size_t i = 0; i < frames_.back().centers.size(); ++i) {
    double displacement = 0;
    const auto [x0, y0] = frames_.back().centers[i];
    for (const auto &old : frames_)
      displacement =
          std::max(displacement, std::hypot(x0 - old.centers[i].first,
                                            y0 - old.centers[i].second) /
                                     extents[i]);
    displacements_[i] = displacement;
    movement_.maximum = std::max(movement_.maximum, displacement);
    sum += displacement * displacement;
  }
  movement_.rms = std::sqrt(sum / frames_.back().centers.size());
  return movement_;
}
double MagnetStability::macro_movement(int node) const {
  const auto found = std::lower_bound(nodes_.begin(), nodes_.end(), node);
  require(found != nodes_.end() && *found == node,
          "missing macro movement observation");
  return displacements_[found - nodes_.begin()];
}
bool MagnetStability::macro_settled(int node) const {
  return movement_.window_complete &&
         macro_movement(node) <= options_.settle_displacement;
}
bool MagnetStability::settled(int iteration) const {
  return iteration >= last_event_ && movement_.window_complete &&
         movement_.maximum <= options_.settle_displacement;
}
bool MagnetStability::gate(int iteration) {
  const bool ready = settled(iteration);
  armed_ |= ready;
  return ready;
}
MagnetDensity::MagnetDensity(const Design &d, int bx, int by)
    : layout_(d.layout), bx_(bx), by_(by), dx_(0), dy_(0) {
  d.validate();
  require(bx >= 1 && bx <= 2048 && by >= 1 && by <= 2048,
          "invalid magnet density grid");
  dx_ = (layout_.xh - layout_.xl) / bx;
  dy_ = (layout_.yh - layout_.yl) / by;
  const double available =
      (layout_.xh - layout_.xl) * (layout_.yh - layout_.yl) -
      macro_union_area(d);
  double cell_area = 0;
  for (int i = 0; i < d.num_movable; ++i) {
    if (!cell(d, i))
      continue;
    const auto b = body(d, i);
    cell_area += (b.xh - b.xl) * (b.yh - b.yl);
  }
  require(available > 0 || cell_area == 0,
          "macro bodies leave no core area for standard cells");
  global_cell_density_ = cell_area == 0 ? 0 : cell_area / available;
  rho_.assign(bx * by, 0);
  for (int i = 0; i < d.num_movable; ++i) {
    if (!cell(d, i))
      continue;
    auto b = body(d, i);
    b.xl = std::clamp(b.xl, layout_.xl, layout_.xh);
    b.xh = std::clamp(b.xh, layout_.xl, layout_.xh);
    b.yl = std::clamp(b.yl, layout_.yl, layout_.yh);
    b.yh = std::clamp(b.yh, layout_.yl, layout_.yh);
    const auto [x0, x1] = detail::clamped_span(b.xl, b.xh, layout_.xl, dx_, bx);
    const auto [y0, y1] = detail::clamped_span(b.yl, b.yh, layout_.yl, dy_, by);
    for (int x = x0; x < x1; ++x)
      for (int y = y0; y < y1; ++y)
        rho_[x * by + y] +=
            overlap_area(b, {layout_.xl + x * dx_, layout_.yl + y * dy_,
                             layout_.xl + (x + 1) * dx_,
                             layout_.yl + (y + 1) * dy_}) /
            (dx_ * dy_);
  }
  prefix_.assign((bx + 1) * (by + 1), 0);
  for (int x = 1; x <= bx; ++x)
    for (int y = 1; y <= by; ++y)
      prefix_[x * (by + 1) + y] =
          rho_[(x - 1) * by + y - 1] + prefix_[(x - 1) * (by + 1) + y] +
          prefix_[x * (by + 1) + y - 1] - prefix_[(x - 1) * (by + 1) + y - 1];
}
double MagnetDensity::integral(double x, double y) const {
  x = std::clamp((x - layout_.xl) / dx_, 0., double(bx_));
  y = std::clamp((y - layout_.yl) / dy_, 0., double(by_));
  const int ix = std::min(int(x), bx_ - 1), iy = std::min(int(y), by_ - 1);
  const double tx = x - ix, ty = y - iy;
  const int stride = by_ + 1;
  return (1 - tx) * ((1 - ty) * prefix_[ix * stride + iy] +
                     ty * prefix_[ix * stride + iy + 1]) +
         tx * ((1 - ty) * prefix_[(ix + 1) * stride + iy] +
               ty * prefix_[(ix + 1) * stride + iy + 1]);
}
double MagnetDensity::density(const MagnetBox &box) const {
  const double xl = std::clamp(box.xl, layout_.xl, layout_.xh),
               xh = std::clamp(box.xh, layout_.xl, layout_.xh);
  const double yl = std::clamp(box.yl, layout_.yl, layout_.yh),
               yh = std::clamp(box.yh, layout_.yl, layout_.yh);
  if (xh <= xl || yh <= yl)
    return 0;
  return std::max(0., integral(xh, yh) - integral(xl, yh) - integral(xh, yl) +
                          integral(xl, yl)) *
         dx_ * dy_ / ((xh - xl) * (yh - yl));
}
double MagnetDensity::region_density(const MagnetRegion &region) const {
  const auto box = bounds(region);
  const double area =
      overlap_area(region, {layout_.xl, layout_.yl, layout_.xh, layout_.yh});
  if (area <= 0)
    return 0;
  // Preserve the constant-time summed-area query for aligned pairs.
  if (std::abs(area - overlap_area(box, {layout_.xl, layout_.yl, layout_.xh,
                                         layout_.yh})) <= area * 1e-12)
    return density(box);
  const auto [x0, x1] =
      detail::clamped_span(box.xl, box.xh, layout_.xl, dx_, bx_);
  const auto [y0, y1] =
      detail::clamped_span(box.yl, box.yh, layout_.yl, dy_, by_);
  double occupied = 0;
  for (int x = x0; x < x1; ++x)
    for (int y = y0; y < y1; ++y) {
      const double rho = rho_[x * by_ + y];
      if (rho == 0)
        continue;
      occupied += rho * overlap_area(region, {layout_.xl + x * dx_,
                                              layout_.yl + y * dy_,
                                              layout_.xl + (x + 1) * dx_,
                                              layout_.yl + (y + 1) * dy_});
    }
  return occupied / area;
}
MagnetSearch magnet_candidates(
    const Design &d, const MagnetOptions &options, const MagnetDensity &density,
    const MagnetStability *stability,
    const std::function<std::string(const MagnetCandidate &)> &check) {
  const auto started = std::chrono::steady_clock::now();
  d.validate();
  options.validate();
  std::map<std::pair<double, double>, int> shapes;
  std::vector<int> classes(d.node_x.size(), -1);
  for (int i = 0; i < d.num_movable; ++i) {
    if (d.node_macro.empty() || !d.node_macro[i])
      continue;
    // Shape is geometry metadata, independent of world-coordinate subtraction.
    // Fractional GP origins can otherwise give equal macros different keys.
    const auto key = std::pair{
        d.node_size_x[i] + (d.node_pad_x.empty() ? 2 * pad(d, i, true) : 0),
        d.node_size_y[i] + (d.node_pad_y.empty() ? 2 * pad(d, i, false) : 0)};
    if (!shapes.contains(key))
      shapes.emplace(key, shapes.size());
    classes[i] = shapes.at(key);
  }
  const auto pairs =
      detail::macro_grid_align_pairs(d, false, &classes, true, true);
  MagnetSearch result;
  result.global_cell_density = density.global_cell_density();
  result.density_limit = options.density_fraction * result.global_cell_density;
  result.exact_density = pairs.size() <= 64 && d.num_movable <= 10000;
  std::vector<std::pair<int, MagnetBox>> obstacles;
  for (int i = 0; i < int(d.node_x.size()); ++i)
    if (!d.node_macro.empty() && d.node_macro[i])
      obstacles.emplace_back(i, outer(d, i));
  for (const auto &pair : pairs) {
    MagnetCandidate c;
    c.a = pair.a;
    c.b = pair.b;
    c.axis = pair.axis == Axis::x ? Axis::y : Axis::x;
    auto a = outer(d, c.a), b = outer(d, c.b);
    const bool x = c.axis == Axis::x;
    const auto center = [&](const MagnetBox &box, bool along_x) {
      return along_x ? (box.xl + box.xh) / 2 : (box.yl + box.yh) / 2;
    };
    if (std::pair{center(b, x), c.b} < std::pair{center(a, x), c.a}) {
      std::swap(c.a, c.b);
      std::swap(a, b);
    }
    const double width =
        d.node_size_x[c.a] + (d.node_pad_x.empty() ? 2 * pad(d, c.a, true) : 0);
    const double height = d.node_size_y[c.a] +
                          (d.node_pad_y.empty() ? 2 * pad(d, c.a, false) : 0);
    const double extent = x ? width : height;
    const double orthogonal = x ? height : width;
    c.gap = center(b, x) - center(a, x) - extent;
    c.gap_region = x ? MagnetRegion{{{a.xh, a.yl},
                                     {std::max(a.xh, b.xl), b.yl},
                                     {std::max(a.xh, b.xl), b.yh},
                                     {a.xh, a.yh}}}
                     : MagnetRegion{{{a.xl, a.yh},
                                     {a.xh, a.yh},
                                     {b.xh, std::max(a.yh, b.yl)},
                                     {b.xl, std::max(a.yh, b.yl)}}};
    c.gap_box = bounds(c.gap_region);
    const auto area = [&](int i) {
      return (d.node_size_x[i] -
              (d.node_pad_x.empty() ? 0 : 2 * pad(d, i, true))) *
             (d.node_size_y[i] -
              (d.node_pad_y.empty() ? 0 : 2 * pad(d, i, false)));
    };
    const double area_a = area(c.a), area_b = area(c.b);
    const double cx = (center(a, true) * area_a + center(b, true) * area_b) /
                      (area_a + area_b);
    const double cy = (center(a, false) * area_a + center(b, false) * area_b) /
                      (area_a + area_b);
    const double w = width * (x ? 2 : 1), h = height * (x ? 1 : 2);
    const double preferred_x = cx - w / 2, preferred_y = cy - h / 2;
    // Preserve the area-weighted center where feasible. Near a core boundary,
    // a small halo overlap can put that ideal pose just outside the core even
    // though the complete union fits. Translate the whole rigid proposal by
    // the minimum amount instead of rejecting such a legal merge.
    const double merged_x =
        w <= d.layout.xh - d.layout.xl
            ? std::clamp(preferred_x, d.layout.xl, d.layout.xh - w)
            : preferred_x;
    const double merged_y =
        h <= d.layout.yh - d.layout.yl
            ? std::clamp(preferred_y, d.layout.yl, d.layout.yh - h)
            : preferred_y;
    c.center_shift_x = merged_x - preferred_x;
    c.center_shift_y = merged_y - preferred_y;
    c.merged_box = {merged_x, merged_y, merged_x + w, merged_y + h};
    c.merge.members = {c.a, c.b};
    c.merge.poses = {{c.merged_box.xl + pad(d, c.a, true),
                      c.merged_box.yl + pad(d, c.a, false)},
                     {c.merged_box.xl + (x ? extent : 0) + pad(d, c.b, true),
                      c.merged_box.yl + (x ? 0 : extent) + pad(d, c.b, false)}};
    for (size_t k = 0; k < 2; ++k) {
      const int i = c.merge.members[k];
      const auto old = body(d, i);
      c.displacement += std::hypot(c.merge.poses[k].first - old.xl,
                                   c.merge.poses[k].second - old.yl);
    }
    if (stability) {
      c.move_a = stability->macro_movement(c.a);
      c.move_b = stability->macro_movement(c.b);
    }
    c.gap_density = result.exact_density ? exact_density(d, c.gap_region)
                                         : density.region_density(c.gap_region);
    if (std::abs(center(a, !x) - center(b, !x)) >= orthogonal / 2)
      c.reason = "orthogonal_offset";
    else if (c.gap > options.max_gap * extent)
      c.reason = "gap_limit";
    else if (c.gap_density > result.density_limit)
      c.reason = "gap_density";
    else if (pad(d, c.a, true) != pad(d, c.b, true) ||
             pad(d, c.a, false) != pad(d, c.b, false))
      c.reason = "outside_halo_mismatch";
    else if (c.merged_box.xl < d.layout.xl || c.merged_box.yl < d.layout.yl ||
             c.merged_box.xh > d.layout.xh || c.merged_box.yh > d.layout.yh)
      c.reason = "outside_core";
    if (c.reason.empty())
      for (int i = d.num_movable; i < int(d.node_x.size()); ++i)
        if (overlap_area(c.merged_box, outer(d, i)) > 0) {
          c.reason = "fixed_overlap";
          break;
        }
    if (c.reason.empty())
      for (const auto &[i, box] : obstacles)
        if (i != c.a && i != c.b && overlap_area(c.gap_box, box) > 0 &&
            overlap_area(c.gap_region, box) > 0) {
          c.reason = "macro_in_gap";
          break;
        }
    if (c.reason.empty() && stability &&
        (!stability->macro_settled(c.a) || !stability->macro_settled(c.b)))
      c.reason = "macro_motion";
    result.candidates.push_back(std::move(c));
  }
  std::sort(result.candidates.begin(), result.candidates.end(),
            [](const auto &a, const auto &b) {
              return std::tie(a.gap_density, a.displacement, a.a, a.b, a.axis) <
                     std::tie(b.gap_density, b.displacement, b.a, b.b, b.axis);
            });
  std::set<int> members;
  std::vector<const MagnetCandidate *> chosen;
  for (auto &c : result.candidates) {
    if (!c.reason.empty())
      continue;
    if (members.contains(c.a) || members.contains(c.b))
      c.reason = "shared_member";
    else if (options.max_per_round &&
             result.merges.size() >= size_t(options.max_per_round))
      c.reason = "round_limit";
    if (c.reason.empty())
      for (const auto *other : chosen)
        if (overlap_area(c.gap_region, other->merged_box) > 0 ||
            overlap_area(other->gap_region, c.merged_box) > 0 ||
            overlap_area(c.merged_box, other->merged_box) > 0) {
          c.reason = "gap_conflict";
          break;
        }
    if (!c.reason.empty())
      continue;
    if (check)
      c.reason = check(c);
    if (!c.reason.empty())
      continue;
    c.selected = true;
    c.reason = "selected";
    members.insert(c.a);
    members.insert(c.b);
    chosen.push_back(&c);
    result.merges.push_back(c.merge);
  }
  result.seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - started)
          .count();
  return result;
}
} // namespace DreamPlace::native
