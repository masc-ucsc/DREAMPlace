// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

#include "dreamplace/native/rudy.hpp"
#include "dreamplace/native/spectral.hpp"

namespace DreamPlace::native {
namespace {
void require(bool condition, const char *message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}
bool finite(double v) { return std::isfinite(v); }
bool aligned(double v) { return std::abs(v - std::round(v)) <= 1e-6; }
int count(size_t n) {
  require(n <= static_cast<size_t>(std::numeric_limits<int>::max()),
          "array exceeds 32-bit kernel indices");
  return static_cast<int>(n);
}
void csr_check(const std::vector<int> &flat, const std::vector<int> &start,
               int objects, int pins, const std::vector<int> &inverse) {
  require(start.size() == static_cast<size_t>(objects) + 1 &&
              start.front() == 0 && start.back() == pins,
          "invalid CSR start/end");
  require(flat.size() == static_cast<size_t>(pins),
          "CSR must contain every pin exactly once");
  std::vector<bool> seen(pins, false);
  for (int n = 0; n < objects; ++n) {
    require(start[n] <= start[n + 1] && start[n] >= 0 && start[n + 1] <= pins,
            "nonmonotonic CSR");
    for (int j = start[n]; j < start[n + 1]; ++j) {
      int p = flat[j];
      require(p >= 0 && p < pins && !seen[p], "invalid or duplicate CSR pin");
      seen[p] = true;
      require(inverse[p] == n, "CSR disagrees with inverse map");
    }
  }
}

// Exact physical overflow, excluding fillers; fixed area reduces bin capacity.
double overflow(const Design &d, int bx, int by, double target_density,
                bool scale_macro = false) {
  const auto &l = d.layout;
  double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  std::vector<double> movable(bx * by), fixed(bx * by);
  double total = 0;
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
           h = d.node_size_y[i];
    const double charge = scale_macro && i < d.num_movable &&
                                  !d.node_macro.empty() && d.node_macro[i]
                              ? target_density
                              : 1;
    if (i < d.num_movable) {
      total += charge * w * h;
    }
    int x0 = std::clamp(static_cast<int>(std::floor((x - l.xl) / dx)), 0, bx);
    int x1 =
        std::clamp(static_cast<int>(std::ceil((x + w - l.xl) / dx)), 0, bx);
    int y0 = std::clamp(static_cast<int>(std::floor((y - l.yl) / dy)), 0, by);
    int y1 =
        std::clamp(static_cast<int>(std::ceil((y + h - l.yl) / dy)), 0, by);
    auto &map = i < d.num_movable ? movable : fixed;
    for (int k = x0; k < x1; ++k) {
      for (int j = y0; j < y1; ++j) {
        double a = std::max(0., std::min(x + w, l.xl + (k + 1) * dx) -
                                    std::max(x, l.xl + k * dx));
        double b = std::max(0., std::min(y + h, l.yl + (j + 1) * dy) -
                                    std::max(y, l.yl + j * dy));
        map[k * by + j] += charge * a * b;
      }
    }
  }
  double excess = 0;
  for (int i = 0; i < bx * by; ++i) {
    excess += std::max(0., movable[i] - target_density *
                                            std::max(0., dx * dy - fixed[i]));
  }
  return total ? excess / total : 0;
}

void options_check(const Options &o) {
  require(o.max_iterations >= 0 && o.bins_x >= 2 && o.bins_x <= 1024 &&
              o.bins_y >= 2 && o.bins_y <= 1024,
          "invalid iteration limit or bin count (2..1024)");
  require(o.num_threads > 0 && o.num_threads <= 128,
          "thread count must be 1..128");
  require(finite(o.target_density) && o.target_density > 0 &&
              o.target_density <= 1,
          "density must be in (0,1]");
  require(finite(o.stop_overflow) && o.stop_overflow >= 0 &&
              o.stop_overflow <= 1,
          "invalid overflow threshold");
  require(finite(o.initial_density_weight) && o.initial_density_weight > 0,
          "invalid density weight");
  require(finite(o.max_step) && o.max_step > 0 && o.max_step <= 1,
          "max_step must be in (0,1]");
  require(!o.detailed_place || o.legalize,
          "detailed placement requires legalization");
  require(finite(o.filler_width) && o.filler_width >= 0 &&
              finite(o.order_weight) && o.order_weight >= 0 &&
              finite(o.order_growth) && o.order_growth > 0 &&
              o.schedule_iterations >= 0,
          "invalid H2D GP hook options");
  require(!o.routability ||
              (finite(o.route_capacity_h) && o.route_capacity_h > 0 &&
               finite(o.route_capacity_v) && o.route_capacity_v > 0),
          "invalid routing capacity");
  require(finite(o.macro_overlap_weight) && o.macro_overlap_weight >= 0 &&
              finite(o.macro_jitter) && o.macro_jitter >= 0 &&
              o.macro_jitter <= .1,
          "invalid macro repulsion controls");
  require(
      finite(o.macro_align_weight) && o.macro_align_weight >= 0 &&
          finite(o.macro_align_sigma) && o.macro_align_sigma >= 0 &&
          finite(o.macro_align_overlap_boost) &&
          o.macro_align_overlap_boost >= 0 &&
          (!o.macro_align_weight || (o.macro_align_start >= 0 &&
                                     o.macro_align_start < o.macro_align_full)),
      "invalid macro alignment controls");
  require(finite(o.macro_snap_tolerance) && o.macro_snap_tolerance >= 0 &&
              o.macro_snap_weight > 0 && o.macro_snap_weight <= 1000000,
          "invalid macro snap controls");
  require(finite(o.cohesion_weight) && o.cohesion_weight >= 0,
          "invalid cohesion weight");
  require(!o.observer || o.observer_interval > 0,
          "observer interval must be positive");
}
bool has_halo(const Design &d) {
  const auto positive = [](double value) { return value > 0; };
  return std::any_of(d.node_halo_x.begin(), d.node_halo_x.end(), positive) ||
         std::any_of(d.node_halo_y.begin(), d.node_halo_y.end(), positive);
}
double halo(const std::vector<double> &values, int i, double pitch) {
  return values.empty() ? 0 : std::ceil(values[i] / pitch - 1e-9) * pitch;
}
Design padded_design(const Design &source) {
  Design d = source;
  for (auto &group : d.cohesion_groups) {
    if (group.areas.empty()) {
      for (int i : group.nodes) {
        group.areas.push_back(source.node_size_x[i] * source.node_size_y[i]);
      }
    }
  }
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    const double x = halo(source.node_halo_x, i, d.layout.site_width);
    const double y = halo(source.node_halo_y, i, d.layout.row_height);
    if (!d.node_boxes.empty()) {
      d.node_boxes[i].xl -= x;
      d.node_boxes[i].xh -= x;
      d.node_boxes[i].yl -= y;
      d.node_boxes[i].yh -= y;
    }
    d.node_x[i] -= x;
    d.node_y[i] -= y;
    d.node_size_x[i] += 2 * x;
    d.node_size_y[i] += 2 * y;
    for (int k = d.flat_node2pin_start[i]; k < d.flat_node2pin_start[i + 1];
         ++k) {
      const int p = d.flat_node2pin[k];
      d.pin_offset_x[p] += x;
      d.pin_offset_y[p] += y;
    }
  }
  for (auto &edge : d.order_edges) {
    const double a =
        edge.axis == Axis::x
            ? halo(source.node_halo_x, edge.from, d.layout.site_width)
            : halo(source.node_halo_y, edge.from, d.layout.row_height);
    const double b =
        edge.axis == Axis::x
            ? halo(source.node_halo_x, edge.to, d.layout.site_width)
            : halo(source.node_halo_y, edge.to, d.layout.row_height);
    // Preserve a relation on physical footprints under origin padding.
    edge.separation -= a + b;
  }
  d.node_halo_x.clear();
  d.node_halo_y.clear();
  return d;
}
} // namespace

void Design::build_maps() {
  int nodes = count(node_x.size()), pins = count(pin2node.size());
  require(!flat_net2pin_start.empty(), "missing net CSR");
  int nets = count(flat_net2pin_start.size() - 1);
  // Build into temporary arrays so a malformed design does not leave partial
  // maps.
  std::vector<int> p2n(pins, -1);
  require(flat_net2pin.size() == static_cast<size_t>(pins) &&
              flat_net2pin_start.front() == 0 &&
              flat_net2pin_start.back() == pins,
          "invalid net CSR");
  for (int n = 0; n < nets; ++n) {
    int a = flat_net2pin_start[n], b = flat_net2pin_start[n + 1];
    require(a >= 0 && a <= b && b <= pins, "nonmonotonic net CSR");
    for (int j = a; j < b; ++j) {
      int p = flat_net2pin[j];
      require(p >= 0 && p < pins && p2n[p] == -1,
              "invalid or duplicate net pin");
      p2n[p] = n;
    }
  }
  std::vector<int> start(nodes + 1), flat(pins);
  for (int i : pin2node) {
    require(i >= 0 && i < nodes, "invalid pin2node");
    ++start[i + 1];
  }
  std::partial_sum(start.begin(), start.end(), start.begin());
  auto next = start;
  for (int p = 0; p < pins; ++p) {
    flat[next[pin2node[p]]++] = p;
  }
  if (net_weights.empty()) {
    net_weights.assign(nets, 1);
  }
  if (net_mask.empty()) {
    net_mask.resize(nets);
    for (int n = 0; n < nets; ++n) {
      net_mask[n] = flat_net2pin_start[n + 1] - flat_net2pin_start[n] > 1;
    }
  }
  pin2net = std::move(p2n);
  flat_node2pin_start = std::move(start);
  flat_node2pin = std::move(flat);
}

void Design::validate() const {
  int n = count(node_x.size()), p = count(pin2node.size()),
      nets = count(net_weights.size());
  require(num_movable >= 0 && num_movable <= n, "invalid movable count");
  const auto &l = layout;
  require(finite(l.xl) && finite(l.yl) && finite(l.xh) && finite(l.yh) &&
              l.xh > l.xl && l.yh > l.yl && finite(l.xh - l.xl) &&
              finite(l.yh - l.yl),
          "invalid core bounds");
  require(finite(l.site_width) && l.site_width > 0 && finite(l.row_height) &&
              l.row_height > 0,
          "invalid site/row pitch");
  require((l.xh - l.xl) / l.site_width < std::numeric_limits<int>::max() &&
              (l.yh - l.yl) / l.row_height < std::numeric_limits<int>::max(),
          "layout exceeds kernel site indices");
  require(aligned((l.xh - l.xl) / l.site_width) &&
              aligned((l.yh - l.yl) / l.row_height),
          "core must contain whole sites/rows");
  require(node_y.size() == node_x.size() &&
              node_size_x.size() == node_x.size() &&
              node_size_y.size() == node_x.size(),
          "node array size mismatch");
  require((node_halo_x.empty() || node_halo_x.size() == node_x.size()) &&
              (node_halo_y.empty() || node_halo_y.size() == node_x.size()),
          "halo array size mismatch");
  require(node_orientation.empty() || node_orientation.size() == node_x.size(),
          "orientation array size mismatch");
  require(node_boxes.empty() || node_boxes.size() == node_x.size(),
          "position box size mismatch");
  for (const auto *capacity :
       {&node_blocked_capacity_h, &node_blocked_capacity_v}) {
    require(capacity->empty() ||
                (capacity->size() == node_x.size() &&
                 std::all_of(capacity->begin(), capacity->end(),
                             [](double v) { return finite(v) && v >= 0; })),
            "invalid macro blocked capacity");
  }
  require(node_align_class.empty() ||
              (node_align_class.size() == node_x.size() &&
               std::all_of(node_align_class.begin(), node_align_class.end(),
                           [](int v) { return v >= -1; })),
          "invalid macro alignment classes");
  require(node_macro.empty() || node_macro.size() == node_x.size(),
          "macro mask size mismatch");
  require(std::all_of(node_macro.begin(), node_macro.end(),
                      [](auto v) { return v <= 1; }),
          "invalid macro mask");
  require(node_inflatable.empty() || node_inflatable.size() == node_x.size(),
          "inflation mask size mismatch");
  for (const auto &e : order_edges) {
    require(e.from >= 0 && e.from < n && e.to >= 0 && e.to < n &&
                e.from != e.to && (e.axis == Axis::x || e.axis == Axis::y) &&
                finite(e.separation) && finite(e.strength.weight) &&
                e.strength.weight >= 0,
            "invalid relative order");
  }
  for (int i = 0; i < n; ++i) {
    double w = node_size_x[i], h = node_size_y[i];
    if (!node_orientation.empty()) {
      require(valid_orientation(node_orientation[i]),
              "invalid node orientation");
      require((w > 0 && h > 0) || node_orientation[i] == Orientation::N,
              "terminal cannot be flipped");
    }
    const double hx = node_halo_x.empty() ? 0 : node_halo_x[i];
    const double hy = node_halo_y.empty() ? 0 : node_halo_y[i];
    require(finite(hx) && finite(hy) && hx >= 0 && hy >= 0,
            "invalid node halo");
    require((hx == 0 && hy == 0) || (w > 0 && h > 0),
            "halo requires a positive-area node");
    require(finite(node_x[i]) && finite(node_y[i]) && finite(w) && finite(h) &&
                w >= 0 && h >= 0 && finite(w * h),
            "nonfinite position or invalid node size");
    require(finite(node_x[i] + w) && finite(node_y[i] + h),
            "node bounds overflow");
    if (!node_boxes.empty()) {
      const auto &b = node_boxes[i];
      require(finite(b.xl) && finite(b.yl) && finite(b.xh) && finite(b.yh) &&
                  b.xl <= b.xh && b.yl <= b.yh,
              "invalid position box");
      if (i < num_movable) {
        require(std::max(l.xl, b.xl) <= std::min(l.xh - w, b.xh) &&
                    std::max(l.yl, b.yl) <= std::min(l.yh - h, b.yh),
                "position box outside core");
      }
    }
    if (!node_inflatable.empty()) {
      require(node_inflatable[i] <= 1, "invalid inflation mask");
    }
    if (i < num_movable) {
      require(w > 0 && h > 0 && w <= l.xh - l.xl && h <= l.yh - l.yl,
              "movable node does not fit core");
      // Macro legalization uses private site/row-rounded occupancy envelopes.
      // Keep the caller's physical dimensions and pin offsets unchanged.
    }
  }
  std::vector<bool> grouped(node_x.size(), false);
  require(centroid_targets.empty() ||
              centroid_targets.size() == centroid_groups.size(),
          "centroid target size mismatch");
  for (auto [x, y] : centroid_targets)
    require(finite(x) && finite(y), "nonfinite frozen centroid target");
  if (!centroid_groups.empty()) {
    require(centroid_node_area.size() == node_x.size() &&
                node_centroid_color.size() == node_x.size(),
            "centroid physical area/color metadata size mismatch");
    std::vector<bool> populated(centroid_color_names.size(), false);
    for (size_t n = 0; n < node_x.size(); ++n) {
      require(finite(centroid_node_area[n]) && centroid_node_area[n] >= 0,
              "invalid centroid physical area");
      const int c = node_centroid_color[n];
      require(c >= -1 && c < static_cast<int>(centroid_color_names.size()) &&
                  (c < 0 || centroid_node_area[n] > 0),
              "invalid centroid color identity");
      if (c >= 0)
        populated[c] = true;
    }
    for (size_t g = 0; g < centroid_groups.size(); ++g) {
      const auto &group = centroid_groups[g];
      require(group.parent >= -1 && group.parent < static_cast<int>(g) &&
                  group.family >= 0,
              "invalid centroid family/parent");
      require(!group.nodes.empty() || !group.colors.empty(),
              "empty centroid group");
      std::set<int> nodes, colors;
      for (int n : group.nodes)
        require(n >= 0 && n < static_cast<int>(node_x.size()) &&
                    nodes.insert(n).second && node_centroid_color[n] < 0 &&
                    centroid_node_area[n] > 0,
                "invalid or duplicate uncolored centroid member");
      for (int c : group.colors)
        require(c >= 0 && c < static_cast<int>(centroid_color_names.size()) &&
                    colors.insert(c).second && populated[c],
                "invalid or duplicate centroid anchor color");
    }
  }
  for (const auto &group : cohesion_groups) {
    require(finite(group.weight) && group.weight >= 0,
            "invalid cohesion group weight");
    require(group.areas.empty() || group.areas.size() == group.nodes.size(),
            "cohesion area size mismatch");
    for (size_t k = 0; k < group.nodes.size(); ++k) {
      int i = group.nodes[k];
      require(i >= 0 && i < n && !grouped[i],
              "invalid or overlapping cohesion membership");
      grouped[i] = true;
      double a = group.areas.empty() ? node_size_x[i] * node_size_y[i]
                                     : group.areas[k];
      require(finite(a) && a > 0 && node_size_x[i] > 0 && node_size_y[i] > 0,
              "cohesion requires physical positive-area nodes");
    }
  }
  require(pin2net.size() == pin2node.size() &&
              pin_offset_x.size() == pin2node.size() &&
              pin_offset_y.size() == pin2node.size() &&
              net_mask.size() == net_weights.size(),
          "pin/net array size mismatch");
  for (int i = 0; i < p; ++i) {
    require(pin2node[i] >= 0 && pin2node[i] < n && pin2net[i] >= 0 &&
                pin2net[i] < nets,
            "invalid pin map index");
    require(finite(pin_offset_x[i]) && finite(pin_offset_y[i]),
            "nonfinite pin offset");
  }
  for (int i = 0; i < nets; ++i) {
    require(finite(net_weights[i]) && net_weights[i] >= 0 && net_mask[i] <= 1,
            "invalid net weight/mask");
  }
  csr_check(flat_net2pin, flat_net2pin_start, nets, p, pin2net);
  csr_check(flat_node2pin, flat_node2pin_start, n, p, pin2node);
  for (int i = 0; i < nets; ++i) {
    require(!net_mask[i] ||
                flat_net2pin_start[i + 1] - flat_net2pin_start[i] >= 2,
            "active net needs at least two pins");
  }
}

double hpwl(const Design &d) {
  d.validate();
  double wl = 0;
  for (int n = 0; n < count(d.net_weights.size()); ++n) {
    if (!d.net_mask[n]) {
      continue;
    }
    double xmin = std::numeric_limits<double>::infinity(), ymin = xmin,
           xmax = -xmin, ymax = -xmin;
    for (int j = d.flat_net2pin_start[n]; j < d.flat_net2pin_start[n + 1];
         ++j) {
      int p = d.flat_net2pin[j], i = d.pin2node[p];
      double x = d.node_x[i] + d.pin_offset_x[p],
             y = d.node_y[i] + d.pin_offset_y[p];
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    wl += d.net_weights[n] * (xmax - xmin + ymax - ymin);
  }
  return wl;
}

bool is_legal(const Design &d) {
  d.validate();
  if (has_halo(d)) {
    return is_legal(padded_design(d));
  }
  const auto &l = d.layout;
  constexpr double eps = 1e-6;
  struct Interval {
    double x, end;
    int id;
  };
  std::vector<std::vector<Interval>> rows(
      static_cast<size_t>(std::round((l.yh - l.yl) / l.row_height)));
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i],
           h = d.node_size_y[i];
    if (!w || !h) {
      continue;
    }
    if (i < d.num_movable &&
        (x < l.xl - eps || y < l.yl - eps || x + w > l.xh + eps ||
         y + h > l.yh + eps || !aligned((x - l.xl) / l.site_width) ||
         !aligned((y - l.yl) / l.row_height))) {
      return false;
    }
    int first = std::max(0., std::floor((y - l.yl + eps) / l.row_height));
    int last = std::min<double>(rows.size(),
                                std::ceil((y + h - l.yl - eps) / l.row_height));
    for (int r = first; r < last; ++r) {
      rows[r].push_back({x, x + w, i});
    }
  }
  for (auto &row : rows) {
    std::sort(row.begin(), row.end(), [](auto a, auto b) { return a.x < b.x; });
    double end_all = -std::numeric_limits<double>::infinity(),
           end_movable = end_all;
    for (auto v : row) {
      if (v.id < d.num_movable ? v.x < end_all - eps
                               : v.x < end_movable - eps) {
        return false;
      }
      end_all = std::max(end_all, v.end);
      if (v.id < d.num_movable) {
        end_movable = std::max(end_movable, v.end);
      }
    }
  }
  return true;
}

const char *stage_name(PlacementStage stage) {
  switch (stage) {
  case PlacementStage::centroid_start:
    return "centroid attraction start";
  case PlacementStage::initial:
    return "initial";
  case PlacementStage::global_placement:
    return "global placement";
  case PlacementStage::legalization:
    return "legalization";
  case PlacementStage::detailed_placement:
    return "detailed placement";
  case PlacementStage::finished:
    return "finished";
  }
  return "unknown";
}

Result place(Design &design, const Options &options, GpState *state) {
  design.validate();
  Options o = options;
  if (o.automatic_bins) {
    const auto bins = detail::automatic_bins(design);
    o.bins_x = bins.first;
    o.bins_y = bins.second;
    o.automatic_bins = false;
  }
  options_check(o);
  require(
      finite(o.centroid_weight) && o.centroid_weight >= 0 &&
          finite(o.centroid_repulsion_weight) &&
          o.centroid_repulsion_weight >= 0 && o.centroid_ramp > 0 &&
          (!o.centroid_repulsion_weight || o.centroid_weight > 0) &&
          finite(o.centroid_stop_overlap) && o.centroid_stop_overlap >= 0 &&
          o.centroid_stop_overlap < 1 &&
          (!o.centroid_adaptive_repulsion || o.centroid_repulsion_weight > 0),
      "invalid centroid attraction/repulsion controls");
  const bool centroids =
      o.centroid_weight > 0 && !design.centroid_groups.empty();
  const bool separate_circles = centroids && o.centroid_adaptive_repulsion;
  require(
      !centroids || !state,
      "centroid placement requires one GP solve without external continuation");
  if (has_halo(design)) {
    Design padded = padded_design(design);
    Options control = o;
    if (o.observer) {
      control.observer = [&](const Design &snapshot, PlacementStage stage,
                             int iteration) {
        Design physical = design;
        for (int i = 0; i < physical.num_movable; ++i) {
          physical.node_x[i] =
              snapshot.node_x[i] +
              halo(design.node_halo_x, i, design.layout.site_width);
          physical.node_y[i] =
              snapshot.node_y[i] +
              halo(design.node_halo_y, i, design.layout.row_height);
        }
        o.observer(physical, stage, iteration);
      };
    }
    const auto result = place(padded, control, state);
    for (int i = 0; i < design.num_movable; ++i) {
      design.node_x[i] = padded.node_x[i] +
                         halo(design.node_halo_x, i, design.layout.site_width);
      design.node_y[i] = padded.node_y[i] +
                         halo(design.node_halo_y, i, design.layout.row_height);
    }
    return result;
  }
  GpState next;
  const bool warm = state && state->signature;
  // FNV over geometry/connectivity. Coordinates, pin offsets and weights may
  // evolve between H2D checkpoints; changing the solver topology cannot reuse
  // it.
  std::uint64_t signature = 14695981039346656037ULL;
  auto salt = [&](std::uint64_t v) {
    signature = (signature ^ v) * 1099511628211ULL;
  };
  salt(design.num_movable);
  salt(o.bins_x);
  salt(o.bins_y);
  salt(o.add_fillers);
  salt(o.seed);
  for (double v :
       {design.layout.xl, design.layout.yl, design.layout.xh, design.layout.yh,
        design.layout.site_width, design.layout.row_height, o.target_density,
        o.filler_width}) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  for (double v : design.node_size_x) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  for (double v : design.node_size_y) {
    salt(std::bit_cast<std::uint64_t>(v));
  }
  for (int v : design.pin2node) {
    salt(v);
  }
  for (int v : design.flat_net2pin) {
    salt(v);
  }
  for (int v : design.flat_net2pin_start) {
    salt(v);
  }
  if (o.cohesion_weight > 0) {
    salt(std::bit_cast<std::uint64_t>(o.cohesion_weight));
    for (const auto &group : design.cohesion_groups) {
      salt(group.nodes.size());
      salt(std::bit_cast<std::uint64_t>(group.weight));
      for (int i : group.nodes) {
        salt(i);
      }
      for (double a : group.areas) {
        salt(std::bit_cast<std::uint64_t>(a));
      }
    }
  }
  if (o.scale_macro_density || o.exact_macro_density || o.adaptive_density ||
      o.macro_overlap_weight || o.macro_jitter || o.improved_routability) {
    salt(o.scale_macro_density);
    salt(o.exact_macro_density);
    salt(o.adaptive_density);
    salt(o.improved_routability);
    if (o.improved_routability) {
      salt(std::bit_cast<std::uint64_t>(o.route_capacity_h));
      salt(std::bit_cast<std::uint64_t>(o.route_capacity_v));
      for (const auto *values :
           {&design.node_blocked_capacity_h, &design.node_blocked_capacity_v})
        for (double v : *values)
          salt(std::bit_cast<std::uint64_t>(v));
      for (auto v : design.node_inflatable)
        salt(v);
    }
    salt(std::bit_cast<std::uint64_t>(o.macro_overlap_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_jitter));
    for (auto m : design.node_macro) {
      salt(m);
    }
  }
  if (o.improved_routability && !o.routability_charge_compensation)
    salt(0x524f555445434847ULL);
  if (o.macro_align_weight > 0) {
    salt(std::bit_cast<std::uint64_t>(o.macro_align_weight));
    salt(std::bit_cast<std::uint64_t>(o.macro_align_sigma));
    salt(std::bit_cast<std::uint64_t>(o.macro_align_overlap_boost));
    salt(o.macro_align_start);
    salt(o.macro_align_full);
    salt(o.macro_align_nearest);
    for (auto c : design.node_align_class)
      salt(c);
    for (auto m : design.node_macro)
      salt(m);
  }
  if (warm && state->signature != signature) {
    throw std::invalid_argument("GP continuation options changed");
  }
  Result r;
  r.bins_x = o.bins_x;
  r.bins_y = o.bins_y;
  r.initial_hpwl = hpwl(design);
  Design d = design;
  const int original = count(d.node_x.size());
  if (!d.node_boxes.empty()) {
    for (int i = 0; i < d.num_movable; ++i) {
      const auto &b = d.node_boxes[i];
      d.node_x[i] = std::clamp(d.node_x[i], std::max(d.layout.xl, b.xl),
                               std::min(d.layout.xh - d.node_size_x[i], b.xh));
      d.node_y[i] = std::clamp(d.node_y[i], std::max(d.layout.yl, b.yl),
                               std::min(d.layout.yh - d.node_size_y[i], b.yh));
    }
  }
  if (o.observer) {
    o.observer(d, PlacementStage::initial, warm ? state->iteration : 0);
  }
  const auto &l = d.layout;
  const double width = l.xh - l.xl, height = l.yh - l.yl;
  double dx = width / o.bins_x, dy = height / o.bins_y;
  // A local engine owns all temporary/filler state. No global RNG or settings.
  std::mt19937 rng(o.seed);
  std::uniform_real_distribution<double> unit(0, 1);
  double area = 0, fixed_area = 0, macro_area = 0;
  for (int i = 0; i < original; ++i) {
    double a = d.node_size_x[i] * d.node_size_y[i];
    if (i < d.num_movable) {
      area += a;
      if (!d.node_macro.empty() && d.node_macro[i]) {
        macro_area += a;
      }
    } else {
      double w = std::max(0., std::min(l.xh, d.node_x[i] + d.node_size_x[i]) -
                                  std::max(l.xl, d.node_x[i]));
      double h = std::max(0., std::min(l.yh, d.node_y[i] + d.node_size_y[i]) -
                                  std::max(l.yl, d.node_y[i]));
      fixed_area += w * h;
    }
  }
  const auto gp_start = std::chrono::steady_clock::now();
  if (d.num_movable && o.max_iterations) {
    auto active = [&](int i) { return i < d.num_movable || i >= original; };
    auto project = [&](int i, double &x, double &y) {
      const PositionBox box =
          i < original && !d.node_boxes.empty()
              ? d.node_boxes[i]
              : PositionBox{l.xl, l.yl, l.xh - d.node_size_x[i],
                            l.yh - d.node_size_y[i]};
      x = std::clamp(x, std::max(l.xl, box.xl),
                     std::min(l.xh - d.node_size_x[i], box.xh));
      y = std::clamp(y, std::max(l.yl, box.yl),
                     std::min(l.yh - d.node_size_y[i], box.yh));
    };
    for (int i = 0; i < d.num_movable; ++i) {
      if (!warm && o.macro_jitter > 0 && !d.node_macro.empty() &&
          d.node_macro[i]) {
        // Index-local RNG avoids changing fillers or other nodes when a mask
        // changes.
        std::mt19937 local(o.seed ^ (0x9e3779b9U * (unsigned(i) + 1)));
        std::uniform_real_distribution<double> jitter(-o.macro_jitter,
                                                      o.macro_jitter);
        d.node_x[i] += jitter(local) * d.node_size_x[i];
        d.node_y[i] += jitter(local) * d.node_size_y[i];
      }
      project(i, d.node_x[i], d.node_y[i]);
    }
    if (centroids) {
      r.centroid.started = true;
      if (o.observer)
        o.observer(d, PlacementStage::centroid_start, 0);
    }
    if (o.add_fillers) {
      double target_area = std::max(
          0., o.target_density * (width * height - fixed_area) - area +
                  (o.scale_macro_density ? (1 - o.target_density) * macro_area
                                         : 0));
      double fw = o.filler_width > 0 ? o.filler_width
                                     : area / d.num_movable / l.row_height;
      double fh = l.row_height;
      fw = std::clamp(fw, l.site_width, width);
      double desired = std::floor(target_area / (fw * fh));
      require(desired <= 1000000,
              "too many fillers; reduce core whitespace or disable fillers");
      int fillers = desired;
      for (int i = 0; i < fillers; ++i) {
        d.node_x.push_back(l.xl + unit(rng) * (width - fw));
        d.node_y.push_back(l.yl + unit(rng) * (height - fh));
        d.node_size_x.push_back(fw);
        d.node_size_y.push_back(fh);
      }
      d.build_maps();
    }
    const int n = count(d.node_x.size());
    if (warm) {
      for (const auto *a :
           {&state->x, &state->y, &state->prev_x, &state->prev_y,
            &state->old_gx, &state->old_gy, &state->old_ex, &state->old_ey,
            &state->inflation}) {
        require(a->size() == static_cast<size_t>(n) &&
                    std::all_of(a->begin(), a->end(),
                                [](double v) { return finite(v); }),
                "invalid GP continuation arrays");
      }
      require(state->iteration >= 0 && finite(state->density_weight) &&
                  state->density_weight > 0 && finite(state->momentum) &&
                  state->momentum >= 1 && finite(state->step) &&
                  state->step >= 0,
              "invalid GP continuation scalars");
    }
    std::vector<double> inflation(n, 1);
    if (warm) {
      inflation = state->inflation;
    }
    std::vector<double> sx(n), sy(n), ox(n), oy(n), ratio(n);
    const auto sharp = [&](int i) {
      return i < original && o.exact_macro_density && !d.node_macro.empty() &&
             d.node_macro[i];
    };
    for (int i = 0; i < n; ++i) {
      sx[i] = active(i) && !sharp(i)
                  ? std::max(d.node_size_x[i], std::sqrt(2.) * dx)
                  : d.node_size_x[i];
      sy[i] = active(i) && !sharp(i)
                  ? std::max(d.node_size_y[i], std::sqrt(2.) * dy)
                  : d.node_size_y[i];
      ox[i] = (d.node_size_x[i] - sx[i]) / 2;
      oy[i] = (d.node_size_y[i] - sy[i]) / 2;
      ratio[i] = sx[i] * sy[i]
                     ? d.node_size_x[i] * d.node_size_y[i] / (sx[i] * sy[i])
                     : 0;
      if (!active(i) || (o.scale_macro_density && i < original &&
                         !d.node_macro.empty() && d.node_macro[i])) {
        ratio[i] *= o.target_density;
      }
    }
    PoissonSolver solver(o.bins_x, o.bins_y, width, height);
    std::vector<double> map, phi, fx, fy;
    auto x = d.node_x, y = d.node_y, prev_x = x, prev_y = y;
    std::vector<double> old_gx(n), old_gy(n), old_ex(n), old_ey(n);
    double density_weight = o.initial_density_weight, t = 1, alpha = 0;
    const int offset = warm ? state->iteration : 0;
    double previous_overflow = warm ? state->previous_overflow : -1;
    double previous_hpwl = warm ? state->previous_hpwl : r.initial_hpwl;
    double route_overflow = 0, inflated_ratio = 1;
    std::vector<AlignmentPair> align_pairs =
        warm ? state->align_pairs : std::vector<AlignmentPair>{};
    double align_residual = 0, align_ramp = 0;
    double centroid_repulsion_weight = o.centroid_repulsion_weight;
    double previous_circle_overlap = 1;
    if (warm && o.improved_routability) {
      double inflated = 0;
      for (int i = 0; i < d.num_movable; ++i)
        inflated += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
      inflated_ratio = area ? inflated / area : 1;
    }
    if (warm) {
      for (const auto &p : state->align_pairs) {
        require(
            p.a >= 0 && p.b >= 0 && p.a < d.num_movable &&
                p.b < d.num_movable && p.a != p.b &&
                (p.axis == Axis::x || p.axis == Axis::y) &&
                !d.node_align_class.empty() && d.node_align_class[p.a] >= 0 &&
                d.node_align_class[p.a] == d.node_align_class[p.b] &&
                !d.node_macro.empty() && d.node_macro[p.a] && d.node_macro[p.b],
            "invalid GP continuation alignment pairs");
      }
      x = state->x;
      y = state->y;
      prev_x = state->prev_x;
      prev_y = state->prev_y;
      old_gx = state->old_gx;
      old_gy = state->old_gy;
      old_ex = state->old_ex;
      old_ey = state->old_ey;
      density_weight = state->density_weight;
      t = state->momentum;
      alpha = state->step;
      // External compaction may have moved physical nodes since the checkpoint.
      for (int i = 0; i < original; ++i) {
        if (i >= d.num_movable) {
          // Refreshed fixed pseudo-terminals must retain the caller's exact
          // coordinates. Translating old state by (new - old) can lose an ulp.
          x[i] = prev_x[i] = old_ex[i] = d.node_x[i];
          y[i] = prev_y[i] = old_ey[i] = d.node_y[i];
          continue;
        }
        const double dx = d.node_x[i] - x[i], dy = d.node_y[i] - y[i];
        x[i] += dx;
        y[i] += dy;
        prev_x[i] += dx;
        prev_y[i] += dy;
        old_ex[i] += dx;
        old_ey[i] += dy;
      }
    }
    for (int iter = 0; iter < o.max_iterations; ++iter) {
      const int global_iter = offset + iter;
      double progress = std::min(
          1., double(global_iter) /
                  std::max(1, (o.schedule_iterations ? o.schedule_iterations
                                                     : o.max_iterations) -
                                  1));
      double gamma = std::max(l.site_width,
                              4 * std::max(dx, dy) * std::pow(0.05, progress));
      double nt = (1 + std::sqrt(1 + 4 * t * t)) / 2, momentum = (t - 1) / nt;
      for (int i = 0; i < n; ++i) {
        if (active(i)) {
          d.node_x[i] = x[i] + momentum * (x[i] - prev_x[i]);
          d.node_y[i] = y[i] + momentum * (y[i] - prev_y[i]);
          project(i, d.node_x[i], d.node_y[i]);
        }
      }
      auto wg = detail::wirelength(d, gamma, o.num_threads);
      const double lambda =
          o.order_weight * std::pow(o.order_growth, global_iter);
      require(finite(lambda), "relative-order schedule overflow");
      if (lambda) {
        for (const auto &e : d.order_edges) {
          const bool horizontal = e.axis == Axis::x;
          const auto &pos = horizontal ? d.node_x : d.node_y;
          const auto &size = horizontal ? d.node_size_x : d.node_size_y;
          const double violation = std::max(0., pos[e.from] + size[e.from] +
                                                    e.separation - pos[e.to]);
          auto &g = horizontal ? wg.x : wg.y;
          const double force = 2 * lambda * e.strength.weight * violation;
          g[e.from] += force;
          g[e.to] -= force;
        }
      }
      if (o.routability && global_iter % 20 == 0 &&
          (!o.improved_routability ||
           (previous_overflow >= 0 && previous_overflow < .3))) {
        Design physical = design;
        std::copy_n(d.node_x.begin(), original, physical.node_x.begin());
        std::copy_n(d.node_y.begin(), original, physical.node_y.begin());
        const auto route =
            rudy(physical,
                 {o.bins_x, o.bins_y, o.route_capacity_h, o.route_capacity_v,
                  o.improved_routability, o.improved_routability});
        route_overflow = route.overflow;
        auto proposed = inflation;
        double added = 0, current = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          current += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
        }
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_inflatable.empty() && d.node_inflatable[i]) {
            const int bx = std::clamp(
                int((d.node_x[i] + d.node_size_x[i] / 2 - l.xl) / dx), 0,
                o.bins_x - 1);
            const int by = std::clamp(
                int((d.node_y[i] + d.node_size_y[i] / 2 - l.yl) / dy), 0,
                o.bins_y - 1);
            const double demand =
                std::clamp(std::max({route.pin_utilization[bx * o.bins_y + by],
                                     route.horizontal[bx * o.bins_y + by],
                                     route.vertical[bx * o.bins_y + by]}),
                           1., 3.);
            proposed[i] = o.improved_routability
                              ? std::max(inflation[i], demand)
                              : demand;
            added += (proposed[i] - inflation[i]) * d.node_size_x[i] *
                     d.node_size_y[i];
          }
        }
        const double budget = o.improved_routability && added > 0
                                  ? std::min(1., .05 * current / added)
                                  : 1.;
        double inflated = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          inflation[i] += budget * (proposed[i] - inflation[i]);
          inflated += d.node_size_x[i] * d.node_size_y[i] * inflation[i];
        }
        inflated_ratio = area ? inflated / area : 1;
      }
      const double target =
          o.target_density / (o.improved_routability ? inflated_ratio : 1);
      double filler_scale = 1;
      if (o.improved_routability && o.routability_charge_compensation &&
          n > original) {
        double filler_area = 0;
        for (int i = original; i < n; ++i)
          filler_area += d.node_size_x[i] * d.node_size_y[i];
        // Keep total charge consistent with the effective capacity. A uniform
        // target shift is the discarded DC mode and cannot change forces.
        const double real_charge =
            area * inflated_ratio -
            (o.scale_macro_density ? (1 - target) * macro_area : 0);
        const double budget =
            std::max(0., target * (width * height - fixed_area) - real_charge);
        filler_scale = filler_area ? budget / filler_area : 0;
        for (int i = original; i < n; ++i)
          ratio[i] = filler_scale * d.node_size_x[i] * d.node_size_y[i] /
                     (sx[i] * sy[i]);
        for (int i = d.num_movable; i < original; ++i)
          ratio[i] = sx[i] * sy[i] > 0 ? target * d.node_size_x[i] *
                                             d.node_size_y[i] / (sx[i] * sy[i])
                                       : 0;
      }
      if (o.routability) {
        for (int i = 0; i < d.num_movable; ++i) {
          const double factor = std::sqrt(inflation[i]);
          sx[i] = std::max(d.node_size_x[i] * factor,
                           sharp(i) ? 0. : std::sqrt(2.) * dx);
          sy[i] = std::max(d.node_size_y[i] * factor,
                           sharp(i) ? 0. : std::sqrt(2.) * dy);
          ox[i] = (d.node_size_x[i] - sx[i]) / 2;
          oy[i] = (d.node_size_y[i] - sy[i]) / 2;
          ratio[i] = d.node_size_x[i] * d.node_size_y[i] * inflation[i] /
                     (sx[i] * sy[i]);
          if (o.scale_macro_density && !d.node_macro.empty() &&
              d.node_macro[i]) {
            ratio[i] *=
                o.routability_charge_compensation ? target : o.target_density;
          }
        }
      }
      detail::density_map(d, o.bins_x, o.bins_y, o.num_threads, sx, sy, ox, oy,
                          ratio, map);
      for (auto &v : map) {
        v = v / (dx * dy) - target;
      }
      const auto poisson_start = std::chrono::steady_clock::now();
      solver.solve(map, phi, fx, fy);
      r.poisson_seconds += std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - poisson_start)
                               .count();
      auto dg = detail::density_force(d, o.bins_x, o.bins_y, o.num_threads, sx,
                                      sy, ox, oy, ratio, fx, fy);
      if (!global_iter) {
        double wsum = 0, dsum = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          wsum += std::abs(wg.x[i]) + std::abs(wg.y[i]);
          dsum += std::abs(dg.x[i]) + std::abs(dg.y[i]);
        }
        density_weight *= std::max(wsum, 1.) / std::max(dsum, 1e-12);
      }
      if (o.macro_align_weight > 0) {
        if (global_iter % 20 == 0)
          align_pairs = detail::macro_align_pairs(d, o.macro_align_nearest);
        align_ramp = detail::macro_align_ramp(o, global_iter);
        double mean_x = 0, mean_y = 0;
        int macros = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          if (!d.node_align_class.empty() && d.node_align_class[i] >= 0 &&
              !d.node_macro.empty() && d.node_macro[i]) {
            mean_x += std::abs(wg.x[i] + density_weight * dg.x[i]);
            mean_y += std::abs(wg.y[i] + density_weight * dg.y[i]);
            ++macros;
          }
        }
        const double scale =
            macros ? o.macro_align_weight * align_ramp / macros : 0;
        const auto ag = detail::macro_alignment(
            d, align_pairs, scale * mean_x, scale * mean_y, o.macro_align_sigma,
            &align_residual);
        for (int i = 0; i < original; ++i) {
          wg.x[i] += ag.x[i];
          wg.y[i] += ag.y[i];
        }
      }
      if (o.cohesion_weight > 0 && progress > 0 && !d.cohesion_groups.empty()) {
        const auto cg = detail::cohesion(d, o.cohesion_weight * progress);
        for (int i = 0; i < original; ++i) {
          wg.x[i] += cg.x[i];
          wg.y[i] += cg.y[i];
        }
      }
      if (o.macro_overlap_weight > 0) {
        const auto mg = detail::macro_overlap(
            d, o.macro_overlap_weight *
                   (1 + o.macro_align_overlap_boost * align_ramp));
        for (int i = 0; i < original; ++i) {
          wg.x[i] += mg.x[i];
          wg.y[i] += mg.y[i];
        }
      }
      if (r.centroid.started) {
        const auto begin = std::chrono::steady_clock::now();
        const auto cg = detail::centroid_attraction(d, 1);
        detail::Gradient rg;
        if (o.centroid_repulsion_weight > 0)
          rg = detail::centroid_repulsion(d, 1, o.target_density);
        double base = 0, attraction = 0, repulsion = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          const bool repels =
              o.centroid_repulsion_weight > 0 && (rg.x[i] != 0 || rg.y[i] != 0);
          if (cg.x[i] != 0 || cg.y[i] != 0 || repels) {
            if (repels)
              repulsion += rg.x[i] * rg.x[i] + rg.y[i] * rg.y[i];
            base += std::pow(wg.x[i] + density_weight * dg.x[i], 2) +
                    std::pow(wg.y[i] + density_weight * dg.y[i], 2);
            attraction += cg.x[i] * cg.x[i] + cg.y[i] * cg.y[i];
          }
        }
        const double ramp =
            std::clamp(double(global_iter + 1 - r.centroid.start_iteration) /
                           o.centroid_ramp,
                       0., 1.);
        // Preserve the existing attraction normalization. A perfectly compact
        // group still needs repulsion, so use its norm only in that zero case.
        const double norm = attraction > 0 ? attraction : repulsion;
        const double scale =
            norm > 0 ? o.centroid_weight * ramp * ramp * std::sqrt(base / norm)
                     : 0;
        for (int i = 0; i < original; ++i) {
          wg.x[i] += scale * cg.x[i];
          wg.y[i] += scale * cg.y[i];
          if (o.centroid_repulsion_weight > 0) {
            wg.x[i] += scale * centroid_repulsion_weight * rg.x[i];
            wg.y[i] += scale * centroid_repulsion_weight * rg.y[i];
          }
        }
        r.centroid.force_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          begin)
                .count();
        r.centroid.effective_repulsion_weight = centroid_repulsion_weight;
        ++r.centroid.iterations;
      }
      double max_grad = 0, ss = 0, yy = 0;
      std::vector<double> gx(n), gy(n);
      for (int i = 0; i < n; ++i) {
        if (active(i)) {
          double degree = i < original ? d.flat_node2pin_start[i + 1] -
                                             d.flat_node2pin_start[i]
                                       : 0;
          double pre = std::max(1., degree + density_weight * d.node_size_x[i] *
                                                 d.node_size_y[i]);
          gx[i] = (wg.x[i] + density_weight * dg.x[i]) / pre;
          gy[i] = (wg.y[i] + density_weight * dg.y[i]) / pre;
          if (!finite(gx[i]) || !finite(gy[i])) {
            throw std::runtime_error("nonfinite placement gradient");
          }
          max_grad =
              std::max({max_grad, std::abs(gx[i]) / dx, std::abs(gy[i]) / dy});
          ss += std::pow(d.node_x[i] - old_ex[i], 2) +
                std::pow(d.node_y[i] - old_ey[i], 2);
          yy += std::pow(gx[i] - old_gx[i], 2) + std::pow(gy[i] - old_gy[i], 2);
        }
      }
      double cap = o.max_step / std::max(max_grad, 1e-12);
      alpha = global_iter && yy > 1e-20
                  ? std::clamp(std::sqrt(ss / yy), cap * 0.01, cap)
                  : cap;
      old_ex = d.node_x;
      old_ey = d.node_y;
      old_gx = gx;
      old_gy = gy;
      prev_x = x;
      prev_y = y;
      for (int i = 0; i < n; ++i) {
        if (active(i)) {
          x[i] = d.node_x[i] - alpha * gx[i];
          y[i] = d.node_y[i] - alpha * gy[i];
          project(i, x[i], y[i]);
        }
      }
      t = nt;
      r.iterations = iter + 1;
      if (o.observer && r.iterations % o.observer_interval == 0) {
        Design physical = design;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        o.observer(physical, PlacementStage::global_placement,
                   offset + r.iterations);
      }
      if ((global_iter + 1) % 20 == 0) {
        Design physical = design;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        r.overflow = overflow(physical, o.bins_x, o.bins_y, o.target_density,
                              o.scale_macro_density);
        const double current_hpwl = hpwl(physical);
        const auto overlaps = detail::macro_overlap_counts(physical);
        if (o.macro_align_weight > 0)
          detail::macro_alignment(physical, align_pairs, 0, 0,
                                  o.macro_align_sigma, &align_residual);
        r.history.push_back(
            {global_iter + 1, r.overflow, current_hpwl, density_weight,
             route_overflow, inflated_ratio, overlaps.first, overlaps.second,
             overflow(physical, o.bins_x, o.bins_y, o.target_density), target,
             filler_scale, align_pairs.size(), align_residual, align_ramp});
        if (separate_circles) {
          const auto circles =
              detail::centroid_circle_overlaps(physical, o.target_density);
          auto &sample = r.history.back();
          sample.centroid_overlaps = circles.overlaps;
          sample.centroid_max_overlap = circles.maximum;
          sample.centroid_repulsion_weight = centroid_repulsion_weight;
          r.centroid.circle_overlaps = circles.overlaps;
          r.centroid.max_circle_overlap = circles.maximum;
          r.centroid.circles_converged =
              circles.maximum <= o.centroid_stop_overlap;
          r.centroid.effective_repulsion_weight = centroid_repulsion_weight;
        }
        if (o.progress_observer) {
          o.progress_observer(r.history.back());
        }
        const bool settled =
            global_iter >= 39 && r.overflow <= o.stop_overflow &&
            (!o.macro_align_weight || global_iter >= o.macro_align_full);
        if (settled && (!centroids || global_iter + 1 >= o.centroid_ramp) &&
            (!separate_circles || r.centroid.circles_converged)) {
          r.converged = true;
          break;
        }
        if (separate_circles && global_iter + 1 >= o.centroid_ramp &&
            !r.centroid.circles_converged) {
          const double improvement =
              (previous_circle_overlap - r.centroid.max_circle_overlap) /
              std::max(previous_circle_overlap, 1e-12);
          if (improvement < .1)
            centroid_repulsion_weight =
                std::min(centroid_repulsion_weight * 2,
                         o.centroid_repulsion_weight * 10000);
          previous_circle_overlap = r.centroid.max_circle_overlap;
        }
        double growth = 1.15;
        if (o.adaptive_density) {
          const double improvement =
              previous_overflow > 0
                  ? (previous_overflow - r.overflow) / previous_overflow
                  : 1;
          growth = improvement < .05 ? 1.4 : 1.05;
          if (previous_hpwl > 0 && current_hpwl > previous_hpwl * 1.2) {
            growth = std::min(growth, 1.05);
          }
        }
        previous_overflow = r.overflow;
        previous_hpwl = current_hpwl;
        // If only circle separation is pending, do not keep escalating an
        // already-satisfied density constraint. Resume if density rises again.
        if (!separate_circles || r.overflow > o.stop_overflow)
          density_weight = std::min(density_weight * growth, 1e12);
      }
    }
    r.centroid.completed = r.centroid.started;
    if (state) {
      next.signature = signature;
      next.iteration = offset + r.iterations;
      next.x = x;
      next.y = y;
      next.prev_x = prev_x;
      next.prev_y = prev_y;
      next.old_gx = old_gx;
      next.old_gy = old_gy;
      next.old_ex = old_ex;
      next.old_ey = old_ey;
      next.density_weight = density_weight;
      next.previous_overflow = previous_overflow;
      next.previous_hpwl = previous_hpwl;
      next.momentum = t;
      next.step = alpha;
      next.inflation = inflation;
      next.align_pairs = align_pairs;
    }
    d.node_x = std::move(x);
    d.node_y = std::move(y);
    d.node_x.resize(original);
    d.node_y.resize(original);
    d.node_size_x.resize(original);
    d.node_size_y.resize(original);
    d.build_maps();
    if (o.observer && r.iterations % o.observer_interval != 0) {
      o.observer(d, PlacementStage::global_placement, r.iterations);
    }
  }
  r.global_placement_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - gp_start)
          .count();
  r.overflow =
      overflow(d, o.bins_x, o.bins_y, o.target_density, o.scale_macro_density);
  if (separate_circles) {
    const auto circles = detail::centroid_circle_overlaps(d, o.target_density);
    r.centroid.circle_overlaps = circles.overlaps;
    r.centroid.max_circle_overlap = circles.maximum;
    r.centroid.circles_converged = circles.maximum <= o.centroid_stop_overlap;
  }
  r.converged =
      r.overflow <= o.stop_overflow &&
      (!separate_circles || r.centroid.circles_converged) &&
      (!o.macro_align_weight ||
       (warm ? state->iteration : 0) + r.iterations > o.macro_align_full);
  if (o.legalize && d.num_movable) {
    if (o.macro_snap) {
      const auto before = d;
      const auto snap = detail::macro_align_snap(d, o.macro_snap_tolerance,
                                                 o.macro_snap_weight);
      AlignmentReport times;
      detail::legalize(d, o.bins_x, o.bins_y, &snap, &times);
      r.alignment = detail::macro_alignment_report(before, d, snap);
      r.alignment.macro_legalization_seconds = times.macro_legalization_seconds;
      r.alignment.cell_legalization_seconds = times.cell_legalization_seconds;
    } else
      detail::legalize(d, o.bins_x, o.bins_y, nullptr, &r.alignment);
    if (!is_legal(d) || !detail::kernel_is_legal(d)) {
      throw std::runtime_error("DREAMPlace legalization failed: insufficient "
                               "row space or unsupported geometry");
    }
    if (o.observer) {
      o.observer(d, PlacementStage::legalization, r.iterations);
    }
    if (o.detailed_place && !d.net_weights.empty()) {
      const auto detail_start = std::chrono::steady_clock::now();
      detail::detailed_place(d, o.bins_x, o.bins_y, o.num_threads,
                             o.macro_snap);
      r.detailed_placement_seconds =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        detail_start)
              .count();
      if (o.observer) {
        o.observer(d, PlacementStage::detailed_placement, r.iterations);
      }
    }
  }
  r.legal = is_legal(d);
  r.hpwl = hpwl(d);
  for (int i = d.num_movable; i < original; ++i) {
    if (d.node_x[i] != design.node_x[i] || d.node_y[i] != design.node_y[i]) {
      throw std::runtime_error("kernel moved a fixed node");
    }
  }
  if (o.observer) {
    o.observer(d, PlacementStage::finished, r.iterations);
  }
  if (state && next.signature) {
    *state = std::move(next);
  }
  design.node_x = std::move(d.node_x);
  design.node_y = std::move(d.node_y);
  return r;
}
} // namespace DreamPlace::native
