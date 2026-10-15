// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

#include "dreamplace/native/spectral.hpp"

namespace DreamPlace::native {
namespace {
void require(bool condition, const char* message) {
  if (!condition) throw std::invalid_argument(message);
}
bool finite(double v) { return std::isfinite(v); }
bool aligned(double v) { return std::abs(v - std::round(v)) <= 1e-6; }
int count(size_t n) {
  require(n <= static_cast<size_t>(std::numeric_limits<int>::max()), "array exceeds 32-bit kernel indices");
  return static_cast<int>(n);
}
void csr_check(const std::vector<int>& flat, const std::vector<int>& start, int objects, int pins,
               const std::vector<int>& inverse) {
  require(start.size() == static_cast<size_t>(objects) + 1 && start.front() == 0 && start.back() == pins,
          "invalid CSR start/end");
  require(flat.size() == static_cast<size_t>(pins), "CSR must contain every pin exactly once");
  std::vector<bool> seen(pins, false);
  for (int n = 0; n < objects; ++n) {
    require(start[n] <= start[n + 1] && start[n] >= 0 && start[n + 1] <= pins, "nonmonotonic CSR");
    for (int j = start[n]; j < start[n + 1]; ++j) {
      int p = flat[j];
      require(p >= 0 && p < pins && !seen[p], "invalid or duplicate CSR pin");
      seen[p] = true;
      require(inverse[p] == n, "CSR disagrees with inverse map");
    }
  }
}

// Exact physical overflow, excluding fillers; fixed area reduces bin capacity.
double overflow(const Design& d, int bx, int by, double target_density) {
  const auto& l = d.layout;
  double dx = (l.xh - l.xl) / bx, dy = (l.yh - l.yl) / by;
  std::vector<double> movable(bx * by), fixed(bx * by);
  double total = 0;
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i], h = d.node_size_y[i];
    if (i < d.num_movable) total += w * h;
    int x0 = std::clamp(static_cast<int>(std::floor((x - l.xl) / dx)), 0, bx);
    int x1 = std::clamp(static_cast<int>(std::ceil((x + w - l.xl) / dx)), 0, bx);
    int y0 = std::clamp(static_cast<int>(std::floor((y - l.yl) / dy)), 0, by);
    int y1 = std::clamp(static_cast<int>(std::ceil((y + h - l.yl) / dy)), 0, by);
    auto& map = i < d.num_movable ? movable : fixed;
    for (int k = x0; k < x1; ++k)
      for (int j = y0; j < y1; ++j) {
        double a = std::max(0., std::min(x + w, l.xl + (k + 1) * dx) - std::max(x, l.xl + k * dx));
        double b = std::max(0., std::min(y + h, l.yl + (j + 1) * dy) - std::max(y, l.yl + j * dy));
        map[k * by + j] += a * b;
      }
  }
  double excess = 0;
  for (int i = 0; i < bx * by; ++i)
    excess += std::max(0., movable[i] - target_density * std::max(0., dx * dy - fixed[i]));
  return total ? excess / total : 0;
}

void options_check(const Options& o) {
  require(o.max_iterations >= 0 && o.bins_x >= 2 && o.bins_x <= 1024 && o.bins_y >= 2 && o.bins_y <= 1024,
          "invalid iteration limit or bin count (2..1024)");
  require(o.num_threads > 0 && o.num_threads <= 128, "thread count must be 1..128");
  require(finite(o.target_density) && o.target_density > 0 && o.target_density <= 1, "density must be in (0,1]");
  require(finite(o.stop_overflow) && o.stop_overflow >= 0 && o.stop_overflow <= 1, "invalid overflow threshold");
  require(finite(o.initial_density_weight) && o.initial_density_weight > 0, "invalid density weight");
  require(finite(o.max_step) && o.max_step > 0 && o.max_step <= 1, "max_step must be in (0,1]");
  require(!o.detailed_place || o.legalize, "detailed placement requires legalization");
  require(!o.observer || o.observer_interval > 0, "observer interval must be positive");
}
}  // namespace

void Design::build_maps() {
  int nodes = count(node_x.size()), pins = count(pin2node.size());
  require(!flat_net2pin_start.empty(), "missing net CSR");
  int nets = count(flat_net2pin_start.size() - 1);
  // Build into temporary arrays so a malformed design does not leave partial maps.
  std::vector<int> p2n(pins, -1);
  require(flat_net2pin.size() == static_cast<size_t>(pins) && flat_net2pin_start.front() == 0 &&
              flat_net2pin_start.back() == pins,
          "invalid net CSR");
  for (int n = 0; n < nets; ++n) {
    int a = flat_net2pin_start[n], b = flat_net2pin_start[n + 1];
    require(a >= 0 && a <= b && b <= pins, "nonmonotonic net CSR");
    for (int j = a; j < b; ++j) {
      int p = flat_net2pin[j];
      require(p >= 0 && p < pins && p2n[p] == -1, "invalid or duplicate net pin");
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
  for (int p = 0; p < pins; ++p) flat[next[pin2node[p]]++] = p;
  if (net_weights.empty()) net_weights.assign(nets, 1);
  if (net_mask.empty()) {
    net_mask.resize(nets);
    for (int n = 0; n < nets; ++n) net_mask[n] = flat_net2pin_start[n + 1] - flat_net2pin_start[n] > 1;
  }
  pin2net = std::move(p2n);
  flat_node2pin_start = std::move(start);
  flat_node2pin = std::move(flat);
}

void Design::validate() const {
  int n = count(node_x.size()), p = count(pin2node.size()), nets = count(net_weights.size());
  require(num_movable >= 0 && num_movable <= n, "invalid movable count");
  const auto& l = layout;
  require(finite(l.xl) && finite(l.yl) && finite(l.xh) && finite(l.yh) && l.xh > l.xl && l.yh > l.yl &&
              finite(l.xh - l.xl) && finite(l.yh - l.yl),
          "invalid core bounds");
  require(finite(l.site_width) && l.site_width > 0 && finite(l.row_height) && l.row_height > 0,
          "invalid site/row pitch");
  require((l.xh - l.xl) / l.site_width < std::numeric_limits<int>::max() &&
              (l.yh - l.yl) / l.row_height < std::numeric_limits<int>::max(),
          "layout exceeds kernel site indices");
  require(aligned((l.xh - l.xl) / l.site_width) && aligned((l.yh - l.yl) / l.row_height),
          "core must contain whole sites/rows");
  require(node_y.size() == node_x.size() && node_size_x.size() == node_x.size() && node_size_y.size() == node_x.size(),
          "node array size mismatch");
  for (int i = 0; i < n; ++i) {
    double w = node_size_x[i], h = node_size_y[i];
    require(finite(node_x[i]) && finite(node_y[i]) && finite(w) && finite(h) && w >= 0 && h >= 0 && finite(w * h),
            "nonfinite position or invalid node size");
    require(finite(node_x[i] + w) && finite(node_y[i] + h), "node bounds overflow");
    if (i < num_movable) {
      require(w > 0 && h > 0 && w <= l.xh - l.xl && h <= l.yh - l.yl, "movable node does not fit core");
      require(aligned(w / l.site_width) && aligned(h / l.row_height), "movable size must be integral sites/rows");
    }
  }
  require(pin2net.size() == pin2node.size() && pin_offset_x.size() == pin2node.size() &&
              pin_offset_y.size() == pin2node.size() && net_mask.size() == net_weights.size(),
          "pin/net array size mismatch");
  for (int i = 0; i < p; ++i) {
    require(pin2node[i] >= 0 && pin2node[i] < n && pin2net[i] >= 0 && pin2net[i] < nets, "invalid pin map index");
    require(finite(pin_offset_x[i]) && finite(pin_offset_y[i]), "nonfinite pin offset");
  }
  for (int i = 0; i < nets; ++i) {
    require(finite(net_weights[i]) && net_weights[i] >= 0 && net_mask[i] <= 1, "invalid net weight/mask");
  }
  csr_check(flat_net2pin, flat_net2pin_start, nets, p, pin2net);
  csr_check(flat_node2pin, flat_node2pin_start, n, p, pin2node);
  for (int i = 0; i < nets; ++i)
    require(!net_mask[i] || flat_net2pin_start[i + 1] - flat_net2pin_start[i] >= 2,
            "active net needs at least two pins");
}

double hpwl(const Design& d) {
  d.validate();
  double wl = 0;
  for (int n = 0; n < count(d.net_weights.size()); ++n) {
    if (!d.net_mask[n]) continue;
    double xmin = std::numeric_limits<double>::infinity(), ymin = xmin, xmax = -xmin, ymax = -xmin;
    for (int j = d.flat_net2pin_start[n]; j < d.flat_net2pin_start[n + 1]; ++j) {
      int p = d.flat_net2pin[j], i = d.pin2node[p];
      double x = d.node_x[i] + d.pin_offset_x[p], y = d.node_y[i] + d.pin_offset_y[p];
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    wl += d.net_weights[n] * (xmax - xmin + ymax - ymin);
  }
  return wl;
}

bool is_legal(const Design& d) {
  d.validate();
  const auto& l = d.layout;
  constexpr double eps = 1e-6;
  struct Interval {
    double x, end;
    int id;
  };
  std::vector<std::vector<Interval>> rows(static_cast<size_t>(std::round((l.yh - l.yl) / l.row_height)));
  for (int i = 0; i < count(d.node_x.size()); ++i) {
    double x = d.node_x[i], y = d.node_y[i], w = d.node_size_x[i], h = d.node_size_y[i];
    if (!w || !h) continue;
    if (i < d.num_movable && (x < l.xl - eps || y < l.yl - eps || x + w > l.xh + eps || y + h > l.yh + eps ||
                              !aligned((x - l.xl) / l.site_width) || !aligned((y - l.yl) / l.row_height)))
      return false;
    int first = std::max(0., std::floor((y - l.yl + eps) / l.row_height));
    int last = std::min<double>(rows.size(), std::ceil((y + h - l.yl - eps) / l.row_height));
    for (int r = first; r < last; ++r) rows[r].push_back({x, x + w, i});
  }
  for (auto& row : rows) {
    std::sort(row.begin(), row.end(), [](auto a, auto b) { return a.x < b.x; });
    double end_all = -std::numeric_limits<double>::infinity(), end_movable = end_all;
    for (auto v : row) {
      if (v.id < d.num_movable ? v.x < end_all - eps : v.x < end_movable - eps) return false;
      end_all = std::max(end_all, v.end);
      if (v.id < d.num_movable) end_movable = std::max(end_movable, v.end);
    }
  }
  return true;
}

const char* stage_name(PlacementStage stage) {
  switch (stage) {
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

Result place(Design& design, const Options& o) {
  design.validate();
  options_check(o);
  Result r;
  r.initial_hpwl = hpwl(design);
  Design d = design;
  const int original = count(d.node_x.size());
  if (o.observer) o.observer(design, PlacementStage::initial, 0);
  const auto& l = d.layout;
  const double width = l.xh - l.xl, height = l.yh - l.yl;
  double dx = width / o.bins_x, dy = height / o.bins_y;
  // A local engine owns all temporary/filler state. No global RNG or settings.
  std::mt19937 rng(1000);
  std::uniform_real_distribution<double> unit(0, 1);
  double area = 0, fixed_area = 0;
  for (int i = 0; i < original; ++i) {
    double a = d.node_size_x[i] * d.node_size_y[i];
    if (i < d.num_movable)
      area += a;
    else {
      double w = std::max(0., std::min(l.xh, d.node_x[i] + d.node_size_x[i]) - std::max(l.xl, d.node_x[i]));
      double h = std::max(0., std::min(l.yh, d.node_y[i] + d.node_size_y[i]) - std::max(l.yl, d.node_y[i]));
      fixed_area += w * h;
    }
  }
  if (d.num_movable && o.max_iterations) {
    auto active = [&](int i) { return i < d.num_movable || i >= original; };
    auto project = [&](int i, double& x, double& y) {
      x = std::clamp(x, l.xl, l.xh - d.node_size_x[i]);
      y = std::clamp(y, l.yl, l.yh - d.node_size_y[i]);
    };
    for (int i = 0; i < d.num_movable; ++i) project(i, d.node_x[i], d.node_y[i]);
    if (o.add_fillers) {
      double target_area = std::max(0., o.target_density * (width * height - fixed_area) - area);
      double fw = area / d.num_movable / l.row_height, fh = l.row_height;
      fw = std::clamp(fw, l.site_width, width);
      double desired = std::floor(target_area / (fw * fh));
      require(desired <= 1000000, "too many fillers; reduce core whitespace or disable fillers");
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
    std::vector<double> sx(n), sy(n), ox(n), oy(n), ratio(n);
    for (int i = 0; i < n; ++i) {
      sx[i] = active(i) ? std::max(d.node_size_x[i], std::sqrt(2.) * dx) : d.node_size_x[i];
      sy[i] = active(i) ? std::max(d.node_size_y[i], std::sqrt(2.) * dy) : d.node_size_y[i];
      ox[i] = (d.node_size_x[i] - sx[i]) / 2;
      oy[i] = (d.node_size_y[i] - sy[i]) / 2;
      ratio[i] = sx[i] * sy[i] ? d.node_size_x[i] * d.node_size_y[i] / (sx[i] * sy[i]) : 0;
      if (!active(i)) ratio[i] *= o.target_density;
    }
    PoissonSolver solver(o.bins_x, o.bins_y, width, height);
    std::vector<double> map, phi, fx, fy;
    auto x = d.node_x, y = d.node_y, prev_x = x, prev_y = y;
    std::vector<double> old_gx(n), old_gy(n), old_ex(n), old_ey(n);
    double density_weight = o.initial_density_weight, t = 1, alpha = 0;
    for (int iter = 0; iter < o.max_iterations; ++iter) {
      double progress = double(iter) / std::max(1, o.max_iterations - 1);
      double gamma = std::max(l.site_width, 4 * std::max(dx, dy) * std::pow(0.05, progress));
      double nt = (1 + std::sqrt(1 + 4 * t * t)) / 2, momentum = (t - 1) / nt;
      for (int i = 0; i < n; ++i)
        if (active(i)) {
          d.node_x[i] = x[i] + momentum * (x[i] - prev_x[i]);
          d.node_y[i] = y[i] + momentum * (y[i] - prev_y[i]);
          project(i, d.node_x[i], d.node_y[i]);
        }
      auto wg = detail::wirelength(d, gamma, o.num_threads);
      detail::density_map(d, o.bins_x, o.bins_y, o.num_threads, sx, sy, ox, oy, ratio, map);
      for (auto& v : map) v = v / (dx * dy) - o.target_density;
      solver.solve(map, phi, fx, fy);
      auto dg = detail::density_force(d, o.bins_x, o.bins_y, o.num_threads, sx, sy, ox, oy, ratio, fx, fy);
      if (!iter) {
        double wsum = 0, dsum = 0;
        for (int i = 0; i < d.num_movable; ++i) {
          wsum += std::abs(wg.x[i]) + std::abs(wg.y[i]);
          dsum += std::abs(dg.x[i]) + std::abs(dg.y[i]);
        }
        density_weight *= std::max(wsum, 1.) / std::max(dsum, 1e-12);
      }
      double max_grad = 0, ss = 0, yy = 0;
      std::vector<double> gx(n), gy(n);
      for (int i = 0; i < n; ++i)
        if (active(i)) {
          double degree = i < original ? d.flat_node2pin_start[i + 1] - d.flat_node2pin_start[i] : 0;
          double pre = std::max(1., degree + density_weight * d.node_size_x[i] * d.node_size_y[i]);
          gx[i] = (wg.x[i] + density_weight * dg.x[i]) / pre;
          gy[i] = (wg.y[i] + density_weight * dg.y[i]) / pre;
          if (!finite(gx[i]) || !finite(gy[i])) throw std::runtime_error("nonfinite placement gradient");
          max_grad = std::max({max_grad, std::abs(gx[i]) / dx, std::abs(gy[i]) / dy});
          ss += std::pow(d.node_x[i] - old_ex[i], 2) + std::pow(d.node_y[i] - old_ey[i], 2);
          yy += std::pow(gx[i] - old_gx[i], 2) + std::pow(gy[i] - old_gy[i], 2);
        }
      double cap = o.max_step / std::max(max_grad, 1e-12);
      alpha = iter && yy > 1e-20 ? std::clamp(std::sqrt(ss / yy), cap * 0.01, cap) : cap;
      old_ex = d.node_x;
      old_ey = d.node_y;
      old_gx = gx;
      old_gy = gy;
      prev_x = x;
      prev_y = y;
      for (int i = 0; i < n; ++i)
        if (active(i)) {
          x[i] = d.node_x[i] - alpha * gx[i];
          y[i] = d.node_y[i] - alpha * gy[i];
          project(i, x[i], y[i]);
        }
      t = nt;
      r.iterations = iter + 1;
      if (o.observer && r.iterations % o.observer_interval == 0) {
        Design physical = design;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        o.observer(physical, PlacementStage::global_placement, r.iterations);
      }
      if ((iter + 1) % 20 == 0) {
        Design physical = design;
        std::copy_n(x.begin(), original, physical.node_x.begin());
        std::copy_n(y.begin(), original, physical.node_y.begin());
        r.overflow = overflow(physical, o.bins_x, o.bins_y, o.target_density);
        if (iter >= 39 && r.overflow <= o.stop_overflow) {
          r.converged = true;
          break;
        }
        density_weight = std::min(density_weight * 1.15, 1e12);
      }
    }
    d.node_x = std::move(x);
    d.node_y = std::move(y);
    d.node_x.resize(original);
    d.node_y.resize(original);
    d.node_size_x.resize(original);
    d.node_size_y.resize(original);
    d.build_maps();
    if (o.observer && r.iterations % o.observer_interval != 0)
      o.observer(d, PlacementStage::global_placement, r.iterations);
  }
  r.overflow = overflow(d, o.bins_x, o.bins_y, o.target_density);
  r.converged = r.overflow <= o.stop_overflow;
  if (o.legalize && d.num_movable) {
    detail::legalize(d, o.bins_x, o.bins_y);
    if (!is_legal(d) || !detail::kernel_is_legal(d))
      throw std::runtime_error("DREAMPlace legalization failed: insufficient row space or unsupported geometry");
    if (o.observer) o.observer(d, PlacementStage::legalization, r.iterations);
    if (o.detailed_place && !d.net_weights.empty()) {
      detail::detailed_place(d, o.bins_x, o.bins_y, o.num_threads);
      if (o.observer) o.observer(d, PlacementStage::detailed_placement, r.iterations);
    }
  }
  r.legal = is_legal(d);
  r.hpwl = hpwl(d);
  for (int i = d.num_movable; i < original; ++i)
    if (d.node_x[i] != design.node_x[i] || d.node_y[i] != design.node_y[i])
      throw std::runtime_error("kernel moved a fixed node");
  if (o.observer) o.observer(d, PlacementStage::finished, r.iterations);
  design.node_x = std::move(d.node_x);
  design.node_y = std::move(d.node_y);
  return r;
}
}  // namespace DreamPlace::native
