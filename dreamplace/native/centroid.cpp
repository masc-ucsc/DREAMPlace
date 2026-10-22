// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/placement.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
namespace DreamPlace::native::detail {
namespace {
bool repel(const CentroidGroup &a, const CentroidGroup &b) {
  if (a.repulsion_area == 0 || b.repulsion_area == 0)
    return false;
  if (a.repulsion_scope >= 0 || b.repulsion_scope >= 0)
    return a.repulsion_scope >= 0 && a.repulsion_scope == b.repulsion_scope;
  return a.family == b.family;
}
struct Moment {
  double area = 0, x = 0, y = 0, m2 = 0;
  int cells = 0, macros = 0;
  void merge(const Moment &b) {
    if (!b.area)
      return;
    const double a = area + b.area, dx = b.x - x, dy = b.y - y;
    m2 += b.m2 + (dx * dx + dy * dy) * area * b.area / a;
    x += dx * b.area / a;
    y += dy * b.area / a;
    area = a;
    cells += b.cells;
    macros += b.macros;
  }
};
Moment node(const Design &d, int n) {
  return {d.centroid_node_area[n],
          d.node_x[n] + d.node_size_x[n] / 2 - d.layout.xl,
          d.node_y[n] + d.node_size_y[n] / 2 - d.layout.yl,
          0,
          1,
          int(!d.node_macro.empty() && d.node_macro[n])};
}
std::vector<Moment> color_moments(const Design &d) {
  std::vector<Moment> result(d.centroid_color_names.size());
  for (size_t n = 0; n < d.node_centroid_color.size(); ++n)
    if (d.node_centroid_color[n] >= 0)
      result[d.node_centroid_color[n]].merge(node(d, n));
  return result;
}
Moment group_moment(const Design &d, const CentroidGroup &g,
                    const std::vector<Moment> &colors) {
  Moment result;
  for (int c : g.colors)
    result.merge(colors[c]);
  for (int n : g.nodes)
    result.merge(node(d, n));
  return result;
}
// Chain rule from a gradient with respect to each group's moving centroid back
// to its members: explicit nodes directly, colored nodes through every color.
void scatter_group_gradient(const Design &d, const std::vector<double> &gx,
                            const std::vector<double> &gy,
                            const std::vector<Moment> &moments, Gradient &out) {
  std::vector<double> cx(d.centroid_color_names.size()), cy(cx.size());
  for (size_t k = 0; k < moments.size(); ++k) {
    const double x = gx[k] / moments[k].area, y = gy[k] / moments[k].area;
    const auto &g = d.centroid_groups[k];
    for (int c : g.colors) {
      cx[c] += x;
      cy[c] += y;
    }
    for (int n : g.nodes) {
      out.x[n] += d.centroid_node_area[n] * x;
      out.y[n] += d.centroid_node_area[n] * y;
    }
  }
  for (size_t n = 0; n < d.node_centroid_color.size(); ++n) {
    const int c = d.node_centroid_color[n];
    if (c >= 0) {
      out.x[n] += d.centroid_node_area[n] * cx[c];
      out.y[n] += d.centroid_node_area[n] * cy[c];
    }
  }
}
} // namespace
std::vector<CentroidSample> centroid_summary(const Design &d) {
  auto colors = color_moments(d);
  std::vector<CentroidSample> out;
  for (const auto &g : d.centroid_groups) {
    auto m = group_moment(d, g, colors);
    out.push_back({m.x + d.layout.xl, m.y + d.layout.yl,
                   std::sqrt(std::max(0., m.m2 / m.area)), m.area, m.cells,
                   m.macros});
  }
  return out;
}
Gradient centroid_attraction(const Design &d, double lambda) {
  Gradient out;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  if (!lambda)
    return out;
  const auto moments = color_moments(d);
  std::vector<double> stiffness(moments.size()), tx(moments.size()),
      ty(moments.size());
  for (size_t k = 0; k < d.centroid_groups.size(); ++k) {
    const auto &g = d.centroid_groups[k];
    auto m = group_moment(d, g, moments);
    const double scale = lambda / m.area;
    out.value += scale * m.m2;
    if (!d.centroid_targets.empty()) {
      const double x = d.centroid_targets[k].first - d.layout.xl;
      const double y = d.centroid_targets[k].second - d.layout.yl;
      out.value += lambda * (std::pow(m.x - x, 2) + std::pow(m.y - y, 2));
      m.x = x;
      m.y = y;
    }
    for (int c : g.colors) {
      stiffness[c] += scale;
      tx[c] += scale * m.x;
      ty[c] += scale * m.y;
    }
    for (int n : g.nodes) {
      const auto p = node(d, n);
      out.x[n] += 2 * scale * p.area * (p.x - m.x);
      out.y[n] += 2 * scale * p.area * (p.y - m.y);
    }
  }
  // A color in n candidates has n terms, each aimed at its candidate's own
  // moving centroid. Aggregate coefficients, not ownership, in one node pass.
  for (size_t n = 0; n < d.node_centroid_color.size(); ++n) {
    const int c = d.node_centroid_color[n];
    if (c < 0)
      continue;
    const auto p = node(d, n);
    out.x[n] += 2 * p.area * (stiffness[c] * p.x - tx[c]);
    out.y[n] += 2 * p.area * (stiffness[c] * p.y - ty[c]);
  }
  return out;
}
namespace {
struct AttractionGroups {
  std::vector<Moment> moments;
  std::vector<double> radii;
  std::vector<std::vector<size_t>> colors;
  AttractionGroups(const Design &d, double density)
      : colors(d.centroid_color_names.size()) {
    if (!std::isfinite(density) || density <= 0 || density > 1)
      throw std::invalid_argument("invalid centroid attraction density");
    const auto cm = color_moments(d);
    for (size_t k = 0; k < d.centroid_groups.size(); ++k) {
      const auto &g = d.centroid_groups[k];
      auto m = group_moment(d, g, cm);
      radii.push_back(
          std::sqrt((g.repulsion_area > 0 ? g.repulsion_area : m.area) /
                    (std::numbers::pi * density)));
      if (!d.centroid_targets.empty()) {
        m.x = d.centroid_targets[k].first - d.layout.xl;
        m.y = d.centroid_targets[k].second - d.layout.yl;
      }
      moments.push_back(m);
      for (int c : g.colors)
        colors[c].push_back(k);
    }
  }
  template <class F> void members(const Design &d, F f) const {
    for (size_t k = 0; k < d.centroid_groups.size(); ++k)
      for (int n : d.centroid_groups[k].nodes)
        f(n, k);
    for (size_t n = 0; n < d.node_centroid_color.size(); ++n)
      if (d.node_centroid_color[n] >= 0)
        for (size_t k : colors[d.node_centroid_color[n]])
          f(n, k);
  }
};
} // namespace
std::vector<double> centroid_attraction_reference(const Design &d,
                                                  double density) {
  const AttractionGroups groups(d, density);
  std::vector<double> out(d.node_x.size());
  groups.members(d, [&](size_t n, size_t k) {
    out[n] +=
        2 * d.centroid_node_area[n] / groups.moments[k].area * groups.radii[k];
  });
  return out;
}
Gradient centroid_circle_attraction(const Design &d, double lambda,
                                    double density, double outside,
                                    bool macros_only) {
  if (macros_only && d.centroid_targets.empty())
    throw std::invalid_argument(
        "macro-only attraction requires frozen targets");
  if (!std::isfinite(lambda) || lambda < 0 || !std::isfinite(outside) ||
      outside < 0)
    throw std::invalid_argument("invalid centroid circle attraction controls");
  const AttractionGroups groups(d, density);
  Gradient out;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  if (!lambda)
    return out;
  std::vector<double> sum_x(groups.moments.size()),
      sum_y(groups.moments.size());
  groups.members(d, [&](size_t n, size_t k) {
    if (macros_only && (n >= size_t(d.num_movable) || d.node_macro.empty() ||
                        !d.node_macro[n]))
      return;
    const auto p = node(d, n);
    const auto &m = groups.moments[k];
    const double dx = p.x - m.x, dy = p.y - m.y;
    const double distance = std::hypot(dx, dy);
    const double excess = std::max(0., distance - groups.radii[k]);
    const double weight = lambda * p.area / m.area;
    out.value += weight * (distance * distance + outside * excess * excess);
    const double scale =
        2 * weight * (distance > 0 ? 1 + outside * excess / distance : 1);
    const double gx = scale * dx, gy = scale * dy;
    out.x[n] += gx;
    out.y[n] += gy;
    sum_x[k] += gx;
    sum_y[k] += gy;
  });
  if (d.centroid_targets.empty()) {
    // The centroid moves with its members. Differentiating it removes the
    // group's resultant force, so the shape spring never pins its translation.
    // Unlike the quadratic law, these terms do not cancel on their own.
    groups.members(d, [&](size_t n, size_t k) {
      const double weight = d.centroid_node_area[n] / groups.moments[k].area;
      out.x[n] -= weight * sum_x[k];
      out.y[n] -= weight * sum_y[k];
    });
  }
  return out;
}
// Explicit hierarchy edges constrain relative centroid positions. The quartic
// excess-distance energy gives a cubic restoring force outside the parent
// circle. Both endpoints move, so the resulting member forces sum to zero.
Gradient centroid_parent_attraction(const Design &d, double lambda,
                                    double density, double outside,
                                    std::uint64_t *edges,
                                    double *max_distance) {
  if (!std::isfinite(lambda) || lambda < 0 || !std::isfinite(outside) ||
      outside < 0 || !std::isfinite(density) || density <= 0 || density > 1 ||
      !d.centroid_targets.empty())
    throw std::invalid_argument(
        "invalid moving parent-centroid attraction controls");
  const auto colors = color_moments(d);
  std::vector<Moment> moments;
  for (const auto &g : d.centroid_groups)
    moments.push_back(group_moment(d, g, colors));
  std::vector<double> gx(moments.size()), gy(moments.size());
  Gradient out;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  if (edges)
    *edges = 0;
  if (max_distance)
    *max_distance = 0;
  for (size_t k = 0; k < moments.size(); ++k) {
    const int p = d.centroid_groups[k].parent;
    if (p < 0)
      continue;
    const auto &parent = d.centroid_groups[p];
    const double area =
        parent.repulsion_area < 0 ? moments[p].area : parent.repulsion_area;
    if (area <= 0)
      continue;
    const double radius = std::sqrt(area / (std::numbers::pi * density));
    const double dx = moments[k].x - moments[p].x,
                 dy = moments[k].y - moments[p].y;
    const double distance = std::hypot(dx, dy),
                 excess = std::max(0., distance - radius);
    const double relative = excess / radius;
    out.value += lambda * (distance * distance +
                           outside * excess * excess * relative * relative);
    const double scale = 2 * lambda *
                         (distance > 0 ? 1 + 2 * outside * relative * relative *
                                                 excess / distance
                                       : 1);
    gx[k] += scale * dx;
    gy[k] += scale * dy;
    gx[p] -= scale * dx;
    gy[p] -= scale * dy;
    if (edges)
      ++*edges;
    if (max_distance)
      *max_distance = std::max(*max_distance, distance / radius);
  }
  scatter_group_gradient(d, gx, gy, moments, out);
  return out;
}
// A circle reserves the anchored physical area at the caller's requested
// density, independently of the GP macro-charge scaling and halo padding.
CircleOverlap centroid_circle_overlaps(const Design &d, double density) {
  if (!std::isfinite(density) || density <= 0 || density > 1)
    throw std::invalid_argument("invalid centroid circle density");
  const auto groups = centroid_summary(d);
  CircleOverlap out;
  for (size_t i = 0; i < groups.size(); ++i)
    for (size_t j = i + 1; j < groups.size(); ++j) {
      if (!repel(d.centroid_groups[i], d.centroid_groups[j]))
        continue;
      ++out.pairs;
      const double reach =
          (std::sqrt(d.centroid_groups[i].repulsion_area < 0
                         ? groups[i].area
                         : d.centroid_groups[i].repulsion_area) +
           std::sqrt(d.centroid_groups[j].repulsion_area < 0
                         ? groups[j].area
                         : d.centroid_groups[j].repulsion_area)) /
          std::sqrt(std::numbers::pi * density);
      if (reach <= 0)
        continue;
      const double depth =
          1 - std::hypot(groups[i].x - groups[j].x, groups[i].y - groups[j].y) /
                  reach;
      if (depth > 0) {
        ++out.overlaps;
        out.maximum = std::max(out.maximum, depth);
      }
    }
  return out;
}
Gradient centroid_repulsion(const Design &d, double lambda, double density) {
  if (!std::isfinite(lambda) || lambda < 0 || !std::isfinite(density) ||
      density <= 0 || density > 1)
    throw std::invalid_argument("invalid centroid repulsion weight/density");
  Gradient out;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  if (!lambda)
    return out;
  const auto colors = color_moments(d);
  std::vector<Moment> moments;
  std::vector<double> radii, gx(d.centroid_groups.size()),
      gy(d.centroid_groups.size());
  for (const auto &g : d.centroid_groups) {
    moments.push_back(group_moment(d, g, colors));
    radii.push_back(std::sqrt(
        (g.repulsion_area < 0 ? moments.back().area : g.repulsion_area) /
        (std::numbers::pi * density)));
  }
  for (size_t i = 0; i < moments.size(); ++i) {
    for (size_t j = i + 1; j < moments.size(); ++j) {
      if (!repel(d.centroid_groups[i], d.centroid_groups[j]))
        continue;
      double dx = moments[i].x - moments[j].x;
      double dy = moments[i].y - moments[j].y;
      const double distance = std::hypot(dx, dy);
      const double depth = radii[i] + radii[j] - distance;
      if (depth <= 0)
        continue;
      if (distance > 0) {
        dx /= distance;
        dy /= distance;
      } else {
        // The radial energy has a cusp here. Choose a deterministic, equal
        // and opposite separating subgradient instead of leaving a stack.
        const std::uint64_t h =
            (i + 1) * 0x9e3779b97f4a7c15ULL ^ (j + 1) * 0xbf58476d1ce4e5b9ULL;
        const double angle = 2 * std::numbers::pi * (h % 1024) / 1024;
        dx = std::cos(angle);
        dy = std::sin(angle);
      }
      out.value += lambda * depth * depth;
      const double scale = -2 * lambda * depth;
      gx[i] += scale * dx;
      gy[i] += scale * dy;
      gx[j] -= scale * dx;
      gy[j] -= scale * dy;
    }
  }
  // Chain rule for c_g = sum(a_i p_i)/A_g. Shared colors receive the sum
  // of all their candidate contributions; no ownership or cell-pair force.
  std::vector<double> cx(colors.size()), cy(colors.size());
  for (size_t k = 0; k < moments.size(); ++k) {
    const double x = gx[k] / moments[k].area, y = gy[k] / moments[k].area;
    const auto &g = d.centroid_groups[k];
    for (int c : g.colors) {
      cx[c] += x;
      cy[c] += y;
    }
    for (int n : g.nodes) {
      out.x[n] += d.centroid_node_area[n] * x;
      out.y[n] += d.centroid_node_area[n] * y;
    }
  }
  for (size_t n = 0; n < d.node_centroid_color.size(); ++n) {
    const int c = d.node_centroid_color[n];
    if (c >= 0) {
      out.x[n] += d.centroid_node_area[n] * cx[c];
      out.y[n] += d.centroid_node_area[n] * cy[c];
    }
  }
  return out;
}
std::pair<double, double> centroid_force_reference(const Design &d,
                                                   const Gradient &wire,
                                                   const Gradient &density,
                                                   double density_weight) {
  std::vector<double> cx(d.centroid_color_names.size()), cy(cx.size());
  const auto force = [&](int n) {
    return std::pair{wire.x[n] + density_weight * density.x[n],
                     wire.y[n] + density_weight * density.y[n]};
  };
  for (int n = 0; n < d.num_movable; ++n) {
    const int c = d.node_centroid_color[n];
    if (c >= 0) {
      const auto [x, y] = force(n);
      cx[c] += x;
      cy[c] += y;
    }
  }
  double x = 0, y = 0;
  for (const auto &g : d.centroid_groups) {
    double gx = 0, gy = 0;
    for (int c : g.colors) {
      gx += cx[c];
      gy += cy[c];
    }
    for (int n : g.nodes) {
      if (n >= d.num_movable)
        continue;
      const auto [fx, fy] = force(n);
      gx += fx;
      gy += fy;
    }
    x += std::abs(gx);
    y += std::abs(gy);
  }
  const auto count = d.centroid_groups.size();
  return {count ? x / count : 0, count ? y / count : 0};
}

Gradient centroid_alignment(const Design &d, double lambda_x, double lambda_y,
                            double sigma, double density, std::uint64_t *count,
                            double *residual) {
  if (!std::isfinite(density) || density <= 0 || density > 1)
    throw std::invalid_argument("invalid centroid alignment density");
  if (!d.centroid_targets.empty())
    throw std::invalid_argument("centroid alignment requires moving centroids");
  const auto colors = color_moments(d);
  std::vector<Moment> moments;
  Design circles;
  circles.layout = d.layout;
  circles.num_movable = d.centroid_groups.size();
  std::map<std::pair<int, int>, int> classes;
  for (const auto &g : d.centroid_groups) {
    const auto m = group_moment(d, g, colors);
    moments.push_back(m);
    const double radius =
        std::sqrt((g.repulsion_area < 0 ? m.area : g.repulsion_area) /
                  (std::numbers::pi * density));
    circles.node_x.push_back(m.x + d.layout.xl - radius);
    circles.node_y.push_back(m.y + d.layout.yl - radius);
    circles.node_size_x.push_back(2 * radius);
    circles.node_size_y.push_back(2 * radius);
    circles.node_macro.push_back(1);
    const auto [it, inserted] =
        classes.try_emplace({g.family, g.parent}, classes.size());
    (void)inserted;
    circles.node_align_class.push_back(radius > 0 ? it->second : -1);
  }
  const auto pairs = macro_grid_align_pairs(circles);
  if (count)
    *count = pairs.size();
  const auto centers =
      macro_alignment(circles, pairs, lambda_x, lambda_y, sigma, residual);
  Gradient out;
  out.value = centers.value;
  out.x.assign(d.node_x.size(), 0);
  out.y.assign(d.node_y.size(), 0);
  scatter_group_gradient(d, centers.x, centers.y, moments, out);
  return out;
}

double centroid_cell_area(const Design &d) {
  double area = 0;
  int count = 0;
  const auto size = d.centroid_node_area.empty() ? d.node_x.size()
                                                 : d.centroid_node_area.size();
  for (size_t n = 0; n < size; ++n) {
    const double a = d.centroid_node_area.empty()
                         ? d.node_size_x[n] * d.node_size_y[n]
                         : d.centroid_node_area[n];
    if ((d.node_macro.empty() || !d.node_macro[n]) && a > 0) {
      area += a;
      ++count;
    }
  }
  return count ? area / count : d.layout.site_width * d.layout.row_height;
}
std::vector<CentroidGroup>
centroid_candidates(const Design &d, const std::vector<std::string> &paths,
                    const std::vector<std::string> &masters, double minimum,
                    int minimum_copies, int minimum_macros) {
  if (paths.size() != d.node_x.size() || masters.size() != paths.size() ||
      !std::isfinite(minimum) || minimum < 0 || minimum_copies < 2 ||
      minimum_macros < 0)
    throw std::invalid_argument(
        "invalid centroid candidate metadata/threshold");
  struct Prefix {
    std::vector<int> macros, nodes;
    std::set<int> colors;
    double area = 0;
    int cells = 0;
  };
  std::map<std::string, Prefix> prefixes;

  for (size_t n = 0; n < paths.size(); ++n) {
    if (d.node_macro.empty() || !d.node_macro[n])
      continue;
    if (n >= static_cast<size_t>(d.num_movable))
      continue;
    for (size_t p = paths[n].find('.'); p != std::string::npos;
         p = paths[n].find('.', p + 1))
      prefixes[paths[n].substr(0, p)].macros.push_back(n);
  }
  const double unit = centroid_cell_area(d);
  for (size_t n = 0; n < paths.size(); ++n) {
    const double a = d.centroid_node_area.empty()
                         ? d.node_size_x[n] * d.node_size_y[n]
                         : d.centroid_node_area[n];
    if (a <= 0)
      continue;
    for (size_t p = paths[n].find('.'); p != std::string::npos;
         p = paths[n].find('.', p + 1)) {
      auto it = prefixes.find(paths[n].substr(0, p));
      if (it == prefixes.end())
        continue;
      auto &v = it->second;
      v.area += a;
      ++v.cells;
      const int c =
          d.node_centroid_color.empty() ? -1 : d.node_centroid_color[n];
      if (c < 0)
        v.nodes.push_back(n);
      else
        v.colors.insert(c);
    }
  }
  // Eligibility is per immediate hierarchy parent, before alias collapse.
  // A deepest unary wrapper may be a singleton while its outer instances are
  // replicated siblings. It must not hide the eligible outer scope.
  std::map<std::string, std::map<std::string, std::vector<std::string>>>
      families;
  for (const auto &[path, v] : prefixes) {
    if (v.macros.size() <= static_cast<size_t>(minimum_macros) ||
        v.area / unit <= minimum)
      continue;
    std::vector<std::string> leaves;
    for (int n : v.macros) {
      std::ostringstream key;
      key << paths[n].substr(path.size() + 1) << ':' << masters[n] << ':'
          << std::hexfloat << d.node_size_x[n] << ':' << d.node_size_y[n];
      leaves.push_back(key.str());
    }
    std::sort(leaves.begin(), leaves.end());
    std::string signature;
    for (const auto &leaf : leaves)
      signature += std::to_string(leaf.size()) + ':' + leaf;
    const auto dot = path.rfind('.');
    const auto parent =
        dot == std::string::npos ? std::string{} : path.substr(0, dot);
    families[signature][parent].push_back(path);
  }
  std::vector<CentroidGroup> groups;
  int family = 0;
  for (const auto &[signature, parents] : families) {
    (void)signature;
    bool admitted = false;
    std::map<std::vector<int>, CentroidGroup> unique;
    for (const auto &[parent, copies] : parents) {
      if (copies.size() < static_cast<size_t>(minimum_copies))
        continue;
      for (const auto &path : copies) {
        const auto &v = prefixes.at(path);
        auto &g = unique[v.macros];
        if (path.size() <= g.name.size())
          continue;
        g.name = path;
        g.family = family;
        g.scope_macros = v.macros.size();
        g.scope_cells = v.cells;
        g.equivalent_cells = v.area / unit;
        g.nodes = v.nodes;
        g.colors.assign(v.colors.begin(), v.colors.end());
        g.copy_parent = parent;
        g.copies_per_parent = copies.size();
      }
    }
    for (auto &[ids, g] : unique) {
      (void)ids;
      groups.push_back(std::move(g));
      admitted = true;
    }
    family += admitted;
  }
  std::sort(groups.begin(), groups.end(),
            [](const auto &a, const auto &b) { return a.name < b.name; });
  for (size_t g = 0; g < groups.size(); ++g)
    for (size_t p = 0; p < g; ++p)
      if (groups[g].name.starts_with(groups[p].name + '.') &&
          (groups[g].parent < 0 ||
           groups[p].name.size() > groups[groups[g].parent].name.size()))
        groups[g].parent = p;
  return groups;
}
} // namespace DreamPlace::native::detail
