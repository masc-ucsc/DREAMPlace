// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/floorplan.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
namespace DreamPlace::native {
namespace {
struct Linear {
  std::vector<std::pair<int, double>> terms;
  double bound = 0;
  bool equality = false;
};
Linear linear(std::initializer_list<std::pair<int, double>> terms,
              double bound = 0, bool equality = true) {
  std::map<int, double> merged;
  for (auto [i, a] : terms)
    merged[i] += a;
  Linear out{{}, bound, equality};
  for (auto [i, a] : merged)
    if (a)
      out.terms.emplace_back(i, a);
  return out;
}
double dot(const Linear &c, const std::vector<double> &x) {
  double v = 0;
  for (auto [i, a] : c.terms)
    v += a * x[i];
  return v;
}
// Dykstra projects the original coordinates onto an intersection of linear
// equalities/halfspaces in the area-weighted metric. Corrections are scalar
// because each projection is along W^-1 a. All scratch state is call-local.
bool fit(const std::vector<double> &origin, const std::vector<double> &weights,
         const std::vector<Linear> &constraints, double tolerance,
         std::vector<double> &x) {
  x = origin;
  std::vector<double> correction(constraints.size()), norm;
  for (const auto &c : constraints) {
    double sum = 0;
    for (auto [i, a] : c.terms)
      sum += a * a / weights[i];
    norm.push_back(sum);
  }
  for (int iteration = 0; iteration < 1200; ++iteration) {
    double movement = 0;
    for (size_t k = 0; k < constraints.size(); ++k) {
      if (!norm[k])
        continue;
      const auto &c = constraints[k];
      const double residual = dot(c, x) + correction[k] * norm[k] - c.bound;
      const double next =
          (c.equality ? residual : std::max(0., residual)) / norm[k];
      const double delta = correction[k] - next;
      for (auto [i, a] : c.terms) {
        const double step = delta * a / weights[i];
        x[i] += step;
        movement = std::max(movement, std::abs(step));
      }
      correction[k] = next;
    }
    if (iteration % 10 == 0 && movement < tolerance) {
      bool valid = true;
      for (const auto &c : constraints) {
        const double error = dot(c, x) - c.bound;
        valid &= (c.equality ? std::abs(error) : error) <= tolerance * 4;
      }
      if (valid)
        return true;
    }
  }
  return false; // bounded solve: never claim infeasible or unconverged fit
                // exact
}
bool ancestor(const std::vector<FloorplanCircle> &v, int a, int b) {
  for (int p = v[b].parent; p >= 0; p = v[p].parent)
    if (p == a)
      return true;
  return false;
}
struct Run {
  std::vector<int> members;
  int axis;
};
std::vector<Run> runs(const std::vector<FloorplanCircle> &v,
                      const std::vector<double> &x,
                      const std::vector<int> &members, int axis,
                      double threshold) {
  std::vector<int> order = members;
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    return std::tie(x[2 * a + axis], a) < std::tie(x[2 * b + axis], b);
  });
  std::vector<Run> out;
  for (size_t begin = 0; begin < order.size();) {
    size_t end = begin + 1;
    double r = v[order[begin]].radius;
    while (end < order.size()) {
      r = std::min(r, v[order[end]].radius);
      // Bounded total span prevents single-link chains collapsing into a line.
      if (x[2 * order[end] + axis] - x[2 * order[begin] + axis] >=
          threshold * r)
        break;
      ++end;
    }
    out.push_back(
        {std::vector<int>(order.begin() + begin, order.begin() + end), axis});
    begin = end;
  }
  return out;
}
struct Gap {
  int a, b, axis;
  double value;
};
// Penalized exact 1-D segmentation of sorted observed neighbor gaps. This
// solves the scalar grouping subproblem, not the joint nonconvex layout.
std::vector<std::vector<Gap>> gap_classes(std::vector<Gap> gaps,
                                          const FloorplanOptions &o) {
  std::stable_sort(gaps.begin(), gaps.end(), [](auto a, auto b) {
    return std::tie(a.value, a.a, a.b) < std::tie(b.value, b.a, b.b);
  });
  const int n = gaps.size();
  std::vector<double> sum(n + 1), square(n + 1),
      cost(n + 1, std::numeric_limits<double>::infinity());
  std::vector<int> previous(n + 1);
  cost[0] = 0;
  for (int i = 0; i < n; ++i) {
    sum[i + 1] = sum[i] + gaps[i].value;
    square[i + 1] = square[i] + gaps[i].value * gaps[i].value;
    for (int j = i; j >= 0; --j) {
      const double mean = (sum[i + 1] - sum[j]) / (i + 1 - j);
      if (gaps[i].value - gaps[j].value > o.spacing_tolerance * mean)
        break;
      const double variance =
          std::max(0., square[i + 1] - square[j] - (i + 1 - j) * mean * mean);
      const double candidate =
          cost[j] + variance / (mean * mean) + o.spacing_penalty;
      if (candidate < cost[i + 1]) {
        cost[i + 1] = candidate;
        previous[i + 1] = j;
      }
    }
  }
  std::vector<std::vector<Gap>> out;
  for (int i = n; i;) {
    const int j = previous[i];
    out.emplace_back(gaps.begin() + j, gaps.begin() + i);
    i = j;
  }
  return out;
}
} // namespace
FloorplanSnap snap_floorplan(const std::vector<FloorplanCircle> &circles,
                             double xl, double yl, double xh, double yh,
                             const FloorplanOptions &o) {
  const auto begin = std::chrono::steady_clock::now();
  if (!std::isfinite(xl) || !std::isfinite(yl) || !std::isfinite(xh) ||
      !std::isfinite(yh) || xh <= xl || yh <= yl)
    throw std::invalid_argument("invalid floorplan snap bounds");
  for (double v : {o.alignment_tolerance, o.spacing_tolerance,
                   o.spacing_penalty, o.symmetry_tolerance, o.max_displacement})
    if (!std::isfinite(v) || v < 0)
      throw std::invalid_argument("invalid floorplan snap control");
  for (size_t i = 0; i < circles.size(); ++i) {
    const auto &c = circles[i];
    if (!std::isfinite(c.x) || !std::isfinite(c.y) ||
        !std::isfinite(c.radius) || !std::isfinite(c.area) || c.radius <= 0 ||
        c.area <= 0 || c.family < 0 || c.parent < -1 ||
        c.parent >= static_cast<int>(i))
      throw std::invalid_argument("invalid floorplan circle");
  }
  FloorplanSnap out;
  out.circles = circles;
  std::map<int, std::vector<int>> families;
  std::vector<int> depths(circles.size());
  for (size_t i = 0; i < circles.size(); ++i) {
    if (circles[i].parent >= 0)
      depths[i] = depths[circles[i].parent] + 1;
    families[circles[i].family].push_back(i);
  }
  std::vector<int> order;
  for (auto &[f, ids] : families)
    order.push_back(f);
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    auto depth = [&](int f) {
      int d = std::numeric_limits<int>::max();
      for (int i : families[f])
        d = std::min(d, depths[i]);
      return d;
    };
    return std::pair(depth(a), a) < std::pair(depth(b), b);
  });
  for (int family : order) {
    const auto &ids = families[family];
    std::vector<FloorplanCircle> local;
    std::vector<double> origin, weights;
    double unit = 0;
    for (int i : ids) {
      local.push_back(out.circles[i]);
      origin.insert(origin.end(), {out.circles[i].x, out.circles[i].y});
      weights.insert(weights.end(), {circles[i].area, circles[i].area});
      unit += circles[i].radius;
    }
    unit /= ids.size();
    for (double &w : weights)
      w /= local.front().area;
    const double tolerance = std::max(1e-9, unit * 1e-8);
    std::vector<Linear> hard;
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
      const double reach = o.max_displacement * local[i].radius;
      const double margin = o.contain_circles ? local[i].radius : 0;
      // Boundary repair takes priority over the local regularity movement cap,
      // including when an ancestor translation moved this child outside.
      const bool fits = 2 * margin <= std::min(xh - xl, yh - yl);
      const double anchor_x =
          o.contain_circles && fits
              ? std::clamp(origin[2 * i], xl + margin, xh - margin)
              : origin[2 * i];
      const double anchor_y =
          o.contain_circles && fits
              ? std::clamp(origin[2 * i + 1], yl + margin, yh - margin)
              : origin[2 * i + 1];
      hard.push_back(linear({{2 * i, -1}},
                            -std::max(xl + margin, anchor_x - reach), false));
      hard.push_back(
          linear({{2 * i, 1}}, std::min(xh - margin, anchor_x + reach), false));
      hard.push_back(linear({{2 * i + 1, -1}},
                            -std::max(yl + margin, anchor_y - reach), false));
      hard.push_back(linear({{2 * i + 1, 1}},
                            std::min(yh - margin, anchor_y + reach), false));
      for (int j = 0; j < i; ++j) {
        if (ancestor(circles, ids[i], ids[j]) ||
            ancestor(circles, ids[j], ids[i]))
          continue;
        double dx = origin[2 * i] - origin[2 * j],
               dy = origin[2 * i + 1] - origin[2 * j + 1];
        const double distance = std::hypot(dx, dy);
        if (distance) {
          dx /= distance;
          dy /= distance;
        } else {
          const double angle = (i * 17 + j * 31) * 2.399963229728653;
          dx = std::cos(angle);
          dy = std::sin(angle);
        }
        hard.push_back(linear(
            {{2 * j, dx}, {2 * j + 1, dy}, {2 * i, -dx}, {2 * i + 1, -dy}},
            -(local[i].radius + local[j].radius), false));
      }
    }
    std::vector<double> x;
    if (!fit(origin, weights, hard, tolerance, x)) {
      out.feasible = false;
      ++out.infeasible_families;
      continue;
    }
    auto accept = [&](std::vector<Linear> extra, FloorplanRelation relation) {
      ++out.attempted;
      const size_t old = hard.size();
      hard.insert(hard.end(), extra.begin(), extra.end());
      std::vector<double> trial;
      if (!fit(origin, weights, hard, tolerance, trial)) {
        hard.resize(old);
        ++out.rejected;
        return false;
      }
      // Later midpoint/gap equalities must not bypass bounded-span admission
      // by transitively collapsing a previously rejected alignment chain.
      for (size_t i = 0; i < local.size(); ++i)
        for (size_t j = 0; j < i; ++j)
          if (local[i].parent == local[j].parent)
            for (int axis = 0; axis < 2; ++axis)
              if (std::abs(trial[2 * i + axis] - trial[2 * j + axis]) <
                      tolerance * 8 &&
                  std::abs(origin[2 * i + axis] - origin[2 * j + axis]) >=
                      o.alignment_tolerance *
                          std::min(local[i].radius, local[j].radius)) {
                hard.resize(old);
                ++out.rejected;
                return false;
              }
      x = std::move(trial);
      for (int &i : relation.nodes)
        i = ids[i];
      out.relations.push_back(std::move(relation));
      return true;
    };
    // Discover regularity only among siblings of one candidate. Keep the
    // family-wide hard circle constraints, including across parent scopes.
    std::map<int, std::vector<int>> scopes;
    for (size_t i = 0; i < local.size(); ++i)
      scopes[local[i].parent].push_back(i);
    for (const auto &[parent, members] : scopes) {
      std::vector<Run> aligned;
      for (int axis = 0; axis < 2; ++axis) {
        auto groups = runs(local, x, members, axis, o.alignment_tolerance);
        for (auto &g : groups) {
          std::vector<Linear> extra;
          for (size_t k = 1; k < g.members.size(); ++k)
            extra.push_back(linear(
                {{2 * g.members[0] + axis, 1}, {2 * g.members[k] + axis, -1}}));
          if (extra.empty() || accept(extra, {"alignment", axis, g.members, 0}))
            aligned.push_back(g);
          else
            for (int m : g.members)
              aligned.push_back({{m}, axis});
        }
      }
      for (int axis = 0; axis < 2; ++axis) {
        std::vector<Gap> gaps;
        for (auto &g : aligned)
          if (g.axis != axis) {
            std::stable_sort(g.members.begin(), g.members.end(),
                             [&](int a, int b) {
                               return std::pair(x[2 * a + axis], a) <
                                      std::pair(x[2 * b + axis], b);
                             });
            for (size_t k = 1; k < g.members.size(); ++k) {
              int a = g.members[k - 1], b = g.members[k];
              gaps.push_back({a, b, axis, x[2 * b + axis] - x[2 * a + axis]});
            }
          }
        for (const auto &cluster : gap_classes(gaps, o))
          if (cluster.size() > 1) {
            const auto base = cluster.front();
            std::vector<Linear> extra;
            std::vector<int> nodes;
            for (auto g : cluster) {
              nodes.insert(nodes.end(), {g.a, g.b});
              if (g.a != base.a || g.b != base.b)
                extra.push_back(linear({{2 * g.b + axis, 1},
                                        {2 * g.a + axis, -1},
                                        {2 * base.b + axis, -1},
                                        {2 * base.a + axis, 1}}));
            }
            accept(extra, {"spacing", axis, nodes, 0});
          }
        struct Mid {
          int a, b;
          double value, radius;
        };
        std::vector<Mid> mids;
        for (const auto &g : aligned)
          if (g.axis != axis) {
            const int a = g.members.front(), b = g.members.back();
            mids.push_back({a, b, (x[2 * a + axis] + x[2 * b + axis]) / 2,
                            std::min(local[a].radius, local[b].radius)});
          }
        std::stable_sort(mids.begin(), mids.end(),
                         [](auto a, auto b) { return a.value < b.value; });
        for (size_t first = 0; first < mids.size();) {
          size_t last = first + 1;
          double r = mids[first].radius;
          while (last < mids.size()) {
            r = std::min(r, mids[last].radius);
            if (mids[last].value - mids[first].value >=
                o.symmetry_tolerance * r)
              break;
            ++last;
          }
          if (last - first > 1) {
            std::vector<Linear> extra;
            std::vector<int> nodes;
            const auto a = mids[first];
            for (size_t k = first; k < last; ++k) {
              const auto b = mids[k];
              nodes.insert(nodes.end(), {b.a, b.b});
              if (k > first)
                extra.push_back(linear({{2 * b.a + axis, .5},
                                        {2 * b.b + axis, .5},
                                        {2 * a.a + axis, -.5},
                                        {2 * a.b + axis, -.5}}));
            }
            accept(extra, {"symmetry", axis, nodes, 0});
          }
          first = last;
        }
      }
    } // sibling scope
    for (size_t k = 0; k < ids.size(); ++k) {
      const int i = ids[k];
      // std::clamp is undefined when lo > hi; a circle wider than the core
      // keeps its fitted center and is reported through the boundary violation.
      if (o.contain_circles &&
          2 * circles[i].radius <= std::min(xh - xl, yh - yl)) {
        x[2 * k] = std::clamp(x[2 * k], xl + circles[i].radius,
                              xh - circles[i].radius);
        x[2 * k + 1] = std::clamp(x[2 * k + 1], yl + circles[i].radius,
                                  yh - circles[i].radius);
      }
      const double dx = x[2 * k] - out.circles[i].x,
                   dy = x[2 * k + 1] - out.circles[i].y;
      out.circles[i].x = x[2 * k];
      out.circles[i].y = x[2 * k + 1];
      // Parent target translations are inherited before child-family fitting.
      for (size_t j = i + 1; j < circles.size(); ++j)
        if (ancestor(circles, i, j)) {
          out.circles[j].x += dx;
          out.circles[j].y += dy;
        }
    }
  }
  for (auto &r : out.relations) {
    auto coordinate = [&](int n) {
      return r.axis ? out.circles[n].y : out.circles[n].x;
    };
    if (r.kind == "alignment")
      r.value = coordinate(r.nodes.front());
    else if (r.kind == "spacing")
      r.value = coordinate(r.nodes[1]) - coordinate(r.nodes[0]);
    else
      r.value = (coordinate(r.nodes[0]) + coordinate(r.nodes[1])) / 2;
  }
  for (size_t i = 0; i < circles.size(); ++i) {
    const auto &c = out.circles[i];
    const double margin = o.contain_circles ? c.radius : 0;
    out.maximum_boundary_violation =
        std::max({out.maximum_boundary_violation, xl + margin - c.x,
                  c.x + margin - xh, yl + margin - c.y, c.y + margin - yh});
    out.maximum_displacement = std::max(
        out.maximum_displacement, std::hypot(out.circles[i].x - circles[i].x,
                                             out.circles[i].y - circles[i].y));
    for (size_t j = 0; j < i; ++j)
      if (circles[i].family == circles[j].family && !ancestor(circles, i, j) &&
          !ancestor(circles, j, i)) {
        const double depth =
            1 - std::hypot(out.circles[i].x - out.circles[j].x,
                           out.circles[i].y - out.circles[j].y) /
                    (circles[i].radius + circles[j].radius);
        out.maximum_overlap = std::max(out.maximum_overlap, depth);
      }
  }
  out.feasible &=
      out.maximum_overlap < 1e-6 &&
      out.maximum_boundary_violation <= 1e-7 * std::max(xh - xl, yh - yl);
  out.seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - begin)
          .count();
  return out;
}
} // namespace DreamPlace::native
