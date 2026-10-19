// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/seed.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
namespace DreamPlace::native::macro {
namespace {
double rounded(double value, double pitch) {
  return std::ceil(value / pitch - 1e-9) * pitch;
}
struct Shape {
  int cols, rows;
  double width, height;
};
} // namespace
Topology seed_topology(const Design &d, std::span<const NodeInfo> info,
                       std::span<const Flow> incoming_flow, const Parameters &p,
                       const Topology &requested) {
  const auto input = remap_unstable(requested, d, info);
  d.validate();
  p.validate();
  if (info.size() != d.node_x.size())
    throw std::invalid_argument("seed metadata size mismatch");
  // Hierarchical HHDS ports can be buses, whereas the physical netlist has
  // one fixed terminal per bit. Expand the bundle to real bit terminals so
  // bus IOs participate in anchors and direct pseudo-nets as scalar IOs do.
  std::map<std::string, std::vector<std::string>> bundles;
  std::set<std::string> exact;
  for (const auto &n : info) {
    exact.insert(n.name);
    if (n.io && n.name.ends_with("]")) {
      auto at = n.name.rfind('[');
      if (at != std::string::npos)
        bundles[n.name.substr(0, at)].push_back(n.name);
    }
  }
  std::map<std::pair<std::string, std::string>, double> logical;
  for (const auto &f : incoming_flow)
    logical[{f.from, f.to}] += f.weight;
  for (const auto &f : input.dataflow)
    logical[{f.from, f.to}] = f.weight;
  auto endpoints = [&](const std::string &name) {
    if (!exact.contains(name) && bundles.contains(name))
      return bundles.at(name);
    return std::vector<std::string>{name};
  };
  std::vector<Flow> flow;
  for (const auto &[ends, weight] : logical)
    for (const auto &from : endpoints(ends.first))
      for (const auto &to : endpoints(ends.second))
        flow.push_back({from, to, weight});
  Topology t;
  t.design = input.design;
  t.dataflow.assign(flow.begin(), flow.end());
  t.diagnostics = input.diagnostics;
  std::map<std::string, int> names;
  for (int i = 0; i < int(info.size()); ++i)
    names[info[i].name] = i;
  std::set<std::string> constrained_fixed;
  for (const auto &f : input.fixed)
    constrained_fixed.insert(f.name);
  std::set<std::string> used;
  for (const auto &g : input.groups)
    for (const auto &m : g.members)
      if (names.contains(m.name))
        used.insert(m.name);
  std::map<std::string, std::vector<int>> buckets;
  for (int i = 0; i < int(info.size()); ++i)
    if (info[i].macro && !used.contains(info[i].name)) {
      const auto &n = info[i];
      std::string key = n.parent + "/" + n.master + "/" +
                        std::to_string(d.node_size_x[i]) + "x" +
                        std::to_string(d.node_size_y[i]);
      if (i >= d.num_movable || constrained_fixed.contains(n.name))
        key = "fixed/" + n.name;
      buckets[key].push_back(i);
      if (!n.stable)
        t.diagnostics.push_back("unstable macro key; signature/ordinal "
                                "fallback within parent: " +
                                n.name);
    }
  const auto &l = d.layout;
  std::map<std::string, std::vector<Shape>> curves;
  bool has_halo = false;
  for (size_t i = 0; i < info.size(); ++i)
    if (info[i].macro)
      has_halo |= (!d.node_halo_x.empty() && d.node_halo_x[i] > 0) ||
                  (!d.node_halo_y.empty() && d.node_halo_y[i] > 0);
  const double channel_x = p.min_channel > 0
                               ? p.min_channel
                               : (has_halo ? 0 : 10 * l.site_width),
               channel_y = p.min_channel > 0
                               ? p.min_channel
                               : (has_halo ? 0 : 2 * l.row_height);
  std::function<void(std::string, std::vector<int>)> add =
      [&](std::string key, std::vector<int> nodes) {
        std::sort(nodes.begin(), nodes.end(),
                  [&](int a, int b) { return info[a].name < info[b].name; });
        int first = nodes.front(), count = nodes.size();
        double hx = d.node_halo_x.empty()
                        ? 0
                        : rounded(d.node_halo_x[first], l.site_width),
               hy = d.node_halo_y.empty()
                        ? 0
                        : rounded(d.node_halo_y[first], l.row_height),
               w = d.node_size_x[first], h = d.node_size_y[first];
        for (int i : nodes) {
          if (!d.node_halo_x.empty())
            hx = std::max(hx, rounded(d.node_halo_x[i], l.site_width));
          if (!d.node_halo_y.empty())
            hy = std::max(hy, rounded(d.node_halo_y[i], l.row_height));
        }
        hx = std::max(hx, rounded(channel_x / 2, l.site_width));
        hy = std::max(hy, rounded(channel_y / 2, l.row_height));
        double px = rounded(w + 2 * hx, l.site_width),
               py = rounded(h + 2 * hy, l.row_height);
        std::vector<Shape> shapes;
        for (int cols = 1; cols <= count; ++cols)
          if (count % cols == 0) {
            int rows = count / cols;
            double width = (cols - 1) * px + w, height = (rows - 1) * py + h;
            if (width + 2 * hx <= l.xh - l.xl && height + 2 * hy <= l.yh - l.yl)
              shapes.push_back({cols, rows, width, height});
          }
        if (shapes.empty()) {
          if (count == 1)
            throw std::runtime_error("H2D macro with halo does not fit core: " +
                                     info[first].name);
          int half = count / 2;
          add(key + "/0",
              std::vector<int>(nodes.begin(), nodes.begin() + half));
          add(key + "/1", std::vector<int>(nodes.begin() + half, nodes.end()));
          return;
        }
        const double aspect = (l.xh - l.xl) / (l.yh - l.yl);
        auto best = std::min_element(
            shapes.begin(), shapes.end(), [&](const auto &a, const auto &b) {
              return std::abs(std::log(a.width / a.height / aspect)) <
                     std::abs(std::log(b.width / b.height / aspect));
            });
        Group g;
        g.id = "hier:" + key;
        g.cols = best->cols;
        g.rows = best->rows;
        g.width = best->width;
        g.height = best->height;
        g.halo_x = hx;
        g.halo_y = hy;
        g.fixed = first >= d.num_movable ||
                  constrained_fixed.contains(info[first].name);
        double cx = 0, cy = 0;
        for (int i : nodes) {
          cx += d.node_x[i] + w / 2;
          cy += d.node_y[i] + h / 2;
        }
        cx /= count;
        cy /= count;
        g.x = l.xl + std::round((cx - g.width / 2 - l.xl) / l.site_width) *
                         l.site_width;
        g.y = l.yl + std::round((cy - g.height / 2 - l.yl) / l.row_height) *
                         l.row_height;
        g.x = std::clamp(g.x, l.xl + hx, l.xh - g.width - hx);
        g.y = std::clamp(g.y, l.yl + hy, l.yh - g.height - hy);
        for (int k = 0; k < count; ++k) {
          int i = nodes[k];
          g.members.push_back(
              {info[i].name, (k % g.cols) * px, (k / g.cols) * py,
               d.node_orientation.empty() ? Orientation::N
                                          : d.node_orientation[i]});
        }
        if (g.fixed) {
          g.x = d.node_x[first];
          g.y = d.node_y[first];
          g.allowed_orientations = 1;
        } else
          g.allowed_orientations = group_orientations(g, info, names);
        if (!(g.allowed_orientations & 1))
          throw std::invalid_argument(
              "seed input orientation violates LEF symmetry");
        curves[g.id] = shapes;
        t.groups.push_back(std::move(g));
      };
  for (auto &[key, nodes] : buckets)
    add(key, std::move(nodes));
  t = merge(t, input, d, info);
  t.dataflow.assign(flow.begin(), flow.end());
  for (auto &g : t.groups)
    for (const auto &fixed : t.fixed)
      if (std::any_of(g.members.begin(), g.members.end(),
                      [&](const Member &m) { return m.name == fixed.name; })) {
        if (g.members.size() != 1)
          throw std::invalid_argument("fixed macro must be a singleton group");
        g.x = fixed.x;
        g.y = fixed.y;
        g.orientation = Orientation::N;
        g.members[0].orientation = fixed.orientation;
        g.fixed = true;
        g.allowed_orientations = 1;
      }

  // Native shape/position annealing over hierarchy clusters. Rigid tilings are
  // generated independently; no RTLMP/Re2MaP implementation is copied.
  std::map<std::string, int> member_group;
  for (int g = 0; g < int(t.groups.size()); ++g)
    for (const auto &m : t.groups[g].members)
      member_group[m.name] = g;
  struct Link {
    int a, b;
    double weight;
  };
  std::map<std::pair<int, int>, double> link_weights;
  for (const auto &f : flow)
    if (member_group.contains(f.from) && member_group.contains(f.to)) {
      int a = member_group[f.from], b = member_group[f.to];
      if (a != b)
        link_weights[{std::min(a, b), std::max(a, b)}] += f.weight;
    }
  std::vector<Link> links;
  for (const auto &[ends, w] : link_weights)
    links.push_back({ends.first, ends.second, w});
  auto score = [&](const Topology &v) {
    double cost = 0;
    for (const auto &e : links) {
      const auto &a = v.groups[e.a];
      const auto &b = v.groups[e.b];
      cost += e.weight * (std::abs(a.x + a.width / 2 - b.x - b.width / 2) +
                          std::abs(a.y + a.height / 2 - b.y - b.height / 2));
    }
    for (size_t i = 0; i < v.groups.size(); ++i)
      for (size_t j = i + 1; j < v.groups.size(); ++j) {
        const auto &a = v.groups[i];
        const auto &b = v.groups[j];
        double x = std::max(0., std::min(a.x + a.width + a.halo_x,
                                         b.x + b.width + b.halo_x) -
                                    std::max(a.x - a.halo_x, b.x - b.halo_x)),
               y = std::max(0., std::min(a.y + a.height + a.halo_y,
                                         b.y + b.height + b.halo_y) -
                                    std::max(a.y - a.halo_y, b.y - b.halo_y));
        cost +=
            p.M * x * y /
            std::max(1., (a.width + a.halo_x * 2) * (a.height + a.halo_y * 2));
      }
    return cost;
  };
  if (!t.groups.empty() && p.seed_sa_evals) {
    std::mt19937 rng(p.seed);
    std::uniform_real_distribution<double> unit(0, 1);
    Topology current = t, best = t;
    double value = score(t), bestvalue = value;
    for (int eval = 0; eval < p.seed_sa_evals; ++eval) {
      auto candidate = current;
      int gi = rng() % candidate.groups.size();
      auto &g = candidate.groups[gi];
      if (g.fixed || g.freeze_topology)
        continue;
      auto found = curves.find(g.id);
      if (found != curves.end() && unit(rng) < .5) {
        const auto &shape = found->second[rng() % found->second.size()];
        double cx = g.x + g.width / 2, cy = g.y + g.height / 2;
        g.cols = shape.cols;
        g.rows = shape.rows;
        g.width = shape.width;
        g.height = shape.height;
        int first = names.at(g.members.front().name);
        double px = rounded(d.node_size_x[first] + 2 * g.halo_x, l.site_width),
               py = rounded(d.node_size_y[first] + 2 * g.halo_y, l.row_height);
        for (int k = 0; k < int(g.members.size()); ++k) {
          g.members[k].dx = (k % g.cols) * px;
          g.members[k].dy = (k / g.cols) * py;
        }
        g.x = cx - g.width / 2;
        g.y = cy - g.height / 2;
      } else {
        g.x += (unit(rng) - .5) * (l.xh - l.xl) * .2;
        g.y += (unit(rng) - .5) * (l.yh - l.yl) * .2;
      }
      g.x = l.xl + std::round((g.x - l.xl) / l.site_width) * l.site_width;
      g.y = l.yl + std::round((g.y - l.yl) / l.row_height) * l.row_height;
      g.x = std::clamp(g.x, l.xl + g.halo_x, l.xh - g.width - g.halo_x);
      g.y = std::clamp(g.y, l.yl + g.halo_y, l.yh - g.height - g.halo_y);
      double trial = score(candidate),
             temperature = std::max(1e-12, p.T0 * std::pow(p.cooling, eval));
      if (trial <= value ||
          unit(rng) < std::exp(std::min(0., (value - trial) / temperature))) {
        current = std::move(candidate);
        value = trial;
        if (value < bestvalue) {
          best = current;
          bestvalue = value;
        }
      }
    }
    t = std::move(best);
  }
  member_group.clear();
  for (int g = 0; g < int(t.groups.size()); ++g)
    for (const auto &m : t.groups[g].members)
      member_group[m.name] = g;
  if (p.anchors) {
    double delta =
        p.delta > 0 ? p.delta : 2 * std::max(l.site_width, l.row_height);
    std::vector<std::array<double, 4>> io_by_group(t.groups.size());
    auto add_io = [&](const std::string &member, const std::string &other,
                      double weight) {
      auto group = member_group.find(member);
      if (group == member_group.end())
        return;
      auto endpoint = names.find(other);
      if (endpoint == names.end() || !info[endpoint->second].io)
        return;
      int i = endpoint->second;
      std::array<double, 4> edge{
          std::abs(d.node_x[i] - l.xl), std::abs(l.xh - d.node_x[i]),
          std::abs(d.node_y[i] - l.yl), std::abs(l.yh - d.node_y[i])};
      io_by_group[group->second][std::min_element(edge.begin(), edge.end()) -
                                 edge.begin()] += weight;
    };
    for (const auto &f : flow) {
      add_io(f.from, f.to, f.weight);
      add_io(f.to, f.from, f.weight);
    }
    for (const auto &g : t.groups)
      if (!g.fixed) {
        std::array<double, 4> distance{
            g.x - g.halo_x - l.xl, l.xh - g.x - g.width - g.halo_x,
            g.y - g.halo_y - l.yl, l.yh - g.y - g.height - g.halo_y};
        int side = std::min_element(distance.begin(), distance.end()) -
                   distance.begin();
        const auto &io = io_by_group[&g - t.groups.data()];
        int ioside = std::max_element(io.begin(), io.end()) - io.begin();
        if (io[ioside] > 0)
          side = ioside;
        if (distance[side] <= delta || io[side] > 0) {
          std::string id = "anchor:" + g.id;
          if (std::none_of(t.anchors.begin(), t.anchors.end(),
                           [&](const auto &a) { return a.id == id; }))
            t.anchors.push_back(
                {id, g.id, static_cast<Side>(side), delta, {false, 1}});
        }
      }
  }
  if (p.order && !links.empty()) {
    std::vector<double> w;
    for (const auto &e : links)
      w.push_back(e.weight);
    std::sort(w.begin(), w.end());
    double threshold =
        w[std::min(w.size() - 1, size_t(p.tau_q * (w.size() - 1)))];
    for (const auto &e : links)
      if (e.weight >= threshold) {
        const auto &a = t.groups[e.a];
        const auto &b = t.groups[e.b];
        double dx = b.x + b.width / 2 - a.x - a.width / 2,
               dy = b.y + b.height / 2 - a.y - a.height / 2;
        bool x = std::abs(dx) / (l.xh - l.xl) >= std::abs(dy) / (l.yh - l.yl);
        bool forward = x ? dx > 0 : dy > 0;
        if ((x ? dx : dy) == 0)
          forward = a.id < b.id;
        std::string from = forward ? a.id : b.id, to = forward ? b.id : a.id,
                    id = "order:" + from + ":" + to;
        if (std::none_of(t.order.begin(), t.order.end(),
                         [&](const auto &r) { return r.id == id; }))
          t.order.push_back(
              {id, from, to, x ? Axis::x : Axis::y, {false, e.weight}});
      }
    // Transitive reduction acts on the relation; surviving edges retain their
    // original maximum dataflow weight rather than summing path weights.
    for (size_t skip = 0; skip < t.order.size();) {
      auto edge = t.order[skip];
      if (edge.strength.hard || !edge.id.starts_with("order:")) {
        ++skip;
        continue;
      }
      std::set<std::string> seen{edge.from};
      std::vector<std::string> todo{edge.from};
      while (!todo.empty()) {
        auto v = todo.back();
        todo.pop_back();
        for (size_t j = 0; j < t.order.size(); ++j)
          if (j != skip && t.order[j].axis == edge.axis &&
              t.order[j].from == v && seen.insert(t.order[j].to).second)
            todo.push_back(t.order[j].to);
      }
      if (seen.contains(edge.to))
        t.order.erase(t.order.begin() + skip);
      else
        ++skip;
    }
  }
  t.macro_fingerprint = fingerprint(d, info);
  t.netlist_hash = fingerprint(d, info, true);
  validate(t, d, info);
  return t;
}
} // namespace DreamPlace::native::macro
