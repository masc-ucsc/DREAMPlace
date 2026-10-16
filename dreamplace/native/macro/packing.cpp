// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/packing.hpp"
#include "dreamplace/native/macro/check.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <optional>
#include <tuple>
namespace DreamPlace::native::macro {
namespace {
using Int = std::int64_t;
struct Rect {
  Int x, y, w, h;
};
bool intersects(Rect a, Rect b) {
  return a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h &&
         b.y < a.y + a.h;
}
bool contains(Rect a, Rect b) {
  return a.x <= b.x && a.y <= b.y && a.x + a.w >= b.x + b.w &&
         a.y + a.h >= b.y + b.h;
}
void occupy(std::vector<Rect> &free, Rect used) {
  std::vector<Rect> next;
  for (auto f : free) {
    if (!intersects(f, used)) {
      next.push_back(f);
      continue;
    }
    if (used.x > f.x)
      next.push_back({f.x, f.y, used.x - f.x, f.h});
    if (used.x + used.w < f.x + f.w)
      next.push_back({used.x + used.w, f.y, f.x + f.w - used.x - used.w, f.h});
    if (used.y > f.y)
      next.push_back({f.x, f.y, f.w, used.y - f.y});
    if (used.y + used.h < f.y + f.h)
      next.push_back({f.x, used.y + used.h, f.w, f.y + f.h - used.y - used.h});
  }
  std::vector<Rect> maximal;
  for (size_t i = 0; i < next.size(); ++i) {
    bool redundant = false;
    for (size_t j = 0; j < next.size(); ++j)
      if (i != j && contains(next[j], next[i]) &&
          (!contains(next[i], next[j]) || j < i)) {
        redundant = true;
        break;
      }
    if (!redundant && next[i].w > 0 && next[i].h > 0)
      maximal.push_back(next[i]);
  }
  free = std::move(maximal);
}
} // namespace
bool pack_rectangles(Design &source) {
  source.validate();
  const auto &l = source.layout;
  auto up = [](double x) { return Int(std::ceil(x - 1e-9)); };
  auto down = [](double x) { return Int(std::floor(x + 1e-9)); };
  if ((l.xh - l.xl) / l.site_width > 1e9 || (l.yh - l.yl) / l.row_height > 1e9)
    return false;
  Rect outline{0, 0, down((l.xh - l.xl) / l.site_width),
               down((l.yh - l.yl) / l.row_height)};
  std::vector<Rect> boxes;
  std::vector<Int> hx, hy;
  std::vector<int> movable, audit;
  std::vector<Rect> initial{outline};
  for (int i = 0; i < int(source.node_x.size()); ++i) {
    for (double value :
         {source.node_size_x[i] / l.site_width,
          source.node_size_y[i] / l.row_height,
          source.node_halo_x.empty() ? 0 : source.node_halo_x[i] / l.site_width,
          source.node_halo_y.empty() ? 0
                                     : source.node_halo_y[i] / l.row_height})
      if (!std::isfinite(value) || value > 1e9)
        return false;
    if (i >= source.num_movable && source.node_size_x[i] > 0 &&
        source.node_size_y[i] > 0 &&
        (std::abs((source.node_x[i] - l.xl) / l.site_width) > 1e9 ||
         std::abs((source.node_y[i] - l.yl) / l.row_height) > 1e9))
      return false;
    Int x = source.node_halo_x.empty()
                ? 0
                : up(source.node_halo_x[i] / l.site_width),
        y = source.node_halo_y.empty()
                ? 0
                : up(source.node_halo_y[i] / l.row_height);
    hx.push_back(x);
    hy.push_back(y);
    Rect box{0, 0, up(source.node_size_x[i] / l.site_width) + 2 * x,
             up(source.node_size_y[i] / l.row_height) + 2 * y};
    boxes.push_back(box);
    if (source.node_size_x[i] <= 0 || source.node_size_y[i] <= 0)
      continue;
    audit.push_back(i);
    if (i < source.num_movable)
      movable.push_back(i);
    else {
      box.x = down((source.node_x[i] - l.xl) / l.site_width) - x;
      box.y = down((source.node_y[i] - l.yl) / l.row_height) - y;
      box.w =
          up((source.node_x[i] + source.node_size_x[i] - l.xl) / l.site_width) +
          x - box.x;
      box.h =
          up((source.node_y[i] + source.node_size_y[i] - l.yl) / l.row_height) +
          y - box.y;
      occupy(initial, box);
    }
  }
  for (int order = 0; order < 4; ++order) {
    auto ids = movable;
    std::stable_sort(ids.begin(), ids.end(), [&](int a, int b) {
      auto u = boxes[a], v = boxes[b];
      auto key = [&](Rect r) {
        switch (order) {
        case 0:
          return std::pair{r.w * r.h, std::max(r.w, r.h)};
        case 1:
          return std::pair{std::max(r.w, r.h), r.w * r.h};
        case 2:
          return std::pair{r.h, r.w};
        default:
          return std::pair{r.w, r.h};
        }
      };
      return key(u) > key(v);
    });
    auto free = initial;
    auto placed = source;
    bool success = true;
    for (int i : ids) {
      auto box = boxes[i];
      using Score = std::tuple<Int, Int, double, Int, Int>;
      std::optional<Score> best;
      Rect selected{};
      for (auto f : free)
        if (box.w <= f.w && box.h <= f.h) {
          double distance =
              std::abs(l.xl + (f.x + hx[i]) * l.site_width - source.node_x[i]) +
              std::abs(l.yl + (f.y + hy[i]) * l.row_height - source.node_y[i]);
          Score score{std::min(f.w - box.w, f.h - box.h),
                      f.w * f.h - box.w * box.h, distance, f.y, f.x};
          if (!best || score < *best) {
            best = score;
            selected = {f.x, f.y, box.w, box.h};
          }
        }
      if (!best) {
        success = false;
        break;
      }
      placed.node_x[i] = l.xl + (selected.x + hx[i]) * l.site_width;
      placed.node_y[i] = l.yl + (selected.y + hy[i]) * l.row_height;
      occupy(free, selected);
    }
    if (success && check(placed, audit).legal()) {
      source = std::move(placed);
      return true;
    }
  }
  return false;
}
} // namespace DreamPlace::native::macro
