// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/render.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
namespace DreamPlace::native::macro {
namespace {
using Bytes = std::vector<unsigned char>;
using Color = std::array<unsigned char, 3>;
void word(Bytes &v, std::uint32_t n) {
  for (int s : {24, 16, 8, 0})
    v.push_back(n >> s);
}
std::uint32_t crc(std::span<const unsigned char> s) {
  std::uint32_t value = ~0U;
  for (unsigned char c : s) {
    value ^= c;
    for (int k = 0; k < 8; ++k)
      value = (value >> 1) ^ (0xedb88320U & -(value & 1));
  }
  return ~value;
}
void chunk(Bytes &out, const char *type, const Bytes &data) {
  word(out, data.size());
  size_t start = out.size();
  for (int k = 0; k < 4; ++k)
    out.push_back(type[k]);
  out.insert(out.end(), data.begin(), data.end());
  word(out, crc(std::span(out).subspan(start)));
}
} // namespace
void write_png(const std::filesystem::path &path, const Design &d,
               std::span<const NodeInfo> info, const Topology &t, int width) {
  validate(t, d, info);
  if (width < 64 || width > 4096)
    throw std::invalid_argument("invalid layout PNG width");
  const auto &l = d.layout;
  int height = std::clamp(
      int(std::round(width * std::min(64., (l.yh - l.yl) / (l.xh - l.xl)))), 64,
      4096);
  Bytes pixels(width * height * 3, 15);
  auto px = [&](double x) {
    return int(std::round((x - l.xl) / (l.xh - l.xl) * (width - 1)));
  };
  auto py = [&](double y) {
    return int(std::round((l.yh - y) / (l.yh - l.yl) * (height - 1)));
  };
  auto dot = [&](int x, int y, Color color) {
    if (x >= 0 && x < width && y >= 0 && y < height)
      for (int c = 0; c < 3; ++c)
        pixels[(y * width + x) * 3 + c] = color[c];
  };
  auto line = [&](int x0, int y0, int x1, int y1, Color c) {
    int steps = std::max(std::abs(x1 - x0), std::abs(y1 - y0));
    if (!steps) {
      dot(x0, y0, c);
      return;
    }
    for (int k = 0; k <= steps; ++k)
      dot(x0 + (x1 - x0) * double(k) / steps,
          y0 + (y1 - y0) * double(k) / steps, c);
  };
  auto rect = [&](double x, double y, double w, double h, Color c, bool fill) {
    int x0 = std::clamp(px(x), 0, width - 1),
        x1 = std::clamp(px(x + w), 0, width - 1),
        y0 = std::clamp(py(y + h), 0, height - 1),
        y1 = std::clamp(py(y), 0, height - 1);
    if (fill)
      for (int yy = y0; yy <= y1; ++yy)
        for (int xx = x0; xx <= x1; ++xx)
          dot(xx, yy, c);
    else {
      line(x0, y0, x1, y0, c);
      line(x0, y1, x1, y1, c);
      line(x0, y0, x0, y1, c);
      line(x1, y0, x1, y1, c);
    }
  };
  for (int i = 0; i < int(info.size()); ++i) {
    if (info[i].io) {
      int x = px(d.node_x[i]), y = py(d.node_y[i]);
      for (int dx = -2; dx <= 2; ++dx)
        for (int dy = -2; dy <= 2; ++dy)
          dot(x + dx, y + dy, {255, 230, 80});
    } else if (d.node_size_x[i] > 0 && d.node_size_y[i] > 0)
      rect(d.node_x[i], d.node_y[i], d.node_size_x[i], d.node_size_y[i],
           info[i].macro ? Color{132, 50, 174} : Color{42, 204, 215}, true);
  }
  std::map<std::string, const Group *> groups;
  for (const auto &g : t.groups) {
    groups[g.id] = &g;
    rect(g.x - g.halo_x, g.y - g.halo_y, g.width + 2 * g.halo_x,
         g.height + 2 * g.halo_y, {90, 90, 110}, false);
    rect(g.x, g.y, g.width, g.height, {220, 180, 240}, false);
  }
  for (const auto &a : t.anchors) {
    const auto &g = *groups.at(a.group);
    double x = g.x + g.width / 2, y = g.y + g.height / 2, ex = x, ey = y;
    switch (a.side) {
    case Side::left:
      ex = l.xl;
      break;
    case Side::right:
      ex = l.xh;
      break;
    case Side::bottom:
      ey = l.yl;
      break;
    case Side::top:
      ey = l.yh;
      break;
    }
    line(px(x), py(y), px(ex), py(ey), {255, 220, 60});
  }
  for (const auto &e : t.order) {
    const auto &a = *groups.at(e.from);
    const auto &b = *groups.at(e.to);
    bool satisfied = e.axis == Axis::x ? a.x + a.width <= b.x + 1e-6
                                       : a.y + a.height <= b.y + 1e-6;
    Color color = satisfied ? Color{90, 240, 120} : Color{255, 65, 65};
    int x = px(a.x + a.width / 2), y = py(a.y + a.height / 2),
        ex = px(b.x + b.width / 2), ey = py(b.y + b.height / 2);
    line(x, y, ex, ey, color);
    double angle = std::atan2(ey - y, ex - x);
    for (double delta : {-.5, .5})
      line(ex, ey, ex - int(8 * std::cos(angle + delta)),
           ey - int(8 * std::sin(angle + delta)), color);
  }
  Bytes raw;
  raw.reserve((width * 3 + 1) * height);
  for (int y = 0; y < height; ++y) {
    raw.push_back(0);
    raw.insert(raw.end(), pixels.begin() + y * width * 3,
               pixels.begin() + (y + 1) * width * 3);
  }
  Bytes z{0x78, 0x01};
  for (size_t start = 0; start < raw.size();) {
    size_t n = std::min<size_t>(65535, raw.size() - start);
    z.push_back(start + n == raw.size() ? 1 : 0);
    z.push_back(n & 255);
    z.push_back(n >> 8);
    unsigned inv = unsigned(n) ^ 65535;
    z.push_back(inv & 255);
    z.push_back(inv >> 8);
    z.insert(z.end(), raw.begin() + start, raw.begin() + start + n);
    start += n;
  }
  std::uint32_t a = 1, b = 0;
  for (unsigned char c : raw) {
    a = (a + c) % 65521;
    b = (b + a) % 65521;
  }
  word(z, (b << 16) | a);
  Bytes png{137, 80, 78, 71, 13, 10, 26, 10}, header;
  word(header, width);
  word(header, height);
  header.insert(header.end(), {8, 2, 0, 0, 0});
  chunk(png, "IHDR", header);
  chunk(png, "IDAT", z);
  chunk(png, "IEND", {});
  auto temporary = path;
  temporary += ".tmp";
  std::ofstream out(temporary, std::ios::binary);
  out.write(reinterpret_cast<const char *>(png.data()), png.size());
  out.close();
  if (!out)
    throw std::runtime_error("cannot write H2D PNG");
  std::filesystem::rename(temporary, path);
}
} // namespace DreamPlace::native::macro
