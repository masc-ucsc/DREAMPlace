// BSD-3-Clause. See LICENSE. Small dependency-free RGB debug renderer.
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <span>
#include <sstream>
#include <vector>

#include "dreamplace/native/synthetic.hpp"

namespace DreamPlace::native::synthetic {
namespace {
using Bytes = std::vector<std::uint8_t>;
struct Color {
  int r, g, b;
};
constexpr Color ink{30, 45, 65};
struct Image {
  int width, height;
  Bytes rgb;
  Image(int w, int h) : width(w), height(h), rgb(3 * w * h, 248) {}
  void pixel(int x, int y, Color c, int alpha = 255) {
    if (x < 0 || x >= width || y < 0 || y >= height) return;
    size_t p = 3 * (y * width + x);
    for (auto v : {c.r, c.g, c.b}) {
      rgb[p] = (v * alpha + rgb[p] * (255 - alpha) + 127) / 255;
      ++p;
    }
  }
  void rect(int x, int y, int w, int h, Color c, int alpha = 255) {
    for (int j = std::max(y, 0); j < std::min(y + h, height); ++j)
      for (int i = std::max(x, 0); i < std::min(x + w, width); ++i) pixel(i, j, c, alpha);
  }
  void line(int x, int y, int xx, int yy, Color c, int alpha = 255) {
    int dx = std::abs(xx - x), sx = x < xx ? 1 : -1, dy = -std::abs(yy - y), sy = y < yy ? 1 : -1, err = dx + dy;
    for (;;) {
      pixel(x, y, c, alpha);
      if (x == xx && y == yy) break;
      int twice = 2 * err;
      if (twice >= dy) {
        err += dy;
        x += sx;
      }
      if (twice <= dx) {
        err += dx;
        y += sy;
      }
    }
  }
};
// Hand-specified 5x7 ASCII column glyphs for readable, portable debug labels.
std::array<unsigned char, 5> glyph(char c) {
  switch (std::toupper(static_cast<unsigned char>(c))) {
    case 'A':
      return {0x7e, 0x11, 0x11, 0x11, 0x7e};
    case 'B':
      return {0x7f, 0x49, 0x49, 0x49, 0x36};
    case 'C':
      return {0x3e, 0x41, 0x41, 0x41, 0x22};
    case 'D':
      return {0x7f, 0x41, 0x41, 0x22, 0x1c};
    case 'E':
      return {0x7f, 0x49, 0x49, 0x49, 0x41};
    case 'F':
      return {0x7f, 0x09, 0x09, 0x09, 0x01};
    case 'G':
      return {0x3e, 0x41, 0x49, 0x49, 0x7a};
    case 'H':
      return {0x7f, 0x08, 0x08, 0x08, 0x7f};
    case 'I':
      return {0x41, 0x41, 0x7f, 0x41, 0x41};
    case 'J':
      return {0x20, 0x40, 0x41, 0x3f, 0x01};
    case 'K':
      return {0x7f, 0x08, 0x14, 0x22, 0x41};
    case 'L':
      return {0x7f, 0x40, 0x40, 0x40, 0x40};
    case 'M':
      return {0x7f, 0x02, 0x0c, 0x02, 0x7f};
    case 'N':
      return {0x7f, 0x04, 0x08, 0x10, 0x7f};
    case 'O':
      return {0x3e, 0x41, 0x41, 0x41, 0x3e};
    case 'P':
      return {0x7f, 0x09, 0x09, 0x09, 0x06};
    case 'Q':
      return {0x3e, 0x41, 0x51, 0x21, 0x5e};
    case 'R':
      return {0x7f, 0x09, 0x19, 0x29, 0x46};
    case 'S':
      return {0x46, 0x49, 0x49, 0x49, 0x31};
    case 'T':
      return {0x01, 0x01, 0x7f, 0x01, 0x01};
    case 'U':
      return {0x3f, 0x40, 0x40, 0x40, 0x3f};
    case 'V':
      return {0x1f, 0x20, 0x40, 0x20, 0x1f};
    case 'W':
      return {0x7f, 0x20, 0x18, 0x20, 0x7f};
    case 'X':
      return {0x63, 0x14, 0x08, 0x14, 0x63};
    case 'Y':
      return {0x07, 0x08, 0x70, 0x08, 0x07};
    case 'Z':
      return {0x61, 0x51, 0x49, 0x45, 0x43};
    case '0':
      return {0x3e, 0x51, 0x49, 0x45, 0x3e};
    case '1':
      return {0x00, 0x42, 0x7f, 0x40, 0x00};
    case '2':
      return {0x62, 0x51, 0x49, 0x49, 0x46};
    case '3':
      return {0x22, 0x41, 0x49, 0x49, 0x36};
    case '4':
      return {0x18, 0x14, 0x12, 0x7f, 0x10};
    case '5':
      return {0x27, 0x45, 0x45, 0x45, 0x39};
    case '6':
      return {0x3c, 0x4a, 0x49, 0x49, 0x30};
    case '7':
      return {0x01, 0x71, 0x09, 0x05, 0x03};
    case '8':
      return {0x36, 0x49, 0x49, 0x49, 0x36};
    case '9':
      return {0x06, 0x49, 0x49, 0x29, 0x1e};
    case '.':
      return {0, 0x60, 0x60, 0, 0};
    case ':':
      return {0, 0x36, 0x36, 0, 0};
    case '-':
      return {0x08, 0x08, 0x08, 0x08, 0x08};
    case '_':
      return {0x40, 0x40, 0x40, 0x40, 0x40};
    case '+':
      return {0x08, 0x08, 0x3e, 0x08, 0x08};
    case '/':
      return {0x20, 0x10, 0x08, 0x04, 0x02};
    case '%':
      return {0x63, 0x13, 0x08, 0x64, 0x63};
    case '|':
      return {0, 0, 0x7f, 0, 0};
    default:
      return {};
  }
}
void text(Image& image, int x, int y, const std::string& str, Color color = ink, int scale = 2) {
  for (char c : str) {
    auto g = glyph(c);
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 7; ++j)
        if ((g[i] >> j) & 1) image.rect(x + i * scale, y + j * scale, scale, scale, color);
    x += 6 * scale;
  }
}
std::string number(double v) {
  std::ostringstream s;
  s << std::fixed << std::setprecision(1) << v;
  return s.str();
}
Color cell_color(int group, int max_group) {
  constexpr std::array<Color, 6> palette{
      {{52, 103, 190}, {47, 150, 190}, {47, 168, 129}, {166, 182, 62}, {235, 163, 53}, {208, 88, 77}}};
  double t = 5.0 * group / std::max(1, max_group);
  int i = std::min(4, static_cast<int>(t));
  double u = t - i;
  auto a = palette[i], b = palette[i + 1];
  return {static_cast<int>(a.r * (1 - u) + b.r * u), static_cast<int>(a.g * (1 - u) + b.g * u),
          static_cast<int>(a.b * (1 - u) + b.b * u)};
}
void panel(Image& image, const HhdsDesign& e, int left) {
  const auto& d = e.design;
  const auto& l = d.layout;
  // Show the die margin when there are IO pins, so perimeter anchors and their
  // connections are visible rather than clipped at the core boundary.
  bool io = !e.terminals.empty();
  double xl = io ? 0 : l.xl, yl = io ? 0 : l.yl;
  double xh = io ? l.xh + 2 : l.xh, yh = io ? l.yh + 2 : l.yh;
  double scale = std::min(720 / (xh - xl), 640 / (yh - yl));
  int w = std::llround((xh - xl) * scale), h = std::llround((yh - yl) * scale);
  int x0 = left + (736 - w) / 2, y0 = 150 + (640 - h) / 2;
  auto xy = [&](double x, double y) {
    return std::pair<int, int>{x0 + std::llround((x - xl) * scale), y0 + h - std::llround((y - yl) * scale)};
  };
  image.rect(x0, y0, w + 1, h + 1, {227, 234, 243});
  auto [cx, cy] = xy(l.xl, l.yh);
  auto [cxh, cyl] = xy(l.xh, l.yl);
  image.rect(cx, cy, cxh - cx + 1, cyl - cy + 1, {255, 255, 255});
  for (double y = l.yl; y <= l.yh; y += l.row_height) {
    auto [x1, y1] = xy(l.xl, y);
    auto [x2, y2] = xy(l.xh, y);
    image.line(x1, y1, x2, y2, {233, 237, 242});
  }
  for (size_t n = 0; n < d.net_weights.size(); ++n) {
    int p = d.flat_net2pin[d.flat_net2pin_start[n]], node = d.pin2node[p];
    auto a = xy(d.node_x[node] + d.pin_offset_x[p], d.node_y[node] + d.pin_offset_y[p]);
    for (int j = d.flat_net2pin_start[n] + 1; j < d.flat_net2pin_start[n + 1]; ++j) {
      int q = d.flat_net2pin[j], other = d.pin2node[q];
      auto b = xy(d.node_x[other] + d.pin_offset_x[q], d.node_y[other] + d.pin_offset_y[q]);
      image.line(a.first, a.second, b.first, b.second, {50, 68, 92}, 42);
    }
  }
  int max_group = 0;
  for (auto n : e.cells) max_group = std::max(max_group, n.attr(cell_state).get().group);
  for (size_t i = 0; i < e.cells.size(); ++i) {
    auto [x, y] = xy(d.node_x[i], d.node_y[i] + d.node_size_y[i]);
    int cw = std::max(2, static_cast<int>(std::llround(d.node_size_x[i] * scale)));
    int ch = std::max(2, static_cast<int>(std::llround(d.node_size_y[i] * scale)));
    image.rect(x, y, cw, ch, ink);
    image.rect(x + 1, y + 1, std::max(1, cw - 2), std::max(1, ch - 2),
               cell_color(e.cells[i].attr(cell_state).get().group, max_group));
  }
  image.line(cx, cy, cxh, cy, ink);
  image.line(cx, cyl, cxh, cyl, ink);
  image.line(cx, cy, cx, cyl, ink);
  image.line(cxh, cy, cxh, cyl, ink);
  image.line(x0, y0, x0 + w, y0, ink);
  image.line(x0, y0 + h, x0 + w, y0 + h, ink);
  image.line(x0, y0, x0, y0 + h, ink);
  image.line(x0 + w, y0, x0 + w, y0 + h, ink);
  for (size_t i = 0; i < e.terminals.size(); ++i) {
    auto [x, y] = xy(d.node_x[e.cells.size() + i], d.node_y[e.cells.size() + i]);
    bool input = e.terminals[i].get_master_node().is_input_node();
    Color color = input ? Color{36, 85, 176} : Color{182, 55, 42};
    image.rect(x - 7, y - 7, 15, 15, ink);
    image.rect(x - 5, y - 5, 11, 11, color);
    bool many = e.terminals.size() > 2;
    std::string label = many ? physical_pin_name(e.terminals[i]) : (input ? "IN" : "OUT");
    int size = many ? 1 : 2;
    text(image, input ? x + 10 : x - 10 - 6 * size * static_cast<int>(label.size()), y - (many ? 13 : 25), label, color,
         size);
  }
  text(image, left + 8, 812, "CORE " + number(l.xh - l.xl) + " X " + number(l.yh - l.yl) + " UM");
}
void big_endian(Bytes& b, std::uint32_t n) {
  for (int shift : {24, 16, 8, 0}) b.push_back((n >> shift) & 255);
}
std::uint32_t crc(std::span<const std::uint8_t> data) {
  std::uint32_t c = 0xffffffff;
  for (auto v : data) {
    c ^= v;
    for (int bit = 0; bit < 8; ++bit) c = (c >> 1) ^ ((c & 1) ? 0xedb88320u : 0);
  }
  return ~c;
}
void chunk(Bytes& png, const std::string& type, std::span<const std::uint8_t> data) {
  big_endian(png, data.size());
  size_t start = png.size();
  png.insert(png.end(), type.begin(), type.end());
  png.insert(png.end(), data.begin(), data.end());
  big_endian(png, crc(std::span<const std::uint8_t>(png).subspan(start)));
}
void encode(const Image& image, const std::filesystem::path& file) {
  Bytes png{137, 80, 78, 71, 13, 10, 26, 10}, header;
  big_endian(header, image.width);
  big_endian(header, image.height);
  for (auto b : {8, 2, 0, 0, 0}) header.push_back(b);
  chunk(png, "IHDR", header);
  // PNG filter 0; zlib with stored DEFLATE blocks. This costs a few MB per
  // debug image, but needs no image, compression, Python or graphics library.
  Bytes raw;
  raw.reserve((3 * image.width + 1) * image.height);
  for (int y = 0; y < image.height; ++y) {
    raw.push_back(0);
    raw.insert(raw.end(), image.rgb.begin() + 3 * y * image.width, image.rgb.begin() + 3 * (y + 1) * image.width);
  }
  Bytes z{0x78, 0x01};
  for (size_t p = 0; p < raw.size();) {
    unsigned n = std::min<size_t>(65535, raw.size() - p), inverse = (~n) & 0xffff;
    z.push_back(p + n == raw.size() ? 1 : 0);
    z.push_back(n & 255);
    z.push_back(n >> 8);
    z.push_back(inverse & 255);
    z.push_back(inverse >> 8);
    z.insert(z.end(), raw.begin() + p, raw.begin() + p + n);
    p += n;
  }
  std::uint32_t a = 1, b = 0;
  for (auto v : raw) {
    a = (a + v) % 65521;
    b = (b + a) % 65521;
  }
  big_endian(z, (b << 16) | a);
  chunk(png, "IDAT", z);
  chunk(png, "IEND", {});
  std::ofstream f(file, std::ios::binary);
  f.exceptions(std::ios::badbit | std::ios::failbit);
  f.write(reinterpret_cast<const char*>(png.data()), png.size());
  f.close();
}
}  // namespace

namespace {
Image render(const HhdsDesign& initial, const HhdsDesign& placed, const std::string& title, const Result& result,
             Initialization initialization, const std::string& phase, const std::string& footer1,
             const std::string& footer2) {
  initial.design.validate();
  placed.design.validate();
  Image image(1600, 900);
  image.rect(0, 0, 1600, 68, ink);
  text(image, 28, 16, "HHDS / " + title, {255, 255, 255}, 3);
  text(image, 28, 46,
       std::to_string(placed.cells.size()) + " CELLS | " + std::to_string(placed.net_drivers.size()) +
           " NETS | UNIFORM WEIGHTS | NO TIMING MODEL",
       {203, 218, 237}, 2);
  text(image, 40, 94, initialization == Initialization::mesh_grid ? "INITIAL MESH GRID" : "INITIAL RANDOM SEED", ink,
       2);
  text(image, 832, 94, phase, {36, 123, 88}, 2);
  text(image, 40, 122, "HPWL " + number(result.initial_hpwl) + " UM");
  text(image, 832, 122, "HPWL " + number(result.hpwl) + " UM | LEGAL " + (result.legal ? "YES" : "NO"));
  panel(image, initial, 32);
  panel(image, placed, 824);
  text(image, 40, 846, title.starts_with("tree") ? "CELL COLOR = TREE STAGE" : "CELL COLOR = LOGICAL GRID COLUMN");
  text(image, 40, 872, "LINES = NETS / SQUARES = IO / HHDS POSITIONS", {92, 108, 128});
  text(image, 832, 846, footer1);
  text(image, 832, 872, footer2, {92, 108, 128});
  return image;
}
}  // namespace

void write_png(const HhdsDesign& initial, const HhdsDesign& placed, const std::filesystem::path& file,
               const std::string& title, const Result& result, Initialization initialization) {
  encode(render(initial, placed, title, result, initialization, "NATIVE GP + LEGALIZATION + DETAIL",
                "GP OVERFLOW " + number(result.overflow * 100) + "% / " + std::to_string(result.iterations) + " STEPS",
                result.converged ? "GP THRESHOLD REACHED" : "GP LIMIT REACHED - CHECK QOR.JSON"),
         file);
}

void write_frame_png(const HhdsDesign& initial, const HhdsDesign& current, const std::filesystem::path& file,
                     const std::string& title, Initialization initialization, const AnimationFrame& frame,
                     size_t overlap_pairs) {
  Result result;
  result.initial_hpwl = hpwl(initial.design);
  result.hpwl = hpwl(current.design);
  result.legal = is_legal(current.design);
  auto full =
      render(initial, current, title, result, initialization, stage_name(frame.stage),
             "GP STEP " + std::to_string(frame.iteration), "OVERLAPPING PAIRS " + std::to_string(overlap_pairs));
  // Half-size frames keep an optional trace small (about 1.1MB per PNG).
  Image small(800, 450);
  for (int y = 0; y < small.height; ++y)
    for (int x = 0; x < small.width; ++x)
      for (int c = 0; c < 3; ++c) {
        int p = 3 * (2 * y * full.width + 2 * x) + c;
        small.rgb[3 * (y * small.width + x) + c] =
            (full.rgb[p] + full.rgb[p + 3] + full.rgb[p + 3 * full.width] + full.rgb[p + 3 * full.width + 3] + 2) / 4;
      }
  encode(small, file);
}
}  // namespace DreamPlace::native::synthetic
