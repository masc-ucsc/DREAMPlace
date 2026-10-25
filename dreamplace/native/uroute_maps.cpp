// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/uroute.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace DreamPlace::native {
namespace {
template <class T> void put(std::ostream &out, T value) {
  auto bytes = std::bit_cast<std::array<char, sizeof(T)>>(value);
  if constexpr (std::endian::native == std::endian::big)
    std::reverse(bytes.begin(), bytes.end());
  out.write(bytes.data(), bytes.size());
}
template <class T> T get(std::istream &in) {
  std::array<char, sizeof(T)> bytes;
  if (!in.read(bytes.data(), bytes.size()))
    throw std::runtime_error("truncated URoute map");
  if constexpr (std::endian::native == std::endian::big)
    std::reverse(bytes.begin(), bytes.end());
  return std::bit_cast<T>(bytes);
}
void shape(const UrouteMaps &m) {
  if (m.channels < 1 || m.channels > 16 || m.bins_x < 1 || m.bins_x > 1024 ||
      m.bins_y < 1 || m.bins_y > 1024)
    throw std::invalid_argument("invalid URoute map shape");
}
void validate(const UrouteMaps &m) {
  shape(m);
  if (m.values.size() != size_t(m.channels) * m.bins_x * m.bins_y ||
      m.ranges.size() != size_t(m.channels))
    throw std::invalid_argument("URoute map shape mismatch");
  const size_t area = size_t(m.bins_x) * m.bins_y;
  for (int c = 0; c < m.channels; ++c) {
    const auto [lo, hi] = m.ranges[c];
    if (!std::isfinite(lo) || !std::isfinite(hi) || hi < lo ||
        !std::isfinite(hi - lo))
      throw std::invalid_argument("invalid URoute normalization range");
    for (size_t i = c * area; i < (c + 1) * area; ++i) {
      const double v = m.values[i];
      if (!std::isfinite(v) ||
          (m.normalized && (v < 0 || v > 1 || (lo == hi && v != 0))) ||
          (!m.normalized &&
           (v < static_cast<float>(lo) || v > static_cast<float>(hi))))
        throw std::invalid_argument("invalid URoute map value");
    }
  }
}
} // namespace
void write_uroute_maps(const std::filesystem::path &path, const UrouteMaps &m) {
  validate(m);
  std::ofstream out(path, std::ios::binary);
  out.write("UROUTEM1", 8);
  put<std::uint32_t>(out, 1);
  put<std::uint32_t>(out, m.channels);
  put<std::uint32_t>(out, m.bins_x);
  put<std::uint32_t>(out, m.bins_y);
  put<std::uint32_t>(out, m.normalized);
  for (const auto &range : m.ranges)
    for (double v : range)
      put(out, v);
  for (float v : m.values)
    put(out, v);
  out.close();
  if (!out)
    throw std::runtime_error("cannot write URoute maps: " + path.string());
}
UrouteMaps read_uroute_maps(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  std::array<char, 8> magic{};
  if (!in.read(magic.data(), magic.size()) ||
      std::string(magic.data(), magic.size()) != "UROUTEM1" ||
      get<std::uint32_t>(in) != 1)
    throw std::runtime_error("unsupported URoute map format: " + path.string());
  UrouteMaps m;
  m.channels = get<std::uint32_t>(in);
  m.bins_x = get<std::uint32_t>(in);
  m.bins_y = get<std::uint32_t>(in);
  const auto normalized = get<std::uint32_t>(in);
  if (normalized > 1)
    throw std::runtime_error("invalid URoute normalization flag");
  m.normalized = normalized;
  shape(m);
  m.ranges.resize(m.channels);
  for (auto &range : m.ranges)
    for (double &v : range)
      v = get<double>(in);
  m.values.resize(size_t(m.channels) * m.bins_x * m.bins_y);
  for (float &v : m.values)
    v = get<float>(in);
  if (in.peek() != std::char_traits<char>::eof())
    throw std::runtime_error("trailing URoute map data");
  validate(m);
  return m;
}
} // namespace DreamPlace::native
