// BSD-3-Clause. See LICENSE.
#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <span>
#include <vector>

namespace DreamPlace::native {
struct Design;
inline constexpr std::uint32_t uroute_feature_version = 1;
inline constexpr std::uint32_t uroute_label_version = 1;
inline constexpr std::uint32_t uroute_arch_version = 1;
// Architecture v1 network contract: support/query grids and channel counts.
inline constexpr int uroute_network_bins = 256;
inline constexpr int uroute_feature_channels = 4, uroute_label_channels = 3;

// Channel-major, then x-major: (channel * bins_x + x) * bins_y + y.
// Features are normalized independently per channel. Labels stay absolute.
// Constant channels normalize to zero; ranges retain their physical value.
struct UrouteMaps {
  int channels = 0, bins_x = 0, bins_y = 0;
  bool normalized = false;
  std::vector<float> values;
  std::vector<std::array<double, 2>> ranges;
};

// Physical geometry only: no WDM/timing weights, obstacles, halos or fillers.
// Compute on the route grid, area-resample to the network grid, then normalize.
// Degenerate net boxes contribute zero, as in legacy bounding-box RUDY.
UrouteMaps uroute_features(const Design &, int route_bins_x, int route_bins_y,
                           int output_bins_x = uroute_network_bins,
                           int output_bins_y = uroute_network_bins);
// Piecewise-constant area averaging on coincident rectangular grids.
std::vector<double> uroute_resample(std::span<const double>, int from_x,
                                    int from_y, int to_x, int to_y);
// Versioned, little-endian format shared with the offline tools. Readers reject
// nonfinite values, oversized shapes, inconsistent ranges and trailing bytes.
void write_uroute_maps(const std::filesystem::path &, const UrouteMaps &);
UrouteMaps read_uroute_maps(const std::filesystem::path &);
} // namespace DreamPlace::native
