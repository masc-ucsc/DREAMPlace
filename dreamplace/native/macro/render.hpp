// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/topology.hpp"
#include <filesystem>
namespace DreamPlace::native::macro {
// One final physical layout: cyan cells, purple macros, group/halo outlines,
// yellow anchors, and green/red satisfied/violated relative-order arrows.
// Pure PNG encoder, no GUI/Torch/HHDS/plotting dependency.
void write_png(const std::filesystem::path &, const Design &,
               std::span<const NodeInfo>, const Topology &, int width = 1200);
} // namespace DreamPlace::native::macro
