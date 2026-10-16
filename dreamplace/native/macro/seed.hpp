// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/dataflow.hpp"
namespace DreamPlace::native::macro {
Topology seed_topology(const Design &, std::span<const NodeInfo>,
                       std::span<const Flow>, const Parameters &,
                       const Topology &input = {});
} // namespace DreamPlace::native::macro
