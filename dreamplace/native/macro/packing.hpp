// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"
namespace DreamPlace::native::macro {
// Deterministic conservative rectangle packing for a group-only fallback.
// No rotations. Fixed positive-area objects are obstacles. Halos and occupancy
// round outwards to sites/rows. A false result leaves the design unchanged;
// failure of this heuristic is not proof that the outline is infeasible.
bool pack_rectangles(Design &);
} // namespace DreamPlace::native::macro
