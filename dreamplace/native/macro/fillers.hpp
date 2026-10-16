// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"
namespace DreamPlace::native::macro {
// Keep the trimmed standard-cell width when possible. On sparse cells-only
// designs, widen private filler rectangles to respect the native node budget.
// Physical dimensions and requested density are unaffected.
double bounded_filler_width(const Design &, const Options &, double preferred);
} // namespace DreamPlace::native::macro
