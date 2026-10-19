// BSD-3-Clause. See LICENSE.
#pragma once
#include <cstdint>
#include "utility/src/namespace.h"
DREAMPLACE_BEGIN_NAMESPACE
struct MacroAlignConstraint {
  int a, b, axis;
  std::int64_t weight;
  // Soft signed displacement: coordinate[b] - coordinate[a] = offset.
  double offset = 0;
};
DREAMPLACE_END_NAMESPACE
