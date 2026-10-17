// BSD-3-Clause. See LICENSE.
#pragma once
#include <cstdint>
#include "utility/src/namespace.h"
DREAMPLACE_BEGIN_NAMESPACE
struct MacroAlignConstraint {
  int a, b, axis;
  std::int64_t weight;
};
DREAMPLACE_END_NAMESPACE
