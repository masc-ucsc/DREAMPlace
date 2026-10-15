// BSD-3-Clause. See LICENSE.
#pragma once
#include "utility/src/namespace.h"
DREAMPLACE_BEGIN_NAMESPACE
enum Direction1DType { kLOW = 0, kHIGH = 1, kX = 0, kY = 1, kLEFT = 0, kRIGHT = 1, kBOTTOM = 0, kTOP = 1 };
enum Direction2DType { kXLOW = 0, kXHIGH = 1, kYLOW = 2, kYHIGH = 3 };

/// extract x/y information from Direction2DType
inline Direction1DType getXY(Direction2DType d) { return Direction1DType(d > 1); }
/// extract low/high information from Direction2DType
inline Direction1DType getLH(Direction2DType d) { return Direction1DType(d & 1); }
/// construct Direction2DType from x/y and low/high information
inline Direction2DType to2D(Direction1DType xy, Direction1DType lh) {
  return Direction2DType(((int)xy << 1) + (int)lh);
}

DREAMPLACE_END_NAMESPACE
