// BSD-3-Clause. See LICENSE.
#pragma once
#include <string>
namespace DreamPlace::native {
enum class Axis { x, y };
enum class Side { left, right, bottom, top };
struct Strength {
  bool hard = false;
  double weight = 1;
};
// Lower-left position bounds (not a footprint box).
struct PositionBox {
  double xl, yl, xh, yh;
};
// Transient solver indices; persistent Topology constraints use names.
struct RelativeOrder {
  int from = 0, to = 0;
  Axis axis = Axis::x;
  double separation = 0;
  Strength strength;
  std::string id;
};
} // namespace DreamPlace::native
