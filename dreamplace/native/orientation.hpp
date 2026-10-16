// BSD-3-Clause. See LICENSE.
#pragma once
#include <utility>

namespace DreamPlace::native {
struct Design;
// Bit 0 mirrors X (FN); bit 1 mirrors Y (FS). No axis exchange or rotation.
enum class Orientation : unsigned char { N = 0, FN = 1, FS = 2, S = 3 };
bool valid_orientation(Orientation orientation);
Orientation compose(Orientation outer, Orientation inner);
std::pair<double, double> orient_pin(double x, double y, double width,
                                     double height, Orientation orientation);
// pin_offset_x/y always contain offsets in the current orientation. Transform
// from the previous orientation, avoiding a second transform on repeated calls.
// Explicit operation, including on fixed nodes; place() itself preserves them.
// Invalid input leaves the design unchanged.
void set_orientation(Design &design, int node, Orientation orientation);
} // namespace DreamPlace::native
