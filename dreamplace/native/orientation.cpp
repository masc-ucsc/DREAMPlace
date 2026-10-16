// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/orientation.hpp"
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <stdexcept>

namespace DreamPlace::native {
bool valid_orientation(Orientation orientation) {
  return static_cast<unsigned>(orientation) <= 3;
}
Orientation compose(Orientation outer, Orientation inner) {
  if (!valid_orientation(outer) || !valid_orientation(inner))
    throw std::invalid_argument("invalid flip orientation");
  return static_cast<Orientation>(static_cast<unsigned>(outer) ^
                                  static_cast<unsigned>(inner));
}
std::pair<double, double> orient_pin(double x, double y, double width,
                                     double height, Orientation orientation) {
  if (!valid_orientation(orientation) || !std::isfinite(x) ||
      !std::isfinite(y) || !std::isfinite(width) || !std::isfinite(height) ||
      width < 0 || height < 0)
    throw std::invalid_argument("invalid oriented pin geometry");
  const auto bits = static_cast<unsigned>(orientation);
  x = bits & 1 ? width - x : x;
  y = bits & 2 ? height - y : y;
  if (!std::isfinite(x) || !std::isfinite(y))
    throw std::invalid_argument("oriented pin offset overflow");
  return {x, y};
}
void set_orientation(Design &d, int node, Orientation orientation) {
  d.validate();
  if (node < 0 || node >= static_cast<int>(d.node_x.size()) ||
      !valid_orientation(orientation))
    throw std::invalid_argument("invalid node orientation");
  if ((d.node_size_x[node] == 0 || d.node_size_y[node] == 0) &&
      orientation != Orientation::N)
    throw std::invalid_argument("terminal cannot be flipped");
  const auto old =
      d.node_orientation.empty() ? Orientation::N : d.node_orientation[node];
  const auto delta = compose(orientation, old);
  auto x = d.pin_offset_x, y = d.pin_offset_y;
  for (int k = d.flat_node2pin_start[node]; k < d.flat_node2pin_start[node + 1];
       ++k) {
    const int p = d.flat_node2pin[k];
    const auto pin =
        orient_pin(x[p], y[p], d.node_size_x[node], d.node_size_y[node], delta);
    x[p] = pin.first;
    y[p] = pin.second;
  }
  auto orientations = d.node_orientation;
  if (orientations.empty())
    orientations.assign(d.node_x.size(), Orientation::N);
  orientations[node] = orientation;
  d.pin_offset_x.swap(x);
  d.pin_offset_y.swap(y);
  d.node_orientation.swap(orientations);
}
} // namespace DreamPlace::native
