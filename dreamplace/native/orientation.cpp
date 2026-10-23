// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/orientation.hpp"
#include "dreamplace/native/placement.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

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
namespace {
// Touches only the node's own pins; every check and transform happens before
// the first write, so invalid input leaves the design unchanged.
void apply_orientation(Design &d, int node, Orientation orientation) {
  if (node < 0 || node >= static_cast<int>(d.node_x.size()) ||
      !valid_orientation(orientation))
    throw std::invalid_argument("invalid node orientation");
  if ((d.node_size_x[node] == 0 || d.node_size_y[node] == 0) &&
      orientation != Orientation::N)
    throw std::invalid_argument("terminal cannot be flipped");
  const auto old =
      d.node_orientation.empty() ? Orientation::N : d.node_orientation[node];
  if (!d.node_orientation_fixed.empty() && d.node_orientation_fixed[node] &&
      orientation != old)
    throw std::invalid_argument("rigid group orientation is frozen");
  const auto delta = compose(orientation, old);
  const int begin = d.flat_node2pin_start[node],
            end = d.flat_node2pin_start[node + 1];
  std::vector<std::pair<double, double>> pins;
  pins.reserve(end - begin);
  for (int k = begin; k < end; ++k) {
    const int p = d.flat_node2pin[k];
    pins.push_back(orient_pin(d.pin_offset_x[p], d.pin_offset_y[p],
                              d.node_size_x[node], d.node_size_y[node],
                              delta));
  }
  if (d.node_orientation.empty())
    d.node_orientation.assign(d.node_x.size(), Orientation::N);
  d.node_orientation[node] = orientation;
  for (int k = begin; k < end; ++k) {
    const int p = d.flat_node2pin[k];
    d.pin_offset_x[p] = pins[k - begin].first;
    d.pin_offset_y[p] = pins[k - begin].second;
  }
}
} // namespace
void set_orientation(Design &d, int node, Orientation orientation) {
  d.validate();
  apply_orientation(d, node, orientation);
}
void set_orientation_unchecked(Design &d, int node, Orientation orientation) {
  apply_orientation(d, node, orientation);
}
} // namespace DreamPlace::native
