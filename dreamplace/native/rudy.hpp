// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"

namespace DreamPlace::native {
enum class RudyModel { bbox, rsmt };
struct RudyOptions {
  int bins_x = 32, bins_y = 32;
  // Routing tracks per layout unit; callers derive these from their technology.
  double horizontal_capacity = 1, vertical_capacity = 1;
  bool macro_blockages = false, pin_access = false;
  RudyModel model = RudyModel::bbox;
  int degree = 40, accuracy = 3;
  // Caller loads FLUTE tables once. Required for degree >= 4 tree nets.
  std::function<void()> initialize;
};
struct RudyResult {
  // X-major bin index = x*bins_y+y. Directional utilization is demand
  // wirelength divided by bin area and directional track density. overflow is
  // the mean excess of max(horizontal,vertical,pin_utilization) over 1 across
  // all bins, plus unreachable-net/blocked-pin penalties in enhanced mode.
  std::vector<double> horizontal, vertical, pin_utilization;
  double overflow = 0, max_utilization = 0;
  int measured_nets = 0, degenerate_nets = 0;
  // Enhanced mode: exact footprint union, no routing capacity inside macros.
  std::vector<double> blocked_fraction;
  int unroutable_nets = 0, blocked_pins = 0;
  double detour_wirelength = 0;
};
// Serial deterministic RUDY measurement; no placement edits or cell inflation.
// Uses DREAMPlace's WDM degree weights and enabled physical hypernets. Like
// upstream's projected-box formula, zero-area net bboxes contribute zero demand
// and are counted separately in legacy mode. macro_blockages enables opaque
// physical macro obstacles (halo whitespace stays available), free-area demand
// redistribution and corridor detours; straight nets consume routing tracks in
// that mode. It remains a proxy, not a detailed router. Enhanced mode clips IOs
// to the core boundary; legacy mode clips bin overlap but keeps the full net
// bbox normalization.
RudyResult rudy(const Design &design, const RudyOptions &options = {});
namespace detail {
struct MacroRoutingPressure {
  // Objective gradients: descent opens congested channels. Fixed nodes stay
  // zero.
  std::vector<double> x, y;
  std::uint64_t channels = 0;
};
// Utilization is an x-major RUDY map. Pressure is quadratic below capacity too,
// capped at utilization 3, and tapers to zero at range_bins. Physical bodies,
// not halos, define channels; third macros occlude non-neighboring edges.
MacroRoutingPressure macro_routing_pressure(const Design &, int bins_x,
                                            int bins_y,
                                            const std::vector<double> &,
                                            double range_bins);
} // namespace detail
} // namespace DreamPlace::native
