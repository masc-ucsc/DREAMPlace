// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/placement.hpp"

namespace DreamPlace::native {
struct RudyOptions {
  int bins_x = 32, bins_y = 32;
  // Routing tracks per layout unit; callers derive these from their technology.
  double horizontal_capacity = 1, vertical_capacity = 1;
  bool macro_blockages = false, pin_access = false;
};
struct RudyResult {
  // X-major bin index = x*bins_y+y. Directional utilization is demand
  // wirelength divided by bin area and directional track density. overflow is
  // the mean excess of max(horizontal,vertical) over 1 across all bins.
  std::vector<double> horizontal, vertical, pin_utilization;
  double overflow = 0, max_utilization = 0;
  int measured_nets = 0, degenerate_nets = 0;
};
// Serial deterministic RUDY measurement; no placement edits or cell inflation.
// Uses DREAMPlace's WDM degree weights and enabled physical hypernets. Like
// upstream's projected-box formula, zero-area net bboxes contribute zero demand
// and are counted separately. IO coordinates outside the core are clipped when
// accumulating overlap, without changing the net's full bbox normalization.
RudyResult rudy(const Design &design, const RudyOptions &options = {});
} // namespace DreamPlace::native
