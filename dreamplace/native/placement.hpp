// BSD-3-Clause. See LICENSE.
#pragma once

#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

#include "dreamplace/native/constraints.hpp"
#include "dreamplace/native/orientation.hpp"

namespace DreamPlace::native {

struct Layout {
  double xl = 0, yl = 0, xh = 0, yh = 0;
  double site_width = 1, row_height = 1;
};

// Owning flat-array boundary. Coordinates/sizes/offsets use one caller-chosen
// unit, normally DBU. Order: movable, fixed, zero-area IO terminals. Fillers
// are internal to global placement and never exposed in the result.
// Disjoint color/region groups. Empty areas use physical width*height.
// Optional areas preserve physical weighting when occupancy is halo padded.
struct CohesionGroup {
  std::vector<int> nodes;
  double weight = 1;
  std::vector<double> areas;
};
struct Design {
  Layout layout;
  std::vector<CohesionGroup> cohesion_groups;
  // Explicit non-standard-cell identity, independent of geometry/row height.
  std::vector<std::uint8_t> node_macro;
  // Sum of track densities on LEF layers obstructed by each macro. RUDY
  // conservatively blocks those layers over its complete occupancy envelope.
  std::vector<double> node_blocked_capacity_h, node_blocked_capacity_v;
  int num_movable = 0;
  std::vector<double> node_x, node_y, node_size_x, node_size_y;
  // Optional per-node clearance on each side, in layout units. Empty means
  // zero. Rounded up to whole sites/rows in private occupancy envelopes;
  // physical dimensions and pin offsets are never changed.
  std::vector<double> node_halo_x, node_halo_y;
  // Empty means N. Pin offsets already reflect this orientation; GP never
  // transforms them implicitly. set_orientation() updates both together.
  std::vector<Orientation> node_orientation;
  std::vector<PositionBox>
      node_boxes; // optional, one physical-origin box per node
  std::vector<RelativeOrder> order_edges;
  std::vector<std::uint8_t> node_inflatable; // optional, standard cells only
  std::vector<int> pin2node, pin2net;
  std::vector<double> pin_offset_x, pin_offset_y;
  std::vector<int> flat_net2pin, flat_net2pin_start{0};
  std::vector<int> flat_node2pin, flat_node2pin_start;
  std::vector<double> net_weights;
  std::vector<std::uint8_t> net_mask;

  // Rebuild inverse maps from net CSR and pin2node. Rejects duplicate/missing
  // pins and malformed CSR before kernels can index raw pointers.
  void build_maps();
  void validate() const;
};

enum class PlacementStage {
  initial,
  global_placement,
  legalization,
  detailed_placement,
  finished
};
const char *stage_name(PlacementStage stage);
// Synchronous, read-only physical snapshots. References are valid only during
// the call; temporary fillers are excluded. Exceptions abort without writeback.
using PlacementObserver =
    std::function<void(const Design &, PlacementStage, int iteration)>;

struct GpSample;
struct Options {
  int max_iterations = 500;
  int bins_x = 32, bins_y = 32;
  // Opt-in repulsion controls. Automatic bins resolve before halo padding.
  bool automatic_bins = false, exact_macro_density = false;
  bool adaptive_density = false, improved_routability = false;
  bool scale_macro_density = false;
  std::function<void(const GpSample &)> progress_observer;
  double macro_overlap_weight = 0,
         macro_jitter = 0; // jitter fraction of macro size
  int num_threads = 1;
  std::uint32_t seed = 1000; // retains the original filler RNG sequence
  double target_density = 0.8;
  double stop_overflow = 0.1;
  double initial_density_weight = 0.01;
  double max_step = 0.5; // fraction of a bin
  // Optional H2D hooks. Zero/default controls preserve the original GP.
  double filler_width = 0;
  // Soft centroid cohesion, inverse coordinate units. Lambda ramps linearly
  // from zero to this value over schedule_iterations (or max_iterations).
  // Zero disables it exactly; no boxes or changes to legalization.
  double cohesion_weight = 0;
  double order_weight = 0, order_growth = 1;
  int schedule_iterations = 0;
  bool routability = false;
  double route_capacity_h = 1, route_capacity_v = 1;
  bool add_fillers = true;
  bool legalize = true;
  bool detailed_place = true;
  PlacementObserver observer;
  int observer_interval = 20; // GP steps between snapshots, when enabled
};

struct GpSample {
  int iteration = 0;
  double overflow = 0, hpwl = 0, density_weight = 0;
  double route_overflow = 0, inflation_area_ratio = 1;
  std::uint64_t macro_overlaps = 0, macro_stacks = 0;
  double raw_overflow = 0;
  double model_density = 0, filler_charge_scale = 1;
};
struct Result {
  int bins_x = 0, bins_y = 0;
  double poisson_seconds = 0;
  std::vector<GpSample> history;
  int iterations = 0;
  double initial_hpwl = 0, hpwl = 0, overflow = 0;
  bool converged = false;
  bool legal = false;
};

// Optional full GP continuation. Includes private fillers, Nesterov/BB history,
// density normalization and global iteration offset. A mismatching design is
// rejected; failed calls leave state untouched. Empty state starts a new solve.
struct GpState {
  std::uint64_t signature = 0;
  int iteration = 0;
  double density_weight = 0, momentum = 1, step = 0;
  double previous_overflow = -1, previous_hpwl = 0;
  std::vector<double> x, y, prev_x, prev_y, old_gx, old_gy, old_ex, old_ey;
  std::vector<double> inflation;
};

// Throws invalid_argument for bad inputs and runtime_error if legalization
// cannot produce a legal result. Failed calls leave the caller's design intact.
// Non-movable nodes remain bit-identical. No Torch, Python, CUDA, or file I/O.
Result place(Design &design, const Options &options = {},
             GpState *state = nullptr);
double hpwl(const Design &design);
bool is_legal(const Design &design);

// Reusable per-instance configuration. The caller owns Design (and any HHDS
// mapping); every call creates fresh solver scratch state and updates only the
// design's coordinates on success. Options are validated when place() runs.
// Independent instances/designs may run concurrently. Configure options before
// launching work; callbacks and shared host/HHDS state require caller ownership
// or synchronization. num_threads controls workers within each call.
class Placer {
public:
  explicit Placer(Options options = {}) : options_(std::move(options)) {}

  Options &options() noexcept { return options_; }
  const Options &options() const noexcept { return options_; }
  Result place(Design &design) const {
    return ::DreamPlace::native::place(design, options_);
  }

private:
  Options options_;
};

namespace detail {
struct Gradient {
  double value = 0;
  std::vector<double> x, y;
};
Gradient cohesion(const Design &design, double lambda);
// Squared exact intersection area, normalized by pair size for unit invariance.
Gradient macro_overlap(const Design &, double weight);
std::pair<std::uint64_t, std::uint64_t> macro_overlap_counts(const Design &);
std::pair<int, int> automatic_bins(const Design &);
Gradient wirelength(const Design &design, double gamma, int threads);
void density_map(const Design &d, int bx, int by, int threads,
                 const std::vector<double> &sx, const std::vector<double> &sy,
                 const std::vector<double> &ox, const std::vector<double> &oy,
                 const std::vector<double> &ratio, std::vector<double> &map);
Gradient
density_force(const Design &d, int bx, int by, int threads,
              const std::vector<double> &sx, const std::vector<double> &sy,
              const std::vector<double> &ox, const std::vector<double> &oy,
              const std::vector<double> &ratio, const std::vector<double> &fx,
              const std::vector<double> &fy);
void legalize(Design &d, int bx, int by);
void detailed_place(Design &d, int bx, int by, int threads);
bool kernel_is_legal(const Design &d);
} // namespace detail
} // namespace DreamPlace::native
