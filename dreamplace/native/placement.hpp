// BSD-3-Clause. See LICENSE.
#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "dreamplace/native/constraints.hpp"
#include "dreamplace/native/orientation.hpp"
#include "macro_legalize/src/alignment.h"

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
// Candidate colors may overlap across any number of groups. Nodes stores only
// uncolored descendants; each color includes its full physical membership.
struct CentroidGroup {
  std::string name;
  int parent = -1, family = 0, scope_macros = 0, scope_cells = 0;
  double equivalent_cells = 0;
  std::vector<int> nodes, colors;
  std::string copy_parent{}; // immediate hierarchy parent used for replication
  int copies_per_parent = 0;
  // -1 keeps same-family repulsion. Nonnegative keys restrict repulsion to
  // one candidate hierarchy scope, independently of snap/type families.
  int repulsion_scope = -1;
};
struct CentroidSample {
  double x = 0, y = 0, radius = 0;
  double area = 0; // physical anchored area; fillers/halo padding excluded
  int cells = 0, macros = 0;
};
struct Design {
  Layout layout;
  std::vector<CohesionGroup> cohesion_groups;
  std::vector<CentroidGroup> centroid_groups;
  // Optional frozen world-coordinate targets, one per candidate. Empty keeps
  // the original moving-centroid attraction exactly. Targets are not fences.
  std::vector<std::pair<double, double>> centroid_targets;
  std::vector<int> node_centroid_color; // -1 for uncolored physical nodes
  std::vector<std::string> centroid_color_names; // scoped host color identities
  std::vector<double> centroid_node_area; // physical weights, not halo area
  // Explicit non-standard-cell identity, independent of geometry/row height.
  std::vector<std::uint8_t> node_macro;
  // Same master x synthesis-color identity; -1 opts a node out.
  std::vector<int> node_align_class;
  // Sum of track densities on LEF layers obstructed by each macro. RUDY
  // conservatively blocks those layers over its complete occupancy envelope.
  std::vector<double> node_blocked_capacity_h, node_blocked_capacity_v;
  int num_movable = 0;
  std::vector<double> node_x, node_y, node_size_x, node_size_y;
  // Optional per-node clearance on each side, in layout units. Empty means
  // zero. Rounded up to whole sites/rows in private occupancy envelopes;
  // physical dimensions and pin offsets are never changed.
  std::vector<double> node_halo_x, node_halo_y;
  // Internal: the padding place() applied to node_size on each side when it
  // solves a halo design. GP ignores it except to recover physical areas after
  // a checkpoint resize. Callers leave it empty; place() never returns it set.
  std::vector<double> node_pad_x, node_pad_y;
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
  centroid_start,
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

// Checkpoints update the objective/standard cells. A replacement changes the
// standard-cell netlist; previous_node maps each replacement node to its old
// index (-1 for an inserted movable standard cell). Existing positions and
// macro/fixed geometry are preserved; momentum restarts after topology changes.
// Empty vectors preserve fields; do not combine them with replacement.
struct GpUpdate {
  // Optional nonnegative setup violation from the host STA at this checkpoint.
  std::optional<double> timing_violation;
  std::vector<double> net_weights, node_size_x, node_size_y;
  std::vector<double> pin_offset_x, pin_offset_y;
  std::optional<Design> replacement;
  std::vector<int> previous_node;
};
using GpCheckpoint = std::function<GpUpdate(const Design &, int)>;
struct GpSample;
struct GpCheckpointSample;
struct Options {
  int checkpoint_interval = 0;
  bool checkpoint_final =
      false; // also sample the iteration cap, without double-calling
  // Defer checkpoints above stop_overflow. Always invoke at the iteration or
  // continuation cap, and never declare convergence before the first callback.
  // Uses the same physical (halo-padded) overflow as GP, before host unpadding.
  bool checkpoint_after_density = false;
  std::function<void(const GpCheckpointSample &)> checkpoint_observer;
  GpCheckpoint checkpoint;

  int max_iterations = 500;
  int bins_x = 32, bins_y = 32;
  // Opt-in repulsion controls. Automatic bins resolve before halo padding.
  bool automatic_bins = false, exact_macro_density = false;
  bool exact_cell_density =
      false; // opt-in physical footprint for movable cells
  bool adaptive_density = false, improved_routability = false;
  bool scale_macro_density = false;
  bool routability_charge_compensation =
      true; // false replays the earlier R4 charge model
  std::function<void(const GpSample &)> progress_observer;
  double wirelength_gamma_end = 0.05; // final/initial smoothing ratio
  double fixed_macro_overlap_weight =
      0; // exact cell/fixed-macro edge penalty; opt-in
  double macro_overlap_weight = 0,
         macro_jitter = 0; // jitter fraction of macro size
  double macro_align_weight = 0, macro_align_sigma = 0;
  double macro_align_overlap_boost = 0;
  int macro_align_start = 100, macro_align_full = 500;
  bool macro_align_nearest = false, macro_snap = false;
  double macro_snap_tolerance = .25;
  std::int64_t macro_snap_weight = 100;
  // Zero disables R8H2 exactly. Attraction starts in the first GP iteration;
  // no preliminary macro placement or extra continuation phase is used.
  double centroid_weight = 0;
  // Relative to centroid attraction. Same-family circles have area A/density.
  // Zero preserves the attraction-only path exactly.
  double centroid_repulsion_weight = 0;
  // Increase the circle penalty until separation and density both settle.
  // Off by default; relative penetration tolerance is measured against r_i+r_j.
  bool centroid_adaptive_repulsion = false;
  double centroid_stop_overlap = .001;
  int centroid_ramp = 100;
  int num_threads = 1;
  std::uint32_t seed = 1000; // retains the original filler RNG sequence
  double target_density = 0.8;
  double stop_overflow = 0.1;
  int min_iterations = 40; // floor before considering density convergence
  // Zero preserves density-only stopping. Otherwise allow factor * first
  // density-feasible iteration ADDITIONAL iterations, bounded by
  // max_iterations.
  int post_density_factor = 0;
  double convergence_tolerance = .01; // relative weighted-HPWL/inflation change
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
  // Experimental filler-free upper-bound penalty. Underfull bins contribute
  // no charge; fixed obstacles alone must not push cells into whitespace.
  bool density_upper_bound = false;
  bool add_fillers = true;
  bool legalize = true;
  bool detailed_place = true;
  int detailed_iterations =
      2; // iteration cap per existing detailed-placement kernel
  PlacementObserver observer;
  int observer_interval = 20; // GP steps between snapshots, when enabled
};

struct AlignmentPair {
  int a, b;
  Axis axis;
};
struct AlignmentGroup {
  int align_class;
  Axis axis;
  double target;
  std::vector<int> nodes;
};
struct AlignmentSnap {
  std::vector<AlignmentGroup> groups;
  std::vector<::DreamPlace::MacroAlignConstraint> constraints;
};
struct AlignmentReport {
  std::vector<AlignmentGroup> groups;
  std::uint64_t x_groups = 0, y_groups = 0, members = 0, both_axes = 0;
  std::uint64_t exact_groups = 0, broken_groups = 0, soft_violations = 0;
  std::vector<double> macro_displacements;
  double macro_legalization_seconds = 0, cell_legalization_seconds = 0;
};
struct GpSample {
  int iteration = 0;
  double overflow = 0, hpwl = 0, density_weight = 0;
  double route_overflow = 0, inflation_area_ratio = 1;
  std::uint64_t macro_overlaps = 0, macro_stacks = 0;
  double raw_overflow = 0;
  double model_density = 0, filler_charge_scale = 1;
  std::uint64_t align_pairs = 0;
  double align_residual = 0, align_ramp = 0;
  std::uint64_t centroid_overlaps = 0;
  double centroid_max_overlap = 0, centroid_repulsion_weight = 0;
};
struct CentroidReport {
  bool started = false, completed = false;
  int start_iteration = 0, iterations = 0;
  double force_seconds = 0;
  bool circles_converged = false;
  std::uint64_t circle_overlaps = 0;
  double max_circle_overlap = 0, effective_repulsion_weight = 0;
};
struct DetailedSample {
  std::string kernel;
  double before = 0, after = 0;
  bool legal = false, macros_preserved = false, accepted = false;
  std::string legality_issue;
};
struct GpConvergence {
  int density_iteration = 0, deadline = 0, stable_samples = 0;
  double previous_hpwl = -1, previous_inflation = -1;
  bool timing_observed = false, timing_satisfied = true;
};
struct GpCheckpointSample {
  int iteration = 0;
  double overflow = 0;
  bool invoked = false, forced = false;
};
struct Result {
  // Populated only with checkpoint_after_density; includes skipped intervals.
  std::vector<GpCheckpointSample> checkpoint_history;
  int density_iteration = 0, continuation_deadline = 0;
  bool objective_converged = false, continuation_limit = false;
  bool timing_observed = false, timing_satisfied = true;
  CentroidReport centroid;
  AlignmentReport alignment;
  int bins_x = 0, bins_y = 0;
  double poisson_seconds = 0;
  double global_placement_seconds = 0, detailed_placement_seconds = 0;
  std::vector<GpSample> history;
  std::vector<DetailedSample> detailed_history;
  int iterations = 0, selected_iteration = 0;
  double initial_hpwl = 0, hpwl = 0, overflow = 0;
  bool converged = false;
  bool legal = false;
};

// Optional full GP continuation. Includes private fillers, Nesterov/BB history,
// density normalization and global iteration offset. A mismatching design is
// rejected; failed calls leave state untouched. Empty state starts a new solve.
struct GpState {
  GpConvergence convergence;
  std::uint64_t signature = 0;
  int iteration = 0;
  double density_weight = 0, momentum = 1, step = 0;
  bool density_initialized = false;
  double previous_overflow = -1, previous_hpwl = 0;
  std::vector<double> x, y, prev_x, prev_y, old_gx, old_gy, old_ex, old_ey;
  std::vector<double> inflation;
  std::vector<AlignmentPair> align_pairs;
};

// Throws invalid_argument for bad inputs and runtime_error if legalization
// cannot produce a legal result. Failed calls leave the caller's design intact.
// Non-movable nodes remain bit-identical. No Torch, Python, CUDA, or file I/O.
Result place(Design &design, const Options &options = {},
             GpState *state = nullptr);
double hpwl(const Design &design);
bool is_legal(const Design &design, std::string *reason = nullptr);

// Reusable per-instance configuration. The caller owns Design (and any HHDS
// mapping); every call creates fresh solver scratch state and updates only the
// design's coordinates on success (also checkpoint-updated geometry/weights
// when a checkpoint is configured). Options are validated when place() runs.
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
std::vector<AlignmentPair> macro_align_pairs(const Design &,
                                             bool nearest = false);
double macro_align_ramp(const Options &, int iteration);
Gradient macro_alignment(const Design &, const std::vector<AlignmentPair> &,
                         double lambda_x, double lambda_y, double sigma,
                         double *residual = nullptr);
AlignmentSnap macro_align_snap(const Design &, double tau, std::int64_t weight);
AlignmentReport macro_alignment_report(const Design &, const Design &,
                                       const AlignmentSnap &);
Gradient cohesion(const Design &design, double lambda);
// Paths/masters are supplied by the caller; identity uses relative macro leaf
// structure, with minimum_copies qualifying siblings under one immediate
// hierarchy parent. Count siblings before collapsing unary wrappers; keep the
// deepest eligible scope for any identical macro set. No connectivity proof
// or template reuse is required for cohesion.
std::vector<CentroidGroup>
centroid_candidates(const Design &, const std::vector<std::string> &paths,
                    const std::vector<std::string> &masters,
                    double minimum_equivalent_cells, int minimum_copies = 4);
Gradient centroid_attraction(const Design &, double lambda);
struct CircleOverlap {
  std::uint64_t pairs = 0, overlaps = 0;
  double maximum = 0;
};
CircleOverlap centroid_circle_overlaps(const Design &, double density);
Gradient centroid_repulsion(const Design &, double lambda, double density);
double centroid_cell_area(const Design &);
std::vector<CentroidSample> centroid_summary(const Design &);
// Squared exact intersection area, normalized by pair size for unit invariance.
Gradient macro_overlap(const Design &, double weight);
Gradient fixed_macro_overlap(const Design &, int bx, int by, double weight);
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
void legalize(Design &d, int bx, int by, const AlignmentSnap *snap = nullptr,
              AlignmentReport *report = nullptr);
void detailed_place(Design &d, int bx, int by, int threads,
                    bool preserve_macros = false, int iterations = 2,
                    std::vector<DetailedSample> *history = nullptr);
bool kernel_is_legal(const Design &d);
} // namespace detail
} // namespace DreamPlace::native
