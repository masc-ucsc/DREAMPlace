// BSD-3-Clause. See LICENSE.
#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "dreamplace/native/constraints.hpp"
#include "dreamplace/native/floorplan.hpp"
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
// Candidate colors may overlap across groups. Explicit nodes may carry a
// color, but cannot also be covered by a color listed in that same group.
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
  // Negative uses member area. Otherwise this is an independent space budget;
  // it changes circle radii only, never centroid moments or density charge.
  double repulsion_area = -1;
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
  // Host alignment identity; -1 opts a node out. LiveHD uses master x color
  // normally, and master alone for global strip-grid alignment.
  std::vector<int> node_align_class;
  // Same original hierarchy scope AND macro dimensions; -1 excludes a node.
  // Used by both macro snap algorithms, never by GP alignment/cohesion.
  std::vector<int> node_snap_class;
  // Legacy LEF blocked-layer metadata. Enhanced RUDY treats every macro
  // body as opaque, independently of these values; halos remain routing space.
  std::vector<double> node_blocked_capacity_h, node_blocked_capacity_v;
  int num_movable = 0;
  std::vector<double> node_x, node_y, node_size_x, node_size_y;
  // Optional timing displacement costs; empty preserves the legacy solver.
  std::vector<double> node_legalize_weights;
  // Defer characterized constant sources until functional users own their
  // sites.
  std::vector<std::uint8_t> node_legalize_last;
  std::vector<double> node_legalize_anchor_x, node_legalize_anchor_y;
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
struct TerminalPosition {
  int node = -1;
  double x = 0, y = 0;
};
enum class RouteModel { off, rudy, rsmt, automatic };
const char *route_model_name(RouteModel model);
struct RoutePolicy {
  RouteModel model = RouteModel::rudy;
  double switch_overflow = .05, target_overflow = 1e-9;
  int hold = 1, degree = 40, accuracy = 3;
  bool probe = false;
  int probe_iterations = 100;
};
// Optimizer-only absolute center targets. No physical nodes/pins/nets are
// added.
struct MacroPull {
  int node = -1;
  double x = 0, y = 0, weight = 0;
};
struct Options {
  // Safe GP boundary independent of timing. nullopt retains targets; an empty
  // vector clears them. A replacement atomically validates and restarts
  // momentum.
  int macro_pull_interval = 0;
  std::function<std::optional<std::vector<MacroPull>>(const Design &, int)>
      macro_pull_update;
  // Cheap terminal-only updates, independent of density-gated STA checkpoints.
  // Only fixed zero-area terminals may move; callbacks see physical geometry.
  int terminal_interval = 0;
  std::function<std::vector<TerminalPosition>(const Design &, int)>
      terminal_update;
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
  // OpenROAD's movable-area/target-density rule, with a configurable resource
  // cap. Unlike legacy automatic_bins (macro pitch), this resolves before halo
  // padding.
  bool dynamic_bins = false;
  int bin_max = 512;
  // Independent auxiliary grid; 0 follows density bins. Keeps RUDY/LG/DP
  // resolution unchanged.
  int auxiliary_bins = 0;
  double local_density_weight = 0, local_density_range = 4,
         local_density_jitter = .01;
  int local_density_bins = 512;
  // Opt-in protection from extra routing/local spreading, using the largest
  // incident active net weight w >= 1: mobility = 1-a+a/w. Hosts supply timing
  // weights; this is not an STA or hold-slack model. Keeps full source charge
  // and the unit-area response in the inflated footprint. Macro forces and
  // wirelength attraction remain unchanged.
  double timing_spread_weight = 0; // a in [0,1]; 0 preserves the original flow
  // Extra local standard-cell pressure on movable macros. The host enables
  // it only during FP/macro GP; fixed macros never receive force.
  double macro_cell_density_weight = 0;
  // Independent weak, local cell-density slope on macro receivers. Weight
  // is relative to ordinary macro force and capped by centroid reference
  // force when anchored. Range is in local-density bin pitches.
  double macro_density_pull_weight = 0, macro_density_pull_range = 4;
  bool exact_cell_density =
      false; // opt-in physical footprint for movable cells
  bool adaptive_density = false, improved_routability = false;
  bool scale_macro_density = false;
  // Direct pressure on facing macro edges from the RUDY channel utilization.
  // Native opt-in; the LiveHD routability flow enables it. Zero is
  // bit-identical.
  double routability_macro_weight = 0,
         routability_macro_range = 2; // range in bins
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
  // Rebuild physical-projection pairs every iteration using minimum-size
  // strips. Host classes must identify master types, independently of color.
  bool macro_align_global = false;
  // Opt-in centroid-style line/gap/midpoint fitting before macro legalization.
  bool macro_snap_regular = false;
  FloorplanOptions macro_snap_options;
  double macro_snap_tolerance = .25;
  std::int64_t macro_snap_weight = 100;
  // Zero disables R8H2 exactly. Attraction starts in the first GP iteration;
  // no preliminary macro placement or extra continuation phase is used.
  double centroid_weight = 0;
  // Calibrate attraction once at one circle radius, retaining restoring force
  // as members drift. The legacy per-iteration normalization stays available.
  bool centroid_distance_attraction = false;
  // Frozen targets only: keep all-member area/calibration, omit cell springs.
  bool centroid_anchor_macros_only = false;
  // Additional stiffness outside the area/density circle; 4 means 5x slope.
  double centroid_attraction_outside = 4;
  // Relative to centroid attraction. Same-family circles have area A/density.
  // Zero preserves the attraction-only path exactly.
  double centroid_repulsion_weight = 0;
  // Increase the circle penalty until separation and density both settle.
  // Off by default; relative penetration tolerance is measured against r_i+r_j.
  bool centroid_adaptive_repulsion = false;
  double centroid_stop_overlap = .001;
  int centroid_ramp = 100;
  // Moving-centroid alignment during FP only. Same family and candidate parent.
  double centroid_align_weight = 0, centroid_align_sigma = .25;
  int centroid_align_start = 100, centroid_align_full = 500;
  // FP-only child-to-parent spring, relative to calibrated member attraction.
  double centroid_parent_weight = 0, centroid_parent_outside = 4;
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
  // Optional virtual-particle budget. Zero retains the legacy size/count rule.
  // A positive cap enlarges both dimensions, never the physical design.
  int filler_max_count = 0;
  double filler_area_fraction = 1; // fraction of whitespace density charge
  // Soft centroid cohesion, inverse coordinate units. Lambda ramps linearly
  // from zero to this value over schedule_iterations (or max_iterations).
  // Zero disables it exactly; no boxes or changes to legalization.
  double cohesion_weight = 0;
  double order_weight = 0, order_growth = 1;
  int schedule_iterations = 0;
  RoutePolicy route;
  std::function<void()> route_initialize;
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
  // Runs on the selected final GP snapshot, before any legalization.
  std::function<void(Design &)> before_legalize;
  // Refine legal macro origins before standard-cell rows are packed.
  std::function<void(Design &)> after_macro_legalize;
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
struct MacroSnapTarget {
  int node;
  double x, y;
};
struct AlignmentSnap {
  std::vector<MacroSnapTarget> targets;
  FloorplanSnap regularity;
  std::vector<AlignmentGroup> groups;
  std::vector<::DreamPlace::MacroAlignConstraint> constraints;
};
struct AlignmentReport {
  FloorplanSnap regularity;
  std::vector<::DreamPlace::MacroAlignConstraint> constraints;
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
  // Maximum applied L1 force / channel count since the preceding sample.
  double route_peak_utilization = 0, macro_routing_force = 0;
  std::uint64_t macro_routing_channels = 0;
  double align_force = 0; // peak applied L1 force since the preceding sample
  std::uint64_t centroid_align_pairs = 0;
  double centroid_align_residual = 0, centroid_align_ramp = 0;
  double centroid_align_force = 0; // peak applied member L1 force
  std::uint64_t centroid_parent_edges = 0;
  double centroid_parent_max_distance = 0, centroid_parent_force = 0;
  RouteModel route_model = RouteModel::rudy;
  double route_seconds = 0, route_switch_rudy = -1, route_switch_rsmt = -1;
  int route_switch_iteration = -1, route_bins = 0;
  bool route_sampled = false, route_probe = false;
  int route_blocked_pins = 0, route_unroutable_nets = 0, route_calls = 0;
  double route_rudy_seconds = 0, route_rsmt_seconds = 0;
};
struct CentroidReport {
  bool started = false, completed = false;
  int start_iteration = 0, iterations = 0;
  double force_seconds = 0, alignment_seconds = 0, parent_seconds = 0;
  double attraction_coefficient = 0;
  int attraction_calibration_iteration = -1;
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
  double poisson_seconds = 0, macro_routing_seconds = 0;
  double rudy_seconds = 0;
  int rudy_calls = 0;
  RouteModel route_model = RouteModel::rudy;
  int route_switch_iteration = -1;
  double route_switch_rudy = -1, route_switch_rsmt = -1;
  bool route_probed = false;
  int filler_count = 0;
  double filler_area = 0;
  double macro_alignment_seconds = 0;
  double global_placement_seconds = 0, detailed_placement_seconds = 0;
  std::vector<GpSample> history;
  std::vector<DetailedSample> detailed_history;
  int iterations = 0, selected_iteration = 0;
  double local_density_seconds = 0, local_density_force = 0;
  double macro_cell_density_seconds = 0, macro_cell_density_force = 0;
  double macro_density_pull_seconds = 0, macro_density_pull_force = 0;
  double macro_density_pull_anchor_ratio = 0;
  double local_density_grid_seconds = 0;
  int timing_spread_cells = 0; // last weight update, eligible cells with m < 1
  double timing_spread_min_mobility = 1;
  double initial_hpwl = 0, hpwl = 0, overflow = 0;
  bool converged = false;
  bool legal = false;
};

// Optional full GP continuation. Includes private fillers, Nesterov/BB history,
// density normalization and global iteration offset. A mismatching design is
// rejected; failed calls leave state untouched. Empty state starts a new solve.
struct GpState {
  GpConvergence convergence;
  RouteModel route_model = RouteModel::rudy;
  int route_hold = 0, route_switch_iteration = -1, route_probe_deadline = 0;
  double route_switch_rudy = -1, route_switch_rsmt = -1;
  bool route_probed = false;
  std::uint64_t signature = 0;
  int iteration = 0;
  double density_weight = 0, momentum = 1, step = 0;
  bool density_initialized = false;
  double previous_overflow = -1, previous_hpwl = 0;
  std::vector<double> x, y, prev_x, prev_y, old_gx, old_gy, old_ex, old_ey;
  std::vector<double> inflation;
  std::vector<double>
      routing_utilization; // last RUDY sample, for macro pressure
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
std::vector<AlignmentPair>
macro_align_pairs(const Design &, bool nearest = false,
                  const std::vector<int> *classes = nullptr);
std::vector<AlignmentPair> macro_grid_align_pairs(const Design &,
                                                  bool nearest = false);
double macro_align_ramp(const Options &, int iteration);
Gradient macro_alignment(const Design &, const std::vector<AlignmentPair> &,
                         double lambda_x, double lambda_y, double sigma,
                         double *residual = nullptr,
                         bool physical_extents = false);
AlignmentSnap macro_align_snap(const Design &, double tau, std::int64_t weight);
AlignmentSnap macro_regular_snap(const Design &, const FloorplanOptions &,
                                 std::int64_t weight);
AlignmentReport macro_alignment_report(const Design &, const Design &,
                                       const AlignmentSnap &);
Gradient cohesion(const Design &design, double lambda);
// Paths/masters are supplied by the caller; identity uses relative macro leaf
// structure, with minimum_copies qualifying siblings under one immediate
// hierarchy parent. Count siblings before collapsing unary wrappers; keep the
// deepest eligible scope for any identical macro set. A scope needs more than
// minimum_macros movable macros. No connectivity proof
// or template reuse is required for cohesion.
std::vector<CentroidGroup>
centroid_candidates(const Design &, const std::vector<std::string> &paths,
                    const std::vector<std::string> &masters,
                    double minimum_equivalent_cells, int minimum_copies = 4,
                    int minimum_macros = 2);
Gradient centroid_attraction(const Design &, double lambda);
Gradient centroid_parent_attraction(const Design &, double lambda,
                                    double density, double outside,
                                    std::uint64_t *edges = nullptr,
                                    double *max_distance = nullptr);
Gradient centroid_circle_attraction(const Design &, double lambda,
                                    double density, double outside,
                                    bool macros_only = false);
// Per-node magnitude of the direct attraction gradient at one circle radius.
// Add contributions when a member participates in multiple candidates.
std::vector<double> centroid_attraction_reference(const Design &,
                                                  double density);
struct CircleOverlap {
  std::uint64_t pairs = 0, overlaps = 0;
  double maximum = 0;
};
CircleOverlap centroid_circle_overlaps(const Design &, double density);
Gradient centroid_repulsion(const Design &, double lambda, double density);
// Reference force is the net wirelength+density force of each candidate's
// movable members. Alignment uses the same bounded well as macro alignment;
// gradients return through the area-weighted moving-centroid chain rule.
std::pair<double, double> centroid_force_reference(const Design &,
                                                   const Gradient &,
                                                   const Gradient &, double);
Gradient centroid_alignment(const Design &, double lambda_x, double lambda_y,
                            double sigma, double density,
                            std::uint64_t *pairs = nullptr,
                            double *residual = nullptr);
double centroid_cell_area(const Design &);
std::vector<CentroidSample> centroid_summary(const Design &);
// Squared exact intersection area, normalized by pair size for unit invariance.
Gradient macro_overlap(const Design &, double weight);
Gradient fixed_macro_overlap(const Design &, int bx, int by, double weight);
std::pair<std::uint64_t, std::uint64_t> macro_overlap_counts(const Design &);
std::pair<int, int> automatic_bins(const Design &d);
std::pair<int, int> dynamic_bins(const Design &d, double target, int maximum);
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
              AlignmentReport *report = nullptr,
              const std::function<void(Design &)> &after_macro_legalize = {});
void detailed_place(Design &d, int bx, int by, int threads,
                    bool preserve_macros = false, int iterations = 2,
                    std::vector<DetailedSample> *history = nullptr);
bool kernel_is_legal(const Design &d);
} // namespace detail
} // namespace DreamPlace::native
