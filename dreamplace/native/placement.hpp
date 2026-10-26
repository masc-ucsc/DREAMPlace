// BSD-3-Clause. See LICENSE.
#pragma once
#define DREAMPLACE_NATIVE_PIN_TIMING_API 3
#define DREAMPLACE_NATIVE_DENSITY_RESERVATION_API 2

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
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
// Disjoint physical color/region groups. Explicit rigid offsets below permit
// logical member aliases after contraction. Empty areas use width*height.
// Optional areas preserve physical weighting when occupancy is halo padded.
struct CohesionGroup {
  std::vector<int> nodes;
  double weight = 1;
  std::vector<double> areas;
  // Optional logical-member centers relative to the owning rectangle center.
  // Rigid contraction can put several logical members on the same node while
  // preserving their original cohesion energy and force by superposition.
  std::vector<std::pair<double, double>> offsets;
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
  // Optional logical samples on explicit nodes after a rigid contraction.
  // Repeated physical nodes are permitted only with these explicit weights
  // and center offsets, so different hierarchy anchors retain their moments.
  std::vector<double> areas;
  std::vector<std::pair<double, double>> offsets;
};
struct CentroidSample {
  double x = 0, y = 0, radius = 0;
  double area = 0; // physical anchored area; fillers/halo padding excluded
  int cells = 0, macros = 0;
};
enum class PlacementMask : std::uint8_t { none, hard, soft };

// Per-sink timing wirelength objective. Arrays use physical pin indices;
// driver=-1 and weight=0 disables a pin. Empty arrays explicitly clear it.
struct PinGapPull {
  int driver = -1;
  // Fixed world-coordinate gap-entry target, in the design's length units.
  double x = 0, y = 0, fraction = 0;
  // A checkpoint may freeze suppression of conflicting straight pairs for a
  // critical gap. No position-dependent branch changes the objective mid-step.
  bool critical = false;
  std::vector<int> suppressed_sinks{};
  bool operator==(const PinGapPull &) const = default;
};
struct PinTiming {
  std::vector<int> driver;
  std::vector<double> weight;
  double coefficient = .25;
  // Fractions share the original per-driver sum of sink weights. This is a
  // bounded part of the original timing weight, never a physical net/pin.
  std::vector<PinGapPull> gap_pulls{};
  double maximum_gap_fraction = .25;
  bool enabled() const;
  bool operator==(const PinTiming &) const = default;
};
// A frozen checkpoint footprint, never a netlist/pin/obstacle node. A bound
// materialized node contributes through its ordinary physical density exactly
// once; its predicted rectangle remains only for virtual-to-real diagnostics.
struct DensityReservation {
  std::string owner;
  double x = 0, y = 0, width = 0, height = 0;
  int materialized_node = -1;
  bool operator==(const DensityReservation &) const = default;
};
struct ReservationCharge {
  double virtual_area = 0, materialized_area = 0;
  double predicted_materialized_area = 0;
  size_t virtual_count = 0, materialized_count = 0;
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
  // Fixed, pinless placement-only geometry, never a routing obstruction.
  // Hard masks represent missing sites; soft masks guide GP but are ignored
  // by legalization/detail, allowing repair cells in reserved channels.
  std::vector<PlacementMask> node_placement_mask;
  PlacementMask placement_mask(size_t i) const {
    return i < node_placement_mask.size() ? node_placement_mask[i]
                                          : PlacementMask::none;
  }
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
  // Rigid meta-macros cannot be flipped independently of their member pins.
  // Empty preserves the ordinary orientation API.
  std::vector<std::uint8_t> node_orientation_fixed;
  // Diagnostic identities only. Expanded member macros retain the group id;
  // no downstream physical constraint is implied by this field.
  std::vector<std::uint64_t> node_group_id;
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
  std::optional<PinTiming> pin_timing;
  std::vector<DensityReservation> density_reservations;

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
enum class GpCheckpointAction { update, retain, restore };
struct GpUpdate {
  // Permanently retire the movable-update hook before applying this update.
  // Existing temporary holds expire normally unless topology is replaced.
  bool stop_movable_updates = false;
  // Optional nonnegative setup violation from the host STA at this checkpoint.
  std::optional<double> timing_violation;
  std::vector<double> net_weights, node_size_x, node_size_y;
  std::vector<double> pin_offset_x, pin_offset_y;
  std::optional<Design> replacement;
  std::vector<int> previous_node;
  // nullopt preserves it; PinTiming{} clears it. Replacement clears old pin
  // identities and requires a fresh timing update on the replacement graph.
  std::optional<PinTiming> pin_timing;
  // nullopt preserves the complete reservation list; an empty list clears it.
  // Replacement must carry its own freshly rebound list, never stale node IDs.
  std::optional<std::vector<DensityReservation>> density_reservations;
  // Retain the accepted full optimizer/filler/objective state, or restore it.
  // Restore cannot carry other updates; spent iterations never rewind.
  GpCheckpointAction action = GpCheckpointAction::update;
};
using GpCheckpoint = std::function<GpUpdate(const Design &, int)>;
struct GpSample;
struct GpCheckpointSample;
struct TerminalPosition {
  int node = -1;
  double x = 0, y = 0;
};
struct MovablePosition {
  int node = -1;
  double x = 0, y = 0;
  int fix_iterations =
      0; // hold for this many subsequent GP steps, then release
};
class UrouteContext;
enum class RouteModel { off, rudy, rsmt, automatic, uroute };
const char *route_model_name(RouteModel model);
struct RoutePolicy {
  RouteModel model = RouteModel::rudy;
  double switch_overflow = .05, target_overflow = 1e-9;
  int hold = 1, degree = 40, accuracy = 3;
  bool probe = false;
  int probe_iterations = 100;
  RouteModel probe_model = RouteModel::rsmt;
  double uroute_switch_overflow = .01;
  int uroute_hold = 2;
  // Only the final cell stage enables the automatic third rung. FP/macro may
  // still request a one-shot URoute convergence probe.
  bool uroute_ladder = false;
};
// Inverse of route_model_name. Throws invalid_argument for an unknown name.
RouteModel route_model_from_name(std::string_view name);
// The checks Options validation applies to Options::route, for hosts that
// resolve stage policies before placement. Throws invalid_argument.
void validate_route_policy(const RoutePolicy &);
// The policy can evaluate the learned model: an explicit URoute rung, the
// automatic ladder or a URoute convergence probe. A missing context is the
// caller's skip/fallback decision.
bool uroute_requested(const RoutePolicy &);
// The active model can change within the stage (automatic switch or probe
// latch), so stage reports and caches record the model actually used.
bool route_latches(const RoutePolicy &);
// Optimizer-only absolute center targets. No physical nodes/pins/nets are
// added. Targets survive magnet-round restarts within one macro stage, but
// never enter GpState or the final cell-placement handoff.
struct MacroPull {
  int node = -1;
  double x = 0, y = 0, weight = 0;
};
struct MagnetOptions {
  bool enabled = false;
  // Let ordinary GP finish, then merge immediately before the macro snap/LP.
  bool before_snap = false;
  double settle_displacement = .05;
  // One-shot global motion gate before the per-object settling checks.
  double start_displacement = .10;
  double max_gap = 1., density_fraction = .25;
  int settle_window = 10;
  int max_per_round = 0, max_rounds = 0;
  void validate() const;
};
struct MagnetReport;
struct MagnetMacroFrame;
namespace detail {
class MagnetRuntime;
}
struct MacroTransparency {
  // Zero cycles with always=false preserves the original GP exactly. Fixed
  // macro bodies (also host-created macro-shaped blockage masks) remain
  // physically immutable.
  int cycles = 0, warmup = 100, period = 20;
  bool always = false; // diagnostic: transparent for the entire GP solve
  bool enabled() const { return always || cycles > 0; }
};
// completed_steps is zero-based: warmup=100 first disables step 101.
bool macro_transparent(const MacroTransparency &, int completed_steps);
// One-based GP step: hold initial size through start, then grow over duration.
double vanished_memory_scale(int iteration, double initial_scale = .6,
                             int growth_iterations = 25,
                             int start_iteration = 70);
double vanished_memory_local_weight(double initial, int iteration,
                                    double floor = 0, int start_iteration = 70);
// Opt-in, within-run speculative final-GP search. All candidates share the
// same timing weights and complete optimizer/filler checkpoint. Host timing
// is scored serially after parallel native trials; it must not change weights.
struct GpTrialTuning {
  double filler_charge = 1, local_weight = 0;
};
struct GpTrialScore {
  double route_overflow = 0, wns = 0, hpwl = 0, tns = 0;
  double physical_overflow = 0;
  // Separate geometric access failures from capacity overflow. Counts are
  // unweighted active nets/pins; overflow already includes their penalties.
  int route_unroutable_nets = 0, route_blocked_pins = 0,
      route_measured_nets = 0;
  double route_congestion_overflow = 0, route_peak_utilization = 0;
  double route_detour_wirelength = 0;
};
struct GpTrialSample {
  int iteration = 0, steps = 0, round = 0, selected = 0;
  bool adjusting_fillers = false;
  std::string reason;
  GpTrialScore before;
  std::vector<GpTrialTuning> settings;
  std::vector<GpTrialScore> scores;
};
struct GpSearch {
  bool enabled = false;
  // Zero rounds keeps adapting until the final settling window.
  int start = 220, interval = 40, steps = 3, rounds = 0;
  double filler_step = .2, local_step = .1;
  // Lexicographic near-best filtering; comparisons always use the best value
  // of the remaining set, avoiding non-transitive pairwise epsilon ordering.
  double route_tolerance = .01, timing_tolerance = .01, hpwl_tolerance = .005;
};
int select_gp_trial(const std::vector<GpTrialScore> &, const GpSearch &,
                    std::string *reason = nullptr);
struct Options {
  MagnetOptions magnet;
  // Host feasibility check before selecting a rigid macro merge. Members use
  // original input-node identities and physical poses relative to the proposed
  // group origin. An empty result accepts; otherwise report the rejection.
  // Called serially, only after geometry, density and motion checks pass.
  std::function<std::string(const std::vector<MagnetMacroFrame> &)>
      macro_merge_check;
  // Private round continuation: a real node contraction starts a fresh solver
  // at current coordinates while preserving the global event/force schedule.
  std::shared_ptr<detail::MagnetRuntime> magnet_runtime;
  // Read-only per-GP-step macro telemetry, independent of optimization and
  // convergence. Coordinates may include temporary fillers; node IDs are real.
  std::function<void(const Design &, const std::vector<double> &,
                     const std::vector<double> &, int)>
      macro_motion_observer;
  // Physical-member-weighted center-step sum in layout units. Invalid steps
  // rebase across merges/snap, rather than counting their packing jumps.
  std::function<void(int, double, bool)> macro_step_observer;
  int iteration_offset = 0;
  GpSearch gp_search;
  // Internal branch continuation only; no host callback runs in worker threads.
  std::optional<GpTrialTuning> trial_tuning;
  std::function<std::pair<double, double>(const Design &)> trial_timing;
  std::function<void(const GpTrialSample &)> trial_observer;
  MacroTransparency macro_transparency;
  bool vanished_memory =
      false; // final GP: grow fixed exclusion while reducing local pressure
  double vanished_memory_initial_scale = .6; // linear width/height fraction
  int vanished_memory_growth_iterations = 25;
  int vanished_memory_start = 70; // growth/local decay begin on the next step
  double vanished_memory_local_floor =
      0; // retain this local weight after the ramp
  // Safe GP boundary independent of timing. nullopt retains targets; an empty
  // vector clears them. A replacement atomically validates and restarts
  // momentum.
  int macro_pull_interval = 0;
  bool macro_pull_after_density = false;
  int macro_pull_settle_iterations =
      0; // minimum GP steps after an anchor change, within the iteration cap
  std::function<std::optional<std::vector<MacroPull>>(const Design &, int)>
      macro_pull_update;
  // Cheap terminal-only updates, independent of density-gated STA checkpoints.
  // Only fixed zero-area terminals may move; callbacks see physical geometry.
  int terminal_interval = 0;
  std::function<std::vector<TerminalPosition>(const Design &, int)>
      terminal_update;
  // One-shot host relocations at safe GP boundaries, in physical coordinates.
  // Only movable standard cells may change. Fillers are hidden, fixed geometry
  // is immutable, and momentum/selected-iterate history restarts after a move.
  int movable_interval = 0;
  std::function<std::vector<MovablePosition>(const Design &, int)>
      movable_update;
  int checkpoint_interval = 0;
  // No in-loop checkpoint or A/B timing evaluation before this global step.
  // Host final-placement reporting callbacks are independent of this gate.
  int checkpoint_start = 0;
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
  double local_density_target = -1; // -1 inherits the global target density
  int local_density_bins = 512;
  // Opt-in protection from extra routing/local spreading, using the largest
  // incident active net weight w >= 1: mobility = 1-a+a/w. Hosts supply timing
  // weights; this is not an STA or hold-slack model. Keeps full source charge
  // and the unit-area response in the inflated footprint. Macro forces and
  // wirelength attraction remain unchanged.
  double timing_spread_weight = 0; // a in [0,1]; 0 preserves the original flow
  // Extra standard-cell pressure on movable macros; host enables macro GP
  // only. Reach is a screening length in largest physical macro side lengths,
  // independent of short-range cell-to-cell spreading.
  double macro_cell_density_weight = 0;
  double macro_cell_density_reach = 2;
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
  // Start a bounded refinement at the existing final smoothing scale. One
  // preserves the original GP schedule exactly; no public host default changes.
  double wirelength_gamma_start = 1;
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
  // Optional particle count guard. Zero derives count from area and resolution.
  // A positive cap enlarges both dimensions, never the physical design.
  int filler_max_count = 0;
  double filler_area_fraction = 1; // whitespace charge multiplier, (0,2]
  bool filler_cell_space =
      false; // budget: (1-target_density) * non-macro space
  double filler_macro_fraction =
      0; // opt-in: smallest physical macro axes, before halos
  double filler_cell_scale = 1; // central-90%-mean standard-cell dimensions
  double filler_bin_fraction =
      .25; // maximum footprint fraction per density-bin axis
  // Soft centroid cohesion, inverse coordinate units. Lambda ramps linearly
  // from zero to this value over schedule_iterations (or max_iterations).
  // Zero disables it exactly; no boxes or changes to legalization.
  double cohesion_weight = 0;
  double order_weight = 0, order_growth = 1;
  int schedule_iterations = 0;
  RoutePolicy route;
  std::shared_ptr<const UrouteContext> uroute;
  std::string uroute_identity, uroute_skipped;
  std::function<void(const Design &, int bins_x, int bins_y)>
      route_probe_observer;
  // Optional expansion of contracted macro geometry for routing measurement.
  // Inflation still consumes the resulting utilization on the same core grid.
  std::function<Design(const Design &)> route_physical_snapshot;
  std::function<void()> route_initialize;
  bool routability = false;
  double route_capacity_h = 1, route_capacity_v = 1;
  // Experimental filler-free upper-bound penalty. Underfull bins contribute
  // no charge; fixed obstacles alone must not push cells into whitespace.
  bool density_upper_bound = false;
  bool add_fillers = true;
  bool legalize = true;
  // Legalize only hard macros; leave movable standard cells at GP coordinates.
  // Requires legalize=true and detailed_place=false. Result::legal still
  // describes the complete design, including the provisional standard cells.
  bool macro_only_legalize = false;
  bool detailed_place = true;
  int detailed_iterations =
      2; // iteration cap per existing detailed-placement kernel
  // Runs on the selected final GP snapshot, before any legalization.
  std::function<void(Design &)> before_legalize;
  // Refine legal macro origins before standard-cell rows are packed.
  std::function<void(Design &)> after_macro_legalize;
  PlacementObserver observer;
  int observer_interval = 20; // GP steps between snapshots, when enabled
  // Fixed-footprint transactions only; unsupported mutable policies fail.
  bool checkpoint_rollback = false;
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
  double route_rudy_seconds = 0, route_rsmt_seconds = 0,
         route_uroute_seconds = 0;
  int route_uroute_switch_iteration = -1;
  double route_switch_uroute_rsmt = -1, route_switch_uroute = -1;
  double route_probe_overflow = -1;
  int magnet_round = 0, magnet_groups = 0;
  double macro_move_max = 0, macro_move_rms = 0;
  int last_instability_iteration = 0;
  double macro_opacity = 1;
  double macro_exclusion_scale = 1, local_density_weight = 0;
  int filler_count = 0;
  double filler_area = 0, filler_size_x = 0, filler_size_y = 0;
  double pin_timing_cost = 0;
  double reserved_area = 0, materialized_reserved_area = 0;
  double effective_filler_area = 0;
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
struct GpRollbackSample {
  int iteration = 0, restored_iteration = 0;
};
struct GpCheckpointSample {
  int iteration = 0;
  double overflow = 0;
  bool invoked = false, forced = false;
};
struct Result {
  std::shared_ptr<MagnetReport> magnet;
  double magnet_seconds = 0;
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
  int route_uroute_switch_iteration = -1;
  double route_switch_uroute_rsmt = -1, route_switch_uroute = -1;
  double route_probe_overflow = -1;
  std::string uroute_skipped;
  std::vector<std::string> uroute_sample_ids;
  int filler_count = 0;
  double filler_area = 0;
  double filler_size_x = 0, filler_size_y = 0;
  bool filler_count_limited = false;
  double reserved_area = 0, materialized_reserved_area = 0;
  double effective_filler_area = 0, physical_overflow = 0;
  double macro_alignment_seconds = 0;
  double global_placement_seconds = 0, detailed_placement_seconds = 0;
  std::vector<GpSample> history;
  std::vector<DetailedSample> detailed_history;
  int iterations = 0, selected_iteration = 0;
  double local_density_seconds = 0, local_density_force = 0;
  double macro_cell_density_seconds = 0, macro_cell_density_force = 0;
  double macro_cell_density_screening_length = 0;
  double macro_cell_density_coefficient = 0;
  double macro_density_pull_seconds = 0, macro_density_pull_force = 0;
  double macro_density_pull_anchor_ratio = 0;
  double local_density_grid_seconds = 0;
  int timing_spread_cells = 0; // last weight update, eligible cells with m < 1
  double timing_spread_min_mobility = 1;
  double initial_hpwl = 0, hpwl = 0, overflow = 0;
  bool converged = false;
  bool legal = false;
  std::vector<GpRollbackSample> rollback_history;
};

// Optional full GP continuation. Includes private fillers, Nesterov/BB history,
// density normalization and global iteration offset. A mismatching design is
// rejected; failed calls leave state untouched. Empty state starts a new solve.
struct GpState {
  GpConvergence convergence;
  RouteModel route_model = RouteModel::rudy;
  int route_uroute_hold = 0;
  int route_hold = 0, route_switch_iteration = -1, route_probe_deadline = 0;
  double route_switch_rudy = -1, route_switch_rsmt = -1;
  bool route_probed = false;
  int route_uroute_switch_iteration = -1;
  double route_switch_uroute_rsmt = -1, route_switch_uroute = -1;
  double route_probe_overflow = -1;
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
  std::optional<PinTiming> pin_timing;
  std::vector<DensityReservation> density_reservations;
  bool reservation_accounting = false;
};

// Validates identities/physical bindings and returns nominal charge. Inflation,
// when supplied, prices bound nodes at their actual current density charge.
ReservationCharge reservation_charge(const Design &,
                                     const std::vector<double> &inflation = {});
namespace detail {
// Deposit only virtual reservations through the ordinary native density kernel.
// No capacity is subtracted and no virtual node enters maps/netlists or RUDY.
void reservation_density_map(const Design &, int bins_x, int bins_y,
                             bool exact_cells, int threads,
                             std::vector<double> &map);
// Full reservation/physical charge displaces equal filler charge, independently
// of the user's filler fraction. Insufficient filler budget is explicit.
double reservation_filler_area(double unreserved_budget, double rounding_loss,
                               const ReservationCharge &);
} // namespace detail

// Throws invalid_argument for bad inputs and runtime_error if legalization
// cannot produce a legal result. Failed calls leave the caller's design intact.
// Non-movable nodes remain bit-identical. No Torch, Python, CUDA, or file I/O.
Result place(Design &design, const Options &options = {},
             GpState *state = nullptr);
double hpwl(const Design &design);
// Exact physical overflow excluding fillers; fixed area reduces bin capacity.
double density_overflow(const Design &, int bins_x, int bins_y,
                        double target_density,
                        bool scale_macro_density = false);
// Solver-compatible overflow with virtual reservations counted as movable
// charge once. Physical-only density_overflow above retains its old meaning.
double reservation_density_overflow(const Design &, int bins_x, int bins_y,
                                    double target_density,
                                    bool scale_macro_density = false);
bool is_legal(const Design &design, std::string *reason = nullptr);
bool macros_are_legal(const Design &design, std::string *reason = nullptr);

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
// Cell of p on count uniform cells of width step starting at lo. Clamps
// before the int conversion, so far-outside points land in an edge cell.
inline int clamped_bin(double p, double lo, double step, int count) {
  return static_cast<int>(
      std::clamp(std::floor((p - lo) / step), 0., double(count - 1)));
}
// Sum coefficient*w[p]*(dx*tanh(dx/(2*gamma))+dy*tanh(dy/(2*gamma)))
// between each weighted sink and its declared driver: exactly the two-pin
// weighted-average wirelength smoothing. Gamma uses layout units;
// weights/coefficient are dimensionless. Gradients include both endpoint nodes.
// No density charge or physical topology is added.
void validate_pin_timing(const Design &);
Gradient pin_timing_force(const Design &, double gamma);
// Raw gradients compensate the GP preconditioner so equal target offsets
// produce equal pull-induced center steps across physical macro sizes.
struct MacroPullForce {
  int node;
  double x, y;
};
std::vector<MacroPullForce> macro_pull_forces(const Design &,
                                              const std::vector<MacroPull> &,
                                              const Gradient &wire,
                                              const Gradient &density,
                                              double density_weight);
std::vector<AlignmentPair>
macro_align_pairs(const Design &, bool nearest = false,
                  const std::vector<int> *classes = nullptr);
std::vector<AlignmentPair>
macro_grid_align_pairs(const Design &, bool nearest = false,
                       const std::vector<int> *classes = nullptr,
                       bool padded_extents = false,
                       bool overlapping_pairs = false);
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
Gradient fixed_macro_overlap(const Design &, int bx, int by, double weight,
                             double exclusion_scale = 1);
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
              const std::function<void(Design &)> &after_macro_legalize = {},
              bool macro_only = false);
void detailed_place(Design &d, int bx, int by, int threads,
                    bool preserve_macros = false, int iterations = 2,
                    std::vector<DetailedSample> *history = nullptr);
bool kernel_is_legal(const Design &d);
} // namespace detail
} // namespace DreamPlace::native
