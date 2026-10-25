// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/rigid.hpp"
#include <array>
#include <deque>

namespace DreamPlace::native {
enum class MagnetEvent { timing, pull, terminal, snap, round, ready };
const char *magnet_event_name(MagnetEvent);
struct MacroMovement {
  double maximum = 0, rms = 0;
  bool window_complete = false;
};
// No position copies for cells, fillers or fixed macros. The metric is the
// largest excursion of each current center from any frame in the window,
// rather than endpoint-only motion that could hide an oscillation.
class MagnetStability {
  struct Frame {
    int iteration;
    std::vector<std::pair<double, double>> centers;
  };
  MagnetOptions options_;
  std::deque<Frame> frames_;
  int last_event_ = 0;
  bool armed_ = false;
  MacroMovement movement_;
  std::vector<int> nodes_;
  std::vector<double> displacements_;

public:
  explicit MagnetStability(MagnetOptions);
  void event(int iteration, MagnetEvent);
  MacroMovement observe(const Design &, const std::vector<double> &x,
                        const std::vector<double> &y, int iteration);
  bool gate(int iteration);
  bool settled(int iteration) const;
  double macro_movement(int node) const;
  bool macro_settled(int node) const;
  int last_event() const { return last_event_; }
  bool armed() const { return armed_; }
  MacroMovement movement() const { return movement_; }
};
struct MagnetBox {
  double xl = 0, yl = 0, xh = 0, yh = 0;
};
// Convex quadrilateral joining the two facing padded edges. With equal shapes
// this is a parallelogram; aligned pairs reduce to the ordinary gap rectangle.
using MagnetRegion = std::array<std::pair<double, double>, 4>;
struct MagnetCandidate {
  int a = -1, b = -1;
  Axis axis = Axis::x; // merge axis, not the orthogonal alignment axis
  double gap = 0, gap_density = 0, displacement = 0;
  double center_shift_x = 0, center_shift_y = 0;
  double move_a = 0, move_b = 0;
  MagnetBox gap_box, merged_box;
  MagnetRegion gap_region;
  MacroMerge merge;
  bool selected = false;
  std::string reason;
  std::vector<int> source_a, source_b;
};
// Physical movable standard-cell footprint density; excludes macros, fillers,
// fixed cells and halo charge. Rectangle queries use the summed-area table;
// offset facing-edge regions integrate density weighted by clipped bin area.
class MagnetDensity {
  Layout layout_;
  int bx_, by_;
  double dx_, dy_, global_cell_density_ = 0;
  std::vector<double> rho_, prefix_;
  double integral(double x, double y) const;

public:
  MagnetDensity(const Design &, int bins_x, int bins_y);
  double global_cell_density() const { return global_cell_density_; }
  double density(const MagnetBox &) const;
  double region_density(const MagnetRegion &) const;
  const std::vector<double> &grid() const { return rho_; }
  int bins_x() const { return bx_; }
  int bins_y() const { return by_; }
};
struct MagnetSearch {
  std::vector<MacroMerge> merges;
  std::vector<MagnetCandidate> candidates;
  double seconds = 0;
  double global_cell_density = 0, density_limit = 0;
  bool exact_density = false;
};
MagnetSearch magnet_candidates(
    const Design &, const MagnetOptions &, const MagnetDensity &,
    const MagnetStability * = nullptr,
    const std::function<std::string(const MagnetCandidate &)> & = {});
struct MagnetEventSample {
  int iteration = 0;
  MagnetEvent kind = MagnetEvent::round;
};
struct MagnetRound {
  int iteration = 0;
  double overflow = 0;
  MacroMovement movement;
  MagnetSearch search;
  std::vector<std::uint64_t> group_ids;
  std::string reason;
};
struct MagnetMacroFrame {
  int source = -1;
  double x = 0, y = 0, width = 0, height = 0;
  Orientation orientation = Orientation::N;
  std::uint64_t group = 0;
  bool fixed = false;
};
struct MagnetFrame {
  int iteration = 0;
  std::string event;
  std::vector<MagnetMacroFrame> macros;
  int bins_x = 0, bins_y = 0;
  std::vector<double> density;
};
struct MagnetReport {
  // Every iteration, including incomplete windows and the pre-gate interval.
  struct MotionSample {
    int iteration = 0;
    MacroMovement movement;
    double center_step_sum = 0;
    int physical_macros = 0, gp_macros = 0;
    bool center_step_valid = false;
  };
  std::vector<MotionSample> motion_history;
  std::vector<std::pair<int, double>> overflow_history;
  int start_iteration = -1;
  double start_observed_displacement = -1;
  std::vector<MagnetRound> rounds;
  std::vector<MagnetEventSample> events;
  RigidMacroMap mapping;
  Layout layout;
  double target_density = 0;
  MagnetOptions options;
  std::vector<MagnetFrame> frames;
  bool stopped_settled = false;
  double seconds = 0;
};
namespace detail {
// Undo halo-only geometry changes using the physical reference for constraints
// and logical cohesion samples. Shared by host callbacks and merge preflight.
Design unpad_snapshot(const Design &, const Design &physical_reference);
struct PreparedMagnetMerge {
  Design design;
  RigidMacroMap mapping;
  MacroMergeResult result;
};
// One lifetime across solver restarts. Only macro centers are sampled every
// iteration; the cell-density grid is built at an eligible round boundary.
class MagnetRuntime {
  Options options_;
  MagnetStability stability_;
  double previous_hpwl_ = -1;
  int stable_hpwl_ = 0;
  bool started_ = false;
  int motion_iteration_ = -1;
  std::vector<int> motion_nodes_, motion_weights_;
  std::vector<std::pair<double, double>> motion_centers_;

public:
  explicit MagnetRuntime(const Options &);
  std::shared_ptr<MagnetReport> report = std::make_shared<MagnetReport>();
  std::vector<MacroMerge> pending;
  std::function<std::string(const MagnetCandidate &)> candidate_check;
  std::function<Design(const Design &)> prepare_snapshot;
  std::optional<PreparedMagnetMerge> prepared;
  bool finalizing = false;
  int groups = 0, round = 0;
  double filler_width = 0, filler_height = 0;
  bool density_initialized = false;
  // Stage-wide routing decisions survive a geometry contraction. The hold
  // counter and geometry-dependent utilization maps intentionally restart.
  struct RouteProgress {
    RouteModel model = RouteModel::rudy;
    int probe_deadline = 0, switch_iteration = -1;
    double switch_rudy = -1, switch_rsmt = -1;
    bool probed = false;
    int uroute_switch_iteration = -1;
    double switch_uroute_rsmt = -1, switch_uroute = -1;
    double probe_overflow = -1;
  };
  std::optional<RouteProgress> routing;

  // Timing observations and the original density-continuation budget belong
  // to the macro stage. Geometry-dependent objective samples restart at each
  // contraction, but neither a merge nor snap grants a fresh budget or timing.
  std::optional<GpConvergence> convergence;
  bool checkpoint_invoked = false;

  // Solver-only force state. A merge drops its members' targets and remaps
  // untouched nodes; nullopt at a later refresh must still retain these.
  std::vector<MacroPull> macro_pulls;
  int macro_pull_changed = -1;
  void remap_pulls(const std::vector<int> &previous_node);
  void event(int iteration, MagnetEvent);
  void observe(const Design &, const std::vector<double> &,
               const std::vector<double> &, int iteration);
  void sample(const Design &, GpSample &);
  bool ready(int iteration) const;
  bool settled(int iteration) const { return stability_.settled(iteration); }
};
Result magnet_place(Design &, const Options &);
} // namespace detail
} // namespace DreamPlace::native
