// BSD-3-Clause. See LICENSE.
#pragma once

#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace DreamPlace::native {

struct Layout {
  double xl = 0, yl = 0, xh = 0, yh = 0;
  double site_width = 1, row_height = 1;
};

// Owning flat-array boundary. Coordinates/sizes/offsets use one caller-chosen
// unit, normally DBU. Order: movable, fixed, zero-area IO terminals. Fillers
// are internal to global placement and never exposed in the result.
struct Design {
  Layout layout;
  int num_movable = 0;
  std::vector<double> node_x, node_y, node_size_x, node_size_y;
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

enum class PlacementStage { initial, global_placement, legalization, detailed_placement, finished };
const char* stage_name(PlacementStage stage);
// Synchronous, read-only physical snapshots. References are valid only during
// the call; temporary fillers are excluded. Exceptions abort without writeback.
using PlacementObserver = std::function<void(const Design&, PlacementStage, int iteration)>;

struct Options {
  int max_iterations = 500;
  int bins_x = 32, bins_y = 32;
  int num_threads = 1;
  double target_density = 0.8;
  double stop_overflow = 0.1;
  double initial_density_weight = 0.01;
  double max_step = 0.5;  // fraction of a bin
  bool add_fillers = true;
  bool legalize = true;
  bool detailed_place = true;
  PlacementObserver observer;
  int observer_interval = 20;  // GP steps between snapshots, when enabled
};

struct Result {
  int iterations = 0;
  double initial_hpwl = 0, hpwl = 0, overflow = 0;
  bool converged = false;
  bool legal = false;
};

// Throws invalid_argument for bad inputs and runtime_error if legalization
// cannot produce a legal result. Failed calls leave the caller's design intact.
// Non-movable nodes remain bit-identical. No Torch, Python, CUDA, or file I/O.
Result place(Design& design, const Options& options = {});
double hpwl(const Design& design);
bool is_legal(const Design& design);

// Reusable per-instance configuration. The caller owns Design (and any HHDS
// mapping); every call creates fresh solver scratch state and updates only the
// design's coordinates on success. Options are validated when place() runs.
// Independent instances/designs may run concurrently. Configure options before
// launching work; callbacks and shared host/HHDS state require caller ownership
// or synchronization. num_threads controls workers within each call.
class Placer {
 public:
  explicit Placer(Options options = {}) : options_(std::move(options)) {}

  Options& options() noexcept { return options_; }
  const Options& options() const noexcept { return options_; }
  Result place(Design& design) const { return ::DreamPlace::native::place(design, options_); }

 private:
  Options options_;
};

namespace detail {
struct Gradient {
  double value = 0;
  std::vector<double> x, y;
};
Gradient wirelength(const Design& design, double gamma, int threads);
void density_map(const Design& d, int bx, int by, int threads, const std::vector<double>& sx,
                 const std::vector<double>& sy, const std::vector<double>& ox, const std::vector<double>& oy,
                 const std::vector<double>& ratio, std::vector<double>& map);
Gradient density_force(const Design& d, int bx, int by, int threads, const std::vector<double>& sx,
                       const std::vector<double>& sy, const std::vector<double>& ox, const std::vector<double>& oy,
                       const std::vector<double>& ratio, const std::vector<double>& fx, const std::vector<double>& fy);
void legalize(Design& d, int bx, int by);
void detailed_place(Design& d, int bx, int by, int threads);
bool kernel_is_legal(const Design& d);
}  // namespace detail
}  // namespace DreamPlace::native
