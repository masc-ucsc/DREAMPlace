# Native CPU placement for LiveHD

The public targets are `//:dreamplace` (flat-array C++ placement) and `//:hhds`
(the optional HHDS adapter). Neither target uses the Python placement loop,
Torch tensors, pybind11, CUDA, Boost, Limbo, Cairo, OpenTimer, or PlaceDB.
The upstream Python/CMake workflow is separate and retains its old dependencies.
No submodule initialization or pip installation is needed for the Bazel build.

## Build and run

Use Bazel 9.3.0 (selected by `.bazelversion`) and a C++23 toolchain.

```sh
bazel test -c opt //...
bazel run -c opt //dreamplace/native:example

# Work against the sibling graphs without editing either repository.
bazel test -c opt //... \
  --override_repository=hhds=../hhds \
  --override_repository=hlop=../hlop

# CPU parallelism, using the registry OpenMP runtime, including on macOS.
bazel test -c opt --config=openmp //...

# Raw-pointer boundary and solver checks with assertions and sanitizers.
bazel test --config=asan --config=ubsan //dreamplace/native:all
```

The default build runs the same kernels serially. OpenMP compiles the runtime
from the registry; it needs Python **at build time** to generate the runtime's
message tables. Placement itself never starts or embeds Python. Serial native
build actions do not use Python. On macOS the standalone build explicitly
selects Command Line Tools instead of discovering an unrelated Xcode install.

## Dependency boundary

| Dependency | Source | Used by |
| --- | --- | --- |
| rules_cc, platforms | Bazel Central Registry | C++ build rules |
| apple_support | Bazel Central Registry, development dependency | Standalone macOS toolchain |
| coin-or-lemon 1.3.1.bcr.1 | Bazel Central Registry | Macro difference-constraint LP, network simplex |
| openmp 21.1.5.bcr.3 | Bazel Central Registry | Optional CPU threading |
| HHDS, HLOP, iassert, emhash | Pinned masc-ucsc/upstream revisions in `MODULE.bazel` | HHDS adapter only |
| RapidJSON | Bazel Central Registry | Optional H2D report/cache target only |
| Abseil | Bazel Central Registry | HHDS/HLOP public API |
| FLUTE | Already vendored, separate `//thirdparty/flute:flute` target | Future parasitic estimator |

`packages/lemon_graph.patch` adds a public `:graph` target to registered LEMON. It
uses the registry's patched headers and `base.cc`, retaining its C++ allocator
compatibility fixes. The default LEMON target also builds the GLPK frontend;
this target avoids compiling or linking GLPK. GLPK still appears in LEMON's
module metadata, but is outside the native placement action/link graph.
The direct macro solver expresses the L1 displacement LP as a min-cost
circulation, so neither Limbo nor a general LP solver is required.

Boost is used by the legacy database/I/O code: Geometry/R-tree, Polygon,
Bind, and small string/sort helpers. None of those files is compiled by these
native targets. The placement kernels and HHDS adapter need no Boost subset.

## Flat arrays

Include `dreamplace/native/placement.hpp` and construct a `Design`:

- Layout: rectangular core bounds, site width and row height. The core contains
  whole sites/rows; arbitrary row gaps, fence regions and orientations are not
  represented in this initial boundary.
- Nodes: `node_x/y`, `node_size_x/y`, ordered movable then fixed then zero-area
  IO terminals. FIXED and LOCKED both belong in the fixed suffix.
- Pins: `pin2node`, `pin_offset_x/y` and net-to-pin CSR
  (`flat_net2pin`, `flat_net2pin_start`). One pin belongs to one net.
- Nets: optional weights and masks. `build_maps()` fills unit weights/default
  masks, `pin2net`, and node-to-pin CSR. Empty/singleton nets are disabled by
  default. Explicitly enabled nets must have at least two pins.

All coordinates, dimensions and pin offsets use the same unit. Use integer
DBU for mixed-size placement: macro LP coefficients and displacement anchors
are rounded to integer units, as in upstream. Pin offsets must already reflect
cell orientation; the placer does not rotate cells. Movable dimensions need
not be integral site/row multiples. Macro legalization uses private conservative
site/row-rounded occupancy envelopes; physical dimensions and pin offsets remain
unchanged. Movable origins are still aligned to sites/rows. Zero-area terminals have no density
or overlap contribution, and their positions are preserved exactly.

Optional `node_halo_x/y` arrays reserve clearance per side in the same layout
units. An empty array means zero clearance; otherwise it must contain one finite,
nonnegative value per node. Zero-area terminals cannot have positive halos.
Clearance rounds up to site/row pitches. Placement uses a private padded design
for density, legalization and detail, shifting pin offsets to preserve physical
pin locations. Physical sizes and pin offsets remain unchanged in the caller's
design and observer snapshots. Fixed coordinates remain exact. `is_legal()` checks
padded occupancy when halos are present, and HPWL uses physical pins. The host
chooses which nodes receive halos; LiveHD applies its global option only to
non-CORE LEF macros. If heuristic macro separation constraints are infeasible,
the native LP leaves both axes unchanged and allows the existing Hannan-grid
fallback to run; infeasibility of that graph is not proof the core cannot fit.

```cpp
DreamPlace::native::Design design;
// Fill layout, node/pin arrays, and net CSR from the mapped physical netlist.
design.build_maps();
DreamPlace::native::Placer placer;
placer.options().num_threads = 4; // requires --config=openmp for actual parallelism
placer.options().max_iterations = 1500;
placer.options().target_density = 0.8;
placer.options().initial_density_weight = 0.01;
const auto result = placer.place(design);
// design.node_x/y now hold the placed coordinates.
```

The H2D foundations are available without changing the existing placement flow:

- `Design::node_orientation` optionally records N/FN/FS/S. Empty means N, and
  `pin_offset_x/y` always hold offsets in the current orientation. Explicit
  `set_orientation()` transforms from the old orientation, preserves dimensions
  and origins, and leaves invalid inputs untouched. `place()` does not select
  orientations. Callers supply LEF symmetry restrictions when choosing flips.
- `Options::seed` controls filler initialization; its default remains the
  historical 1000. Changing it does not change fixed coordinates.
- `macro/check.hpp` checks an explicit macro subset, including fixed-fixed
  overlap, site/row pitch and halo-expanded core bounds. This complements the
  existing general placement checker; it does not silently move illegal objects.
- `rudy.hpp` measures directional routing utilization and mean bin overflow
  without Torch, graph dependencies, placement edits or cell inflation. Track
  capacities use tracks per layout unit. It uses the existing WDM degree weights,
  physical pins and net masks. Zero-area net bboxes contribute zero under the
  upstream projected-box formula and are counted separately, so this estimate
  does not capture demand from those nets.

The full opt-in H2D engine is available through `//:macro_place`; the separate
`//:macro_hhds` target supplies hierarchy, array tiling, native register-hop
dataflow and the seed. `//:macro_report` and `//:macro_render` handle JSON and
PNG output. The core target remains HHDS- and JSON-free.

```cpp
#include "dreamplace/native/macro/h2d.hpp"
namespace macro = DreamPlace::native::macro;
macro::Parameters parameters;
parameters.K = 4;
parameters.search = "lite"; // or "none", "sa"
std::vector<macro::NodeInfo> names; // one entry per physical node
// Fill names from stable hierarchical paths, geometry/master classification,
// sequential/IO classification and technology orientation masks.
macro::Topology guidance; // optional typed groups/constraints
const auto result = macro::place_h2d(design, names, parameters, guidance);
```

H2D runs grouped GP, GroupLegalize, physical HPWL Compact, optional sequence-pair
search and a fixed-macro cells pass. It keeps the best legal iteration, honors
fixed geometry and transforms pins on flips. Pseudo-nets use hierarchy dataflow
and refreshed cluster centroids; they never enter the final physical netlist.
RUDY-based standard-cell density inflation is separately opt-in. The trimmed
standard-cell filler width is H2D-only; regular placement defaults are unchanged.

`Design::node_boxes` and `order_edges` are shared GP hooks. The optional third
`place(design, options, &state)` argument accepts a `GpState` continuation with
full momentum, BB/density history and filler state. Use `schedule_iterations`
for a fixed schedule when splitting a GP run into checkpoints. Incompatible
solver geometry/connectivity is rejected without changing the design or state.

LiveHD enables this with `synth.place.macro=h2d`, exposes every parameter under
`synth.place.macro.*` and writes `macro_topology.json` plus `h2d.png` under its
placement workdir. JSON is output/cache only; future hierarchy attributes will
supply user guidance through the typed seam. See [`h2d_todo.md`](../../h2d_todo.md)
for stages, controls, units, validation and pending external QoR evaluation.

`Placer` owns its typed `Options`; hosts can configure an instance once, change
its options between calls, or create multiple independently configured instances.
An existing `Options` value can also be passed to `Placer(options)` or assigned
through `placer.options()`. The caller retains ownership of `Design` and any
HHDS mapping. For HHDS, pass `exported.design` to `placer.place()` and then
call `apply_positions(exported, writer)` as below.

This is a reusable configuration object: each call starts fresh solver scratch
state from the supplied coordinates. It does not retain optimization momentum
or provide incremental-placement caches. Explicit `GpState` continuation is a
separate opt-in API. The free function
`place(design, options)` remains available. Both APIs expose all current native
controls in `Options`; upstream Python-only knobs require a native implementation
before they can be added here. No CLI, JSON, Python or file I/O is required.

Independent `Placer` instances can run concurrently on independent `Design`
objects with different options. The linked native kernels keep their placement,
random-generator and LP state local to each call. Configure each instance before
launching it; host callbacks must own their mutable data or synchronize it.
HHDS export/writeback is a separate host operation: exporting independent arrays
does not make concurrent edits to a shared HHDS graph/library safe.

```cpp
// a_design and b_design are independent owning arrays, e.g. two HHDS exports.
DreamPlace::native::Placer a, b;
a.options().num_threads = b.options().num_threads = 1;
a.options().target_density = 0.7;
b.options().target_density = 0.8;
auto first = std::async(std::launch::async, [&] { return a.place(a_design); });
auto second = std::async(std::launch::async, [&] { return b.place(b_design); });
auto a_result = first.get();
auto b_result = second.get();
```

Include `<future>` for this example. `num_threads` is the OpenMP worker count
inside a placement, separate from the host's concurrent placement jobs. Use 1
when the host already fills the available cores, or budget jobs times workers
when using both levels. Parallel tests compare four simultaneous instances
against serial references and exercise simultaneous two-worker OpenMP teams,
including the macro LP path.

`place()` validates all maps before entering raw-pointer kernels. It uses
private work arrays and internal fillers; exceptions leave the caller's arrays
unchanged. Fixed nodes are checked after placement and never written back.
Malformed inputs raise `invalid_argument`; failed legalization raises
`runtime_error`. The returned metrics distinguish global-placement convergence
from final legality. `overflow` is the physical bin overflow at the end of GP,
excluding fillers, normalized by movable area and measured against
`target_density` times available bin area. `hpwl` is weighted, masked HPWL after
the requested stages. A legal result does not imply the GP overflow goal was met.

The global loop uses the upstream merged weighted-average wirelength and
stretched-density/electric-force CPU kernels, with C++ Nesterov extrapolation,
preconditioning, gamma/density-weight schedules and filler insertion. The
Poisson solve uses a separable orthonormal DCT, with no FFT dependency. It costs
`O(nx * ny * (nx + ny))`; the default grid is 32 by 32. This plain DCT is a
portable starting point, not a large-design performance claim.

Full mode applies Hannan/LP macro legalization, greedy and Abacus cell
legalization, then k-reorder, auction-based independent-set matching, and global
swap. Both a separate site/row/overlap checker and the upstream checker validate
legalization. The detailed kernels (k-reorder, independent-set matching and
global swap) use the net weights, and a pass is kept only if it stays legal and
does not increase the weighted, masked HPWL. `Options::detailed_iterations`
(default 2, range 1 to 100) caps the iterations of each kernel, and
`Result::detailed_history` records every pass with its HPWL before and after, a
legal flag and an accepted flag. The native direct LP reuses DREAMPlace's
constraint graph for both initial macro legalization and LP refinement.

Further opt-in GP controls: `min_iterations` (default 40) is the floor before
density convergence is considered; `wirelength_gamma_end` (default 0.05, in
(0, 1]) is the final/initial wirelength smoothing ratio;
`fixed_macro_overlap_weight` (default 0) adds a bounded exact-footprint penalty
between movable cells and fixed macros; `exact_cell_density` charges movable
cells by their physical footprint; and `density_upper_bound` is an experimental
filler-free penalty where underfull bins contribute no charge. A
`GpUpdate::replacement` swaps the standard-cell netlist at a checkpoint:
`previous_node` maps each replacement node to its old index (-1 for inserted
cells), and the other `GpUpdate` vectors must then be empty. Without halos the
GP checkpoint contract is unchanged; with halos the replacement is validated in
physical units before it is padded.

For GP only, set both `legalize=false` and `detailed_place=false`. Coordinates
may overlap, and `Result::legal` is an observation, not a promise. To legalize
an existing placement without GP, set `max_iterations=0`.

## HHDS adapter

Include `dreamplace/native/hhds.hpp`. `from_hhds()` accepts a flat, mapped HHDS
body and three caller-owned geometry functions. This keeps the fork independent
of LiveHD's LEF/Liberty database and placement attribute definitions.

```cpp
hhds::GraphLibrary library;
library.load(database_directory);
auto io = library.find_io(top_name);
if (!io || !io->has_graph()) throw std::runtime_error("missing top body");
auto graph = io->get_graph();

DreamPlace::native::HhdsGeometry geometry;
geometry.cell = [&](hhds::Node node) {
  // Look up the mapped cell's dimensions, current DBU position and status.
  return physical_cell_geometry(node);
};
geometry.pin_offset = [&](hhds::Pin pin) { return physical_pin_offset(pin); };
geometry.terminal_position = [&](hhds::Pin pin) { return top_pin_position(pin); };
// Optional: geometry.net = [&](hhds::Pin driver) { return NetProperties{weight, enabled}; };
auto exported = DreamPlace::native::from_hhds(*graph, layout, geometry);
DreamPlace::native::Placer placer;
placer.options().target_density = 0.8;
auto result = placer.place(exported.design);
DreamPlace::native::apply_positions(exported, [&](hhds::Node node, double x, double y) {
  // Round as appropriate and update LiveHD place_t/status here.
  write_physical_position(node, x, y);
});
```

Every physical driver pin becomes one hypernet, including multi-output cells
and fanout. Connected top pins become separate fixed zero-area terminals.
Constants have no physical driver and their nets are excluded; insert physical
tie cells first when those wires must be placed. Leaf black-box GraphIO stubs
are accepted; subnodes with bodies are rejected so that hierarchy cannot be
silently interpreted as flat connectivity. Hierarchy flattening belongs to
LiveHD. The graph/library must outlive the exported handles. Do not mutate
structure while exporting or writing positions. Writeback rejects a graph
whose structural epoch changed since export and only visits movable cells.

## Import into a host that already uses HHDS (LiveHD)

Use a pinned `http_archive` for this fork and its shipped BUILD files. This
keeps `@hhds` and `@hlop` bound to the **host's existing repositories**, avoiding
a second HHDS implementation or attribute system in the same process. The host
supplies its normal C++23 toolchain and dependencies; archive `.bazelrc` files
are not imported.

```starlark
# Host MODULE.bazel (insert the published commit/archive checksum).
http_archive = use_repo_rule("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
http_archive(
    name = "dreamplace",
    urls = ["https://github.com/masc-ucsc/DREAMPlace/archive/<commit>.tar.gz"],
    strip_prefix = "DREAMPlace-<commit>",
    sha256 = "<archive sha256>",
)
bazel_dep(name = "coin-or-lemon", version = "1.3.1.bcr.1")
bazel_dep(name = "openmp", version = "21.1.5.bcr.3")
single_version_override(
    module_name = "coin-or-lemon",
    version = "1.3.1.bcr.1",
    patches = ["//packages:dreamplace_lemon_graph.patch"],
    patch_strip = 1,
)
# Host BUILD: deps = ["@dreamplace//:hhds"] or ["@dreamplace//:dreamplace"].
```

Copy `packages/lemon_graph.patch` into the host as
`packages/dreamplace_lemon_graph.patch` and export it from that package. Bazel
requires module-override patches to live in the root workspace, so a label
pointing into `@dreamplace` is not accepted here. No overlay
`packages/dreamplace.BUILD` is needed. The host retains its existing
rules_cc, platforms, HHDS, HLOP, Abseil, iassert and emhash declarations. For a
core-only consumer those HHDS dependencies are not needed. Enable CPU threading
in a host with `--define=dreamplace_openmp=true`; only the standalone repository
owns the convenient `--config=openmp` alias. If importing as a Bzlmod module
instead, the root module must still apply the LEMON graph-target patch, since
Bazel only honors version overrides from the root module.

## Synthetic HHDS fixtures and visual debugging

```sh
# Generate and place all five fixtures, then export physical files and images.
bazel run -c opt //dreamplace/native:synthetic_place -- --out results/synthetic

# Compare against the random-start mesh using a mesh-shaped initial placement.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case mesh-grid --out results/synthetic

# Record actual iterates and create a local playback/seek viewer.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case mesh-grid --animate --frame-step 20 --out results/synthetic

# Increase the core area while keeping cell sizes, connectivity and solver settings.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case mesh-grid --utilization 0.5 --animate --out results/synthetic-u50

# Add row-aligned input/output pairs on the left/right sides of the grid mesh.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case mesh-grid-io --utilization 0.45 --animate --out results/synthetic-u45

# Select a fixture, change the seed, or increase the 4-way expansion depth.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case tree --tree-depth 4 --seed 7 --out results/synthetic

# Read an existing HHDS checkpoint instead of regenerating the connectivity.
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --load results/synthetic/mesh_16x16/initial.hhds --out results/reloaded

# Optional parallel CPU run.
bazel run -c opt --config=openmp //dreamplace/native:synthetic_place -- \
  --threads 4 --out results/synthetic-openmp
```

| Default fixture | Cells | Nets | Connectivity |
| --- | ---: | ---: | --- |
| `mesh_16x16` | 256 | 960 | Four nearest neighbors; bidirectional links; open grid boundaries |
| `mesh_16x16_grid` | 256 | 960 | Same mesh, initialized as a regular grid instead of random positions |
| `mesh_16x16_grid_io` | 256 | 992 | Same grid-start mesh plus 16 left inputs and 16 right outputs, one pair per logical row |
| `torus_12x8` | 96 | 384 | Four nearest neighbors; bidirectional links; wraparound in both dimensions |
| `tree_4to2_d3` | 148 | 212 | One input, three 4-way expansion levels (64 leaves), six pairwise reduction levels, one output |

Each directed link is a separate two-pin net with unit weight. The tree's count
includes its two IO nets. These fixtures exercise connectivity and placement;
no timing model or STA dependency is needed. The tree has fixed zero-area IO
anchors at the left/right die-edge midpoints (input/output), outside the 2um
core margin. For the default tree these are (0, 11) and (32, 11) microns.
The PNG shows the margin and labels these anchors IN/OUT; DEF exports fixed
metal1 pin shapes extending inward from the die edges.
`--mesh-width/--mesh-height`,
`--torus-width/--torus-height`, and `--tree-depth` configure graph size. The
original fixtures start at deterministic random positions independent of logical
grid/tree coordinates. The additional `--case mesh-grid` uses logical row/column
coordinates to seed an evenly spaced grid over the same core, with identical
cell sizes and connectivity. Every cell remains movable; there are no persistent
relative constraints. Its 16 columns of 2um cells slightly overlap in the 29um
core, so final legalization can increase HPWL over that seed. Compare final
legal placements between the two mesh cases. `--mesh-width/--mesh-height`
configure both meshes; `_grid` keeps their output directories separate.
`--case mesh-grid-io` adds one input net to each leftmost cell and one output
net from each rightmost cell. The original mesh links, cell sizes, core and
starting coordinates are preserved. Fixed zero-area pins `in_0..in_15` and
`out_0..out_15` sit on the left/right die edges at matching initial row-center
heights (row 0 at the bottom). There are no top/bottom pins. Dimensions control
the pin count; each added two-pin IO net has unit weight. The `_grid_io`
directory keeps this variant separate. The image labels each pin, and DEF
preserves its declared HHDS name. All mesh cells remain movable.
`--utilization` sets requested cell area / core area for generated fixtures
(0.1..0.8, default 0.6), with site/row rounding. Lower values give more space;
the default mesh grows from 29x30um to 32x32um at 0.5. Use a separate `--out`
directory to retain comparisons. This is separate from the solver's 0.8 target
density, which remains unchanged. Saved `--load` fixtures retain their floorplan.

The fixture flow defaults to a 16x16 density grid, target density 0.8 and up to
1,500 GP iterations, followed by legalization and detailed placement. This
iteration budget lets all five default fixtures reach the 0.1 overflow goal;
the library's general defaults are unchanged. `--iterations`, `--bins`, and
`--threads` override the fixture settings. The output path is relative to the
workspace when invoked with `bazel run`.

The fixture CLI also exposes the library's existing GP controls:

| Flag | Default | Effect |
| --- | ---: | --- |
| `--stop-overflow` | 0.1 | Lower values require less physical bin-density overflow before GP stops |
| `--target-density` | 0.8 | Sets bin capacity and filler area; lower values encourage more whitespace |
| `--density-weight` | 0.01 | Initial density-force multiplier, normalized against wirelength forces at startup |
| `--max-step` | 0.5 | Bounds the step-size calculation relative to bin dimensions; Nesterov extrapolation also moves cells |

`--iterations` is a budget, not a convergence guarantee; check `gp_converged`
and `gp_overflow` in `qor.json`. These flags affect global placement. The native
legalizer currently runs macro legalization, greedy cell packing, then Abacus
within the rows chosen by greedy. Its packing preferences are inherited from
upstream and are not exposed as a local repulsion or displacement-limit knob.
The legacy Python flow separately exposes `abacus_legalize_flag`.

Bin-density convergence does not imply nonoverlapping rectangles or legal
row/site coordinates. The animation now reports physical overlapping pairs
and their summed pairwise intersection area in its frames, viewer and CSV.
Pairwise area is not union area; zero-area terminals and fixed-fixed pairs are
excluded. This separates the remaining overlap from the density stopping test.
A stricter example, retaining the baseline files, is:

```sh
bazel run -c opt //dreamplace/native:synthetic_place -- \
  --case mesh-grid-io --utilization 0.45 --stop-overflow 0.02 \
  --iterations 4000 --animate --frame-step 40 --out results/synthetic-stricter-long
```

Each fixture directory contains:

- `initial.hhds/` and `placed.hhds/`: real HHDS libraries and a `fixture.txt`
  manifest. A fixture-specific persisted node attribute stores coordinates,
  logical provenance and placed status. It is separate from LiveHD's attrs.
  The V3 manifest records side IO along with layout and initialization; saved
  V1/V2 fixtures remain readable.
- `cells.lef`: synthetic two-layer technology, a 1um site and one 2um by 1um
  `SYNTH4` macro with physical `I0..I3` / `O0..O3` pin rectangles.
- `placed.def`: die, rows, placed components, fixed IO and net connectivity at
  1,000 DBU/um. LEF describes the library geometry; DEF carries the placement.
- `placement.png`: 1600x900 comparison of initial seed and final placement,
  with colors for logical grid columns or tree stages, unrouted net lines and
  HPWL/convergence labels. The image uses the reloaded, placed HHDS checkpoint.
- `positions.csv` and `qor.json`: individual coordinates and metrics/settings,
  including the initialization mode and whether the initial placement was legal,
  distinguishing GP overflow/convergence from final legality. `internal_hpwl`
  excludes nets connected to top IO, and `io_hpwl` reports the remainder so the
  additional IO nets do not obscure comparisons of the original mesh links.

Before physical export, the flow saves/reloads the initial graph, places its
arrays, writes coordinates back through the HHDS adapter, and saves/reloads
the final graph. It checks coordinates, connectivity and legality after reload.
The adapter now retains original pin handles alongside the flat pin arrays,
allowing physical exporters to preserve port identities exactly.

PNG rendering is native C++ with a small RGB rasterizer and PNG encoder using
stored DEFLATE blocks. It introduces no image/compression library, Python,
display server, or GPU dependency. Images are about 4.3 MB each. The fixture
LEF/DEF writers deliberately support this one synthetic technology; a general
LEF database or real technology export remains a separate integration task.

With `--animate`, the flow also writes `animation.html`, `animation.csv`,
`animation.ffconcat` and an `animation_frames/` directory of 800x450 PNGs
(about 1.1 MB per frame). Open the HTML locally to play, pause, step or scrub
through the actual iterates alongside the initial placement. Keep the frames
directory with the HTML. GP is sampled every `--frame-step` iterations (20 by
default), plus its last iteration and the initial, legalized, detailed and final
states. Fillers are hidden. Phase changes pause longer; there is no position
interpolation, and legalization/detailed placement each appear as one transition.
The native trace/viewer has no new dependencies and observation does not change
the placement result. The fixture's core observer API is also available to a host:

```cpp
options.observer_interval = 20;
options.observer = [](const Design& physical, PlacementStage stage, int iteration) {
  // Read-only physical state; copy positions here if retaining the snapshot.
};
```

An already-installed FFmpeg can optionally turn these frames into a portable
video or GIF; it is not a Bazel or native runtime dependency:

```sh
ffmpeg -safe 0 -i results/synthetic/mesh_16x16_grid/animation.ffconcat \
  -vf fps=20 -c:v libx264 -pix_fmt yuv420p -movflags +faststart \
  results/synthetic/mesh_16x16_grid/placement.mp4
ffmpeg -safe 0 -i results/synthetic/mesh_16x16_grid/animation.ffconcat \
  -filter_complex 'fps=10,split[a][b];[a]palettegen[p];[b][p]paletteuse' \
  -loop 0 results/synthetic/mesh_16x16_grid/placement.gif
```

`synthetic_mesh_test`, `synthetic_mesh_grid_test`, `synthetic_mesh_grid_io_test`,
`synthetic_torus_test`, and `synthetic_tree_test` run
the complete flow. They check exact topology (including torus wraparound),
reachability, checkpoint persistence, GP convergence, legal placement, HPWL
improvement over random seeds, grid seed ordering/alignment, per-row IO connectivity,
fixed side anchors and exported component/pin/net identities.
The DEF subset is independently read back in the tests. The initial generated
LEF/DEF files have also been parsed with the independent Si2 5.8 readers in the
sibling LiveHD checkout; these readers are not build dependencies here.

## Validation and remaining integration work

The tests cover weighted wirelength finite differences, Poisson eigenmodes,
exact HPWL, malformed CSR, locked/fixed obstacles and terminals, deterministic
serial runs, threaded runs, GP-only and legalization-only modes, multi-row
macros, negative core origins, failed legalization without partial writeback,
and export/place/writeback of a persisted HHDS graph. The direct macro LP is
checked against exhaustive integer optima for small problems.

This fork supplies the native engine, adapter and synthetic debug flow. LiveHD
still owns production technology LEF/DEF integration, floorplan generation, IO
slot assignment, orientations, STA/RC integration and incremental cache policy.
No stock-DREAMPlace benchmark QoR parity or large-design scalability is claimed
by these tests.

### Soft color cohesion

`Design::cohesion_groups` describes disjoint groups of physical nodes (for
example, one synthesis color in one instance scope). A group's `weight`
defaults to one. `Options::cohesion_weight` defaults to zero, preserving ordinary
GP exactly. Positive values enable a soft, translation-invariant objective:

```
E = lambda(t) * sum_g weight_g / area_g * sum_i area_i * |center_i - centroid_g|^2
```

Centroids are area weighted and recomputed at the current Nesterov evaluation
positions each iteration. Lambda ramps linearly from zero to `cohesion_weight`
over the existing GP schedule. The coefficient has inverse-coordinate units
(use coefficient/DBU-per-micron when specifying it in inverse microns). Node
centers and physical areas are used; halo padding preserves both the physical
centers and area weights. Fillers and zero-area terminals are excluded. Fixed
members can contribute to a centroid but remain fixed. Singleton groups have
zero energy. A node can belong to at most one group.

There are no group boxes, projections, fixed centroid terminals, or changes to
macro/cell legalization. Group processing is linear in membership count (a
centroid accumulation and a gradient pass). Continuation validates membership,
weights and the coefficient along with the existing GP state identity.

### Repulsion controls (all opt-in)

`Design::node_macro` identifies non-standard-cell macros independently of row
height. `automatic_bins` resolves a 32–64 power-of-two grid from the smallest
physical macro dimension before halo padding (four mean cell pitches when no
macros are identified). `exact_macro_density` removes macro smoothing; standard
cells retain the sqrt(2)-bin floor. The Poisson solver still uses the existing
separable direct DCT loops; no FFT dependency is introduced.

`scale_macro_density` charges movable macros at `target_density`, uses that
same charge in exact overflow, and adjusts filler area consistently. Without
this correction an isolated resolved macro at density 1 intrinsically violates
any target below 1, so raw density overflow is not a macro-overlap count.
`GpSample::raw_overflow` retains the legacy metric for comparison.

`macro_overlap_weight` adds `weight * intersection_area^2 / min_area^1.5`
for spatially hashed macro pairs on halo-padded occupancy boxes. Normalization
makes the weight independent of caller coordinate units. Equal edges choose a
deterministic separating subgradient. `macro_jitter` is a seeded fractional
macro-size initial offset, projected into bounds, and does not alter filler RNG.
Fixed objects never move. As with an exact intersection-area objective generally,
strict containment has locally flat directions; the term is not a legalizer.

`adaptive_density` uses 20-iteration overflow progress to choose 1.4× growth
when improvement is below 5%, otherwise 1.05×; 20% window HPWL growth acts only
as a brake. `max_iterations` remains the explicit safety cap. `history` and an
optional synchronous `progress_observer` expose progress, macro overlaps,
near-coincident center pairs, density weights and inflation. `poisson_seconds`
measures direct transform solve time, excluding wirelength and legalization.
Continuation state retains feedback history and new objective controls are
included in the continuation signature.

For an opt-in replacement of the original routability update, combine
`routability` with `improved_routability`. It starts below 0.3 physical overflow,
retains inflation monotonically, limits each update to 5% of movable area, and
reduces the effective target density with accumulated inflation. Macro/fixed charges and the filler-charge budget are rebalanced at that density; shifting only the uniform target would not change Poisson forces. Only nodes in
`node_inflatable` can inflate. Directional capacity is tracks per layout unit.
`node_blocked_capacity_h/v` supplies the sum of capacities on each macro's OBS
layers; the model conservatively blocks each such layer across the full macro
occupancy envelope. A 1% capacity floor keeps fully blocked bins finite.
`RudyOptions::pin_access` estimates macro signal-pin demand against available
tracks at bin edges. These estimates are not a detailed router or signoff model.

### Optional same-class macro alignment and soft snapping

`Design::node_align_class` uses one integer per node, with `-1` opting out.
The host assigns a distinct class to each LEF master × fully scoped synthesis
color. Only movable macros in `node_macro` participate. Fixed nodes, fillers
and standard cells are excluded.

`macro_align_weight=0` disables GP alignment exactly. Positive weights scale
an equal/opposite axis force by the current mean macro wirelength+density force.
`macro_align_sigma=0` is quadratic; positive sigma gives a bounded well, with
0.25 used in the Eyeriss experiment. Degree normalization reduces pair-count
bias. Pair projections must overlap on exactly one axis. The pair list refreshes
every20 global iterations and is preserved across GP continuation. Changed
active alignment options/classes reject continuation. `macro_align_start=100`
and `macro_align_full=500` define the squared ramp; enabled runs cannot stop
before full strength. `macro_align_overlap_boost=0` and
`macro_align_nearest=false` keep the ordinary overlap weight and all-pairs rule.

`macro_snap=false` preserves stock legalization. When enabled,
`macro_snap_tolerance=0.25` groups nearby same-class macros per axis and rounds
the median origin onto the site/row lattice. Macro anchors change before the
LP; standard-cell anchors do not. Chain equalities use finite-capacity arcs
with `macro_snap_weight=100` against unit displacement. Core/halo constraints
remain uncapacitated hard arcs. A group can break rather than violate geometry;
`Result::alignment` records exact/broken groups, violated equalities and macro
movement from GP. The opt-in macro mask includes small SRAMs below the legacy
five-row classifier and protects their positions through cell LG/detail.
No flips or second cell-global-placement pass are added.

The new tests cover pair selection and color isolation, finite-difference
quadratic/bounded gradients, ramps, continuation mismatches, disabled controls,
small SRAMs, halo spacing and soft-vs-hard conflicts. The difference LP is also
checked against an exhaustive bounded-integer oracle with soft capacities.

`routability_charge_compensation=true` retains the existing charge budget.
False is an explicit compatibility setting for the earlier R4 inflation-only
experiment, not a changed default. Direct DCT loops and the32-bin default are
retained; automatic refinement remains capped at64 per axis.

### Opt-in candidate/color centroids (R8H2)

The rejected R8H template engine is removed. R8H2 starts its first GP from scratch
with moving-centroid attraction, then exports legal macro locations for the usual
fixed-macro second placement pass. No preliminary macro GP or continuation phase
is used. The default centroid weight remains zero.

`detail::centroid_candidates(design, paths, masters, minimum_area,
minimum_copies = 4)` finds relative macro path/master/size structures with more
than two movable macros and physical logic-plus-macro area strictly above the
caller's equivalent-cell threshold. The unit is mean mapped standard-cell area,
or site × row area if none exist. At least `minimum_copies` qualifying copies must
share one immediate hierarchy parent; copies across parents never pool. The
copy threshold must be at least two. Count siblings before collapsing unary
wrappers, so an ineligible singleton inner wrapper does not hide repeated outer
instances. Reports expose `copy_parent` and `copies_per_parent`; numeric `parent`
still identifies the nearest retained candidate. The family identity remains
structural across qualifying parent scopes for same-type repulsion. Host color
keys are scoped synthesis identities, not numeric colors reused across unrelated
instances. LiveHD exposes the copy threshold as `synth.place.fp.min_copies` (4).

Each candidate includes every full color seen inside its hierarchy, including
colored cells outside that hierarchy, plus its uncolored descendants. The current
physical-area-weighted centroid uses all those macros and standard cells. A color
in n candidates has n anchors toward the n candidate centroids. Membership overlap
is intentional. Parent ids describe hierarchy; no ownership or hard bounds remove
any anchor. Identical memberships can have identical centroids.

The variance energy per group is lambda/A times sum of a*|p-centroid|^2. Color
moments are computed once per iteration, candidate coefficients are accumulated
per color, and a final node pass applies their sum. This is exactly the sum of all
anchors without scanning a large shared color once per candidate. Physical-area
metadata survives halo padding. The centroid derivative cancels, forces sum to
zero, and translation does not introduce an absolute anchor.

With nonzero `Options::centroid_weight`, the attraction ramps over
`centroid_ramp` iterations (default100) in ordinary GP, sharing its normal stop/cap,
repulsion, density and routability controls. `centroid_start` observes iteration0;
`Result::centroid` reports iterations and force time. External GP continuations with
active groups are rejected. `centroid_test` checks shared n-anchor additivity,
outside standard cells, finite-difference gradients, physical halo weights,
default-off identity, fixed positions and from-scratch initialization.


### Same-type centroid-circle repulsion (R8H3)

`Options::centroid_repulsion_weight` is opt-in and defaults to zero, preserving
attraction-only and ordinary GP. It requires active centroid attraction. The
LiveHD R8H3 experiment supplies weight 5 and a >1000-equivalent-cell eligibility
threshold; the native candidate factory accepts the caller's threshold.

`detail::centroid_cell_area` supplies the common mean physical standard-cell
area unit. Both candidate eligibility and circle sizing count physical cells and
macro-equivalent area. For each candidate, its full anchored physical area A
(including shared colors outside the hierarchy) reserves circle area A/rho at
the requested target density: radius = sqrt(A/(pi*rho)). Each node is counted
once within a candidate; its area contributes to every candidate attracting it.
Halo padding and temporary fillers do not change this area.

Only equal-family pairs interact, with energy
`lambda * max(0, r_i+r_j-||c_i-c_j||)^2`. Finite-range gradients are propagated
through the moving centroids using a_i/A_g and summed for shared colors. Equal
centers get a deterministic separating subgradient. Repulsion shares the
attraction ramp and normalization; the requested density, rather than the
routability-adjusted density or macro-charge scaling, sizes the circles. This
adds no hard bounds, fixed centroid terminals, or individual-cell pair forces.
The attraction-only and fixed-penalty paths retain density-based stopping. Tests cover the analytic
gradient, type filtering, cutoff, area/density scaling, halo invariance, shared
memberships, force balance, deterministic ties and fixed nodes.


R8H3 now enables `centroid_adaptive_repulsion` in LiveHD. The geometry is unchanged:
repulsion persists until r_i+r_j (two radii for equal circles). Every 20 iterations
after the attraction ramp, the penalty doubles when maximum relative circle
penetration improves by less than 10%, up to 10000 times its initial weight.
The density weight is held when its own target is satisfied and only circle
separation remains pending. GP can terminate successfully only when both density
and `centroid_stop_overlap` (default 0.001 of r_i+r_j) are satisfied, or stops
unconverged at the existing iteration cap. The option defaults off in the native
API; ordinary placement and attraction-only modes preserve their prior behavior.
`GpSample` and `CentroidReport` expose actual overlap, effective penalty and circle
convergence. The equal-radius regression checks nonzero force at 1.5r, zero
outside 2r, and optimization of a feasible overlapping pair.


For a nested candidate, parent attraction acts through every child member as well
as the child's own attraction. Summing the parent gradient over the child gives
2*lambda*A_child/A_parent*(c_child-c_parent), so the child as a whole is pulled
toward the parent. The child remains deformable. Different-type parent/child
pairs have no repulsion; same-type copies do. The nested regression checks this
net force and verifies that adding the parent creates no repulsion on its child.


### Circle centroid snapping and frozen targets

`floorplan.hpp` provides the pure native `snap_floorplan(circles, bounds, options)`
API. It performs no file I/O and does not run GP. Lines, neighbor spacings and
symmetry axes are fitted from the supplied positions; no global lattice is given.
The default alignment threshold is strictly one smaller-circle radius, permitting
half-pitch staggered rows. Scalar penalized dynamic programming groups similar
neighbor gaps. Alignment, spacing and symmetry relations require the same type
and immediate parent candidate; roots share the implicit design scope. Different
parent candidates retain independent snap lines, gap classes and symmetry axes.
Hard same-type circle separation still applies across parents.
Accepted exact equalities are fitted by area-weighted Dykstra
projection with same-type circle-separation halfspaces along initial pair normals,
core center bounds, and local movement budgets. Parent translations carry child
initial targets before child types are fitted. Complete circles may extend outside
the core; these are interaction ranges, not cell fences.

The relationship selector is greedy and the projection solver is bounded; rejection
does not prove mathematical infeasibility. Failed families are retained and counted
in `infeasible_families`, and `feasible` is false if separation is unresolved.
`FloorplanSnap` reports target geometry, accepted relations, rejected attempts,
maximum penetration/displacement and time. This is an experimental local heuristic,
not a globally optimal symmetry or packing algorithm.

`Design::centroid_targets` optionally supplies frozen world-coordinate targets,
one per `centroid_groups` entry. With `Options::centroid_weight>0`, attraction uses
those targets. Empty targets retain the original moving-centroid energy and gradient.
Frozen-target gradients include translation; moving-centroid gradients still sum to
zero per group. The caller disables centroid repulsion for refinement. Targets do
not enforce exact member means or constrain hard region bounds. Shared full-color
memberships can couple several targets; independent point feasibility does not
establish that every target can be realized by those members.

LiveHD exposes this as opt-in `synth.place.fp=true`: floorplan GP without LG/DP,
snap, then warm-start frozen-target GP and native LG/DP. The Eyeriss runner separately
fixes SRAMs for the normal final cell GP without anchors. All evaluation GP grids
are 64×64; native/default options remain unchanged. `floorplan_test` covers a mesh,
staggered triangle, centered tree, projection threshold, bounded-span chains,
independent families, inherited parent translations, unequal radii, hard-conflict
reporting, unit scaling, and frozen-target finite differences. `floorplan_demo`
emits the three pattern experiments as JSON for external plots.

Eyeriss: all 204 targets separate in 1.15 s, accepting 123 line, 36 spacing and
32 midpoint relationships. Weak-target refinement then reintroduces measured-mean
circle overlaps. Final legal HPWL is 6.578 m versus 4.595 m for the reused R8H3
reference (+43.15%). Extra refinement GP and snapping change together in this
comparison. Keep off by default; no new routing/timing result is claimed.

### Bounded continuation after density convergence

`Options::post_density_factor` is opt-in (zero preserves the legacy path). After
geometry/density first meet their stopping conditions at iteration N, the
solver allows factor × N additional iterations, still bounded by
`max_iterations`. Weighted HPWL and inflation must settle over two consecutive
20-iteration samples within `convergence_tolerance` (default 1%), and current
RUDY overflow must be zero when enabled. Host checkpoints can supply a
nonnegative `GpUpdate::timing_violation`; a positive setup violation prevents
early convergence. Without host timing, weighted HPWL is only a timing proxy.

Weight/topology updates reset stability but never the deadline. External
`GpState` continuation retains the deadline and sample history. The result
distinguishes objective convergence, host timing availability and continuation
budget exhaustion. LiveHD enables this policy with factor four by default.

### Density-gated host timing checkpoints

`Options::checkpoint_after_density` is opt-in (false by default). With an active
checkpoint interval, it defers the host callback while physical overflow exceeds
`stop_overflow`. The measurement uses the same halo-padded geometry and macro
charge scaling as GP convergence, before a halo wrapper unpads the host snapshot.
GP cannot declare early convergence before its first eligible callback. The
iteration cap and bounded-continuation deadline force a callback even above the
threshold, including when the cap is off cadence. Interval zero disables all
callbacks. Initial caller-supplied timing weights remain active.

`Result::checkpoint_history` records iterations, measured overflow, invocation
and forced-cap flags. An optional `checkpoint_observer` receives these decisions
before the host callback, allowing a host movement filter to honor forced final
samples. The gate itself neither calls a timing engine nor resizes cells.
