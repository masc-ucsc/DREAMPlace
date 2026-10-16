# H2D macro placement

The current native stages are implemented in this checkout and integrated
with LiveHD. Eyeriss evaluation has exposed the quality and performance gaps
listed below. There is no separate rev2 specification. The user approved the
full implementation, including TopoSearch, pseudo-nets and opt-in routability,
without waiting for a real-design quality gate. Real-design evaluation remains
external in HighTide; passing synthetic tests is not a QoR win over OpenROAD.

## Confirmed boundaries

- Native C++23 and Bazel; no Python, Torch, CUDA or OpenROAD runtime dependency.
- `//:dreamplace` stays HHDS-free. `//:macro_hhds` reads the caller's hierarchical
  HHDS graph before flattening. `//:macro_place` is the flat realization engine.
- Hierarchy grouping, array tiling, shape curves and annealing reimplement the
  ideas independently. No OpenROAD placer code was ported. No unlicensed
  Re2MaP code enters the tree; any RM comparison is external and run-only.
- C++ `Topology` is the source of truth. Every stage accepts/merges typed
  topology, keyed by hierarchical names. JSON is automatic output/cache under
  the LiveHD workdir, never user input. There is no topology CLI input or emit
  slot and no standalone ORFS/DEF driver.
- Future HHDS attributes will supply hierarchy/leaf constraints, with leaf
  precedence. Reading those attributes is explicitly future work.
- Orientations are N/FN/FS/S, with real pin transforms and LEF symmetry masks.
  No 90-degree rotations.
- LiveHD selects `synth.place.macro=none|h2d`; the existing default stays `none`.
  H2D requires full placement. All native H2D controls are available under
  `synth.place.macro.<name>`.
- Existing placement/filler defaults stay unchanged. The trimmed standard-cell
  filler width is used only within the opt-in H2D calls. If it would exceed the
  existing one-million-filler limit, widen only those private rectangles to
  preserve target area within the limit.
- Tests use synthetic inputs declared in their own repository. Real designs,
  placement/pre-route timing and post-route comparison belong in HighTide.

## Implemented stages

| Stage | Code and behavior |
| --- | --- |
| Geometry | `orientation.*`, `macro/check.*`: four flips, composition, physical pin transforms, strict macro/fixed-object overlap, pitch, bounds and halo checks. |
| Routing proxy | `rudy.*`: native directional projected-box demand, WDM degree factors, masks, explicit track capacities and mean positive bin overflow. |
| Typed topology | `macro/topology.*`: groups, members, anchors, relative order, fixed macros, blockages, dataflow, realized sequence pair/graph, validation and merge. Constraints carry hard/soft strength and weight. |
| Hierarchy | `macro/hierarchy.*`: contextual traversal and boundary-pin connectivity, including repeated shared module definitions; instance hierarchy is captured before LiveHD flattens. No shortcut through a module from its input to its output. |
| Dataflow | `macro/dataflow.*`: bounded register-hop BFS, combinational-cycle deduplication, fanout cap, macro/IO stopping endpoints and hierarchy register-cluster endpoints. Weight is `1/2^h` for the shortest connection to an endpoint. |
| Seed | `macro/seed.*`: parent/master/dimension buckets, fixed singletons, divisor-grid shape curves, splitting when no shape fits, native position/shape annealing, IO/nearness anchors, weighted order and transitive reduction. |
| Abstract/Expand | `macro/abstract.*`: rigid group nodes, remapped physical pins, disabled intra-group nets, fixed suffix and blockages. Expand preserves exact physical sizes and fixed inputs and transforms pins for selected flips. |
| GP hooks | `placement.*`, `constraints.hpp`: boxes projected at every GP snapshot; weighted squared-hinge order force; schedule/seed; full continuation state including fillers, momentum, density normalization and BB history. Topology mismatch rejects stale state. |
| Difference LP | `macro/diff_lp.*`: integer pitch units, weighted L1/linear terms, hard/soft arcs, bounded numeric quantization, LEMON network simplex dual and named hard negative-cycle diagnostics. No general LP/GLPK dependency. |
| GroupLegalize | `macro/compact.*`: geometry graph and L1 realization, explicit hard order, optional soft order arcs, lowest-weight soft-edge relaxation, existing native macro/Hannan fallback plus deterministic conservative rectangle-packing recovery, exact fixed objects and halo envelopes. Hard guidance is never silently discarded by fallback. |
| Compact | Physical weighted HPWL objective, collapsed one-group hinge terms and multigroup net bounds; outline/obstacle/soft-anchor slack; actual realized H/V graph retained. Coarse pitch preserves off-coarse-grid fixed coordinates. |
| TopoSearch | `macro/toposearch.*`: geometry SP and prior SP seeds, greedy lite and serial SA chains, plus/minus/both swaps, symmetry-allowed flips, frozen relations, O(G^2) packing screen, delayed RUDY acceptance, normalized cost and best legal fully scored candidate. |
| Pseudo-nets | Direct group/group and group/IO endpoints; fixed zero-area standard-cluster centroid terminals refreshed each outer iteration; threshold and beta controls. They are removed before physical Compact and the cells pass. |
| H2D driver | `macro/h2d.*`: continued GP, GroupLegalize, optional search, Compact, fixed-group cells legalization/detail, common best-iterate merit, freeze/hysteresis/early stop and L22 coarse-snap trial. A failed later cells/snap trial preserves the prior legal best. |
| Routability | Opt-in GP density occupancy inflation for marked movable standard cells only. Physical sizes/pins remain unchanged; legalization uses original dimensions. No macro or filler inflation. |
| Reports | `macro/report.*`: versioned integer-DBU topology cache, fingerprints, resolved knobs, per-round HPWL/overflow/GP state, constraint violations/slack, accepted moves, LP/screen counts, selection and phase times. `macro/render.*` writes native PNG with macros/cells, halos, anchors, order arrows and IOs. |
| LiveHD | `pass/place/h2d.*`: technology-independent hierarchy/Liberty classification, LEF symmetry and routing pitch capacities, automatic macro cache, FIXED frame before ordinary block-cache/ECO logic, JSON/PNG and QoR block. Timing-seed snapshots skip H2D. |

## Algorithm and units

1. Read the actual hierarchical graph and classify macro, sequential, logic and
   IO endpoints using host callbacks. Build a stable physical-name mapping.
2. Seed macro arrays and topology, merging typed input guidance. Fixed input
   coordinates/orientations take precedence. Missing old macros are pruned and
   new leaves are seeded. Unnamed leaves can match by parent, pin/size signature
   and ordinal, with a diagnostic; this is a fallback, not a stable naming claim.
3. Abstract groups, refresh pseudo terminals, and continue GP with boxes/order.
   Preserve full GP history when solver geometry/connectivity matches; reset it
   when abstraction shape/connectivity changes.
4. Legalize group geometry, search topology when enabled, and minimize physical
   HPWL with Compact. Ordinary movable standard cells are never LP obstacles.
5. Freeze groups, remove pseudo-nets, run cells-only GP/LG/DP and expand to the
   physical netlist. Strictly audit macro and fixed-object legality.
6. Retain the best legal iteration by
   `physical_HPWL + alpha * (core_width + core_height) * RUDY_overflow`.
   Search additionally uses displacement, order, anchor and feasibility costs;
   RUDY is counted once in its delayed-acceptance cost.
7. Optionally snap through a coarse-pitch Compact/cells trial; accept only legal
   candidates with no worse common merit. Publish the selected result atomically
   at the API boundary. Impossible packing throws without changing the input.

All native geometry uses one caller-chosen layout unit; LiveHD uses DBU. Native
track capacities are tracks per layout unit, hence tracks per DBU in LiveHD.
LiveHD derives horizontal capacity from horizontal layers' perpendicular Y
pitches and vertical capacity from vertical layers' perpendicular X pitches.
Missing pitches require explicit `route_cap_h/v`; no technology names are guessed.
The existing `synth.place.halo` option remains in microns. H2D `min_channel`,
`delta` and `pitch_x/y` are layout units (DBU). Pitch overrides must be whole
site/row multiples. A minimum channel reserves at least the requested gap, rounded
outward; explicit halos and channels use the maximum reservation rather than
adding the same margin twice. With neither, seeding reserves 10 sites/2 rows.

RUDY is a pre-route estimate. Zero-area net bboxes produce zero projected-box
demand and are counted separately. It does not model pin access, vias, clock
routing, detailed congestion or routed delay. Routability inflation is opt-in
and is not evidence that physical routing will complete.

The cache uses labeled noncryptographic FNV1a fingerprints. Macro geometry,
pins, hierarchy names, fixed input state and technology identify the frame;
standard-cell connectivity can change while the frame is reused. Input guidance,
all resolved H2D parameters, IO geometry and placement implementation sources
participate in host cache identity. Macro pin signatures use canonical N offsets,
so replaying a flipped output does not spuriously invalidate the input. Generated
FIXED macros in a persisted floorplan remain cache data; explicit user FIXED
constraints stay authoritative. A malformed/stale frame causes a fresh placement
and records its miss reason. JSON is not an agent editing interface.

## Controls and ablations

`macro/params.def` is the single parameter registry used by C++ parsing, reporting
and LiveHD CLI registration. Use `lhd describe synth.place.macro.<name>` for current defaults/help.

| Area | `synth.place.macro.*` suffixes |
| --- | --- |
| Stages | `anchors`, `order`, `loop`, `pseudonets`, `toposearch`, `routability` |
| Loop | `K`, `iters_per_outer`, `eps`, `freeze_after`, `hysteresis` |
| Seed/dataflow | `seed`, `seed_sa_evals`, `min_channel`, `delta`, `tau_q`, `hops_max`, `fanout_cap`, `theta_q`, `beta` |
| GP | `lambda0`, `gamma`, `bins`, `target_density`, `threads`, `route_cap_h`, `route_cap_v` |
| Realization/search | `alpha`, `mu`, `rho`, `kappa`, `M`, `anchor_mode=exact|soft`, `order_in_lp`, `normalize_at_seed`, `pitch_x`, `pitch_y` |
| Search | `search=none|lite|sa`, `chains`, `evals`, `T0`, `cooling`, `move_plus`, `move_minus`, `move_both`, `move_flip`, `orient_moves` |

Density zero inherits the caller's placement density. LiveHD inherits global
placement bins/threads unless their macro-specific options are explicit.
One-thread runs are deterministic for the same input, build and seed. SA chains
run serially with independent derived seeds. Do not assume bit identity across
OpenMP worker counts or toolchains.

## Validation checkpoint (2026-10-16)

- All 20 native optimized test targets passed. Existing native placement and
  synthetic mesh/tree/torus tests remain green.
- Six debug targets passed: GP hooks, Abstract/Expand, Compact, H2D, H2D scenarios
  and HHDS hierarchy. This is focused debug validation, not the full debug suite.
- Eight final focused targets passed with AddressSanitizer and
  UndefinedBehaviorSanitizer together.
- Independent LP tests compare exhaustive optima, legacy L1 LP, soft penalties,
  hard cycles, fixed coarse-grid residuals and a 50-group/5,000-net case.
- Native fixtures cover a pipeline, register hops, repeated hierarchy, a 4x4 SRAM
  array at 64% physical utilization, edge IO anchors, ring/scrambled-ring,
  tight soft/hard anchors, flips, fixed inputs and impossible packing.
- Exact serial determinism, 40 iterations versus 20+20 continued GP, best K>1
  merit versus K=1, topology round-trip and stage ablations are tested.
- PNG chunks/CRC/zlib data were independently decoded and the image inspected.
- All 27 optimized targets passed through LiveHD: 20 native tests plus
  placement, repair, physical LEF/DEF, normal placement CLI, native-memory CLI,
  Liberty-memory CLI and the H2D CLI test.
- The H2D CLI test keeps four supplied SRAM black boxes, exports legal DEF/PNG,
  reproduces warm DEF bytes, reuses the macro frame after a standard-cell cone
  edit, and checks parameter invalidation and invalid-option errors. It caught
  and now covers generated-FIXED input classification, IO frame preservation,
  pin-order-independent fingerprints and changed-block seed overwrites.
- No real-design quality, OpenROAD timing advantage, routing completion or
  2,000-evaluation throughput claim has been established by these tests.

## Eyeriss checkpoint (2026-10-16)

- Full H2D now completes on the corrected HighTide Eyeriss: 1040 supplied SRAMs,
  193 groups, three rounds, selected round 1. Native and independent OpenROAD
  placement checks pass, as does the strict per-SRAM halo audit. Baseline and
  H2D share exact connectivity, cell masters, area, core and IO positions.
- Placement/export: none 77.0 s, H2D 592.5 s; common pin-center HPWL 6201 versus
  12565 mm. Preliminary placed-RC OpenSTA WNS improves from -77.15 to -34.13 ns;
  both miss the 2500 ps clock. No timing weighting, sizing, CTS or routing.
- Fresh OpenROAD geometry GP/DP plus previously measured preparation totals
  848.7 s, with 8941 mm HPWL and -16.84 ns WNS. Its strict per-SRAM halo audit
  fails; explicit macro-halo settings do not demonstrate clearance parity.
  The standard timing-driven HighTide reference remains incomplete.
- All 401 TopoSearch proposals are screened, with zero accepted moves. An
  independent longest-path audit of the legal selected group frame finds that
  the center x-y/x+y sequence pair needs 1805.502 um height in a 1412.64 um core.
  Conversion does not preserve this packing's outline feasibility.
- 2,388,229 named dataflow edges produce a 736 MiB topology JSON. Candidate
  copies include these immutable edges. Cells phases total 234.9 s, search
  132.1 s, abstract GP 19.9 s; another 192.0 s is outside the round timers.
  Full same-binary sampling completed separately in 720.2 s, with 64,657 main
  thread samples and byte-identical DEF. Driver topology copy/destruction is
  19.4% of samples, TopoSearch 18.9%, detailed placement 16.1% (inclusive inside
  the cells path), dataflow 2.25%, and host reporting 0.72%.
  Final violation reporting copies all dataflow for each of 193 anchors and
  848 order constraints. Remove these 1041 copies first using const references
  and a cached group index. Sampling overhead/variation is excluded from
  unprofiled quality times.
- Bring-up fixes: one-pass IO weight indexing; conservative group packing
  recovery; bounded private filler widths; exact copying of refreshed fixed
  terminals during continued GP. A small fixed-terminal reproducer failed
  before its correction. Seven final focused optimized targets pass; this
  does not replace the earlier broader synthetic validation with a new full run.
- LiveHD artifacts: `repros/eyeriss_h2d_20261016/comparison.html` and its JSON,
  CSV, matched-scale PNGs, exact commands, hashes and independent audits.

## Pending: external quality evaluation and future work

- [ ] Repair geometry-to-sequence-pair conversion and validate a feasible
  incumbent before searching; retain the screen and hard guidance. Recompute
  legal fallback slack rather than retaining a failed LP's outline slack.
- [ ] Complete the standard HighTide reference and align actual per-SRAM halo
  semantics before a controlled cross-backend QoR claim. Keep the completed
  geometry-only comparison and shared preparation explicitly labeled.
- [ ] Broaden the measured designs after fixing Eyeriss quality/performance.
  Keep explicit options/core/IOs controlled; report cells/area, HPWL, OpenSTA,
  PNGs and phase times. Stop before CTS/routing until placement justifies them.
- [ ] Sweep baseline `none` versus H2D and ablations: groups only, anchors,
  order/loop, pseudo-nets, lite/SA, optional inflation. Label synthetic robustness
  tests separately from demonstrated quality discrimination.
- [ ] Evaluate post-route TNS/DRC/congestion through the same downstream OpenROAD
  flow when placement results justify that cost. Any external RM experiment
  must record its revision/environment and import none of its source code.
- [ ] Remove full-topology copying in final constraint reporting (1041 copies,
  19.4% of main-thread samples); evaluate constraints against const references
  and one group index, preserving exact arithmetic/output.
- [ ] Measure an opt-in cells-DP-on-final-frame experiment: repeated detailed
  placement accounts for 16.1% of samples. Preserve current defaults pending
  quality/legality/halo comparisons.
- [ ] Apply the Eyeriss profile: share immutable dataflow, names and membership
  while copying only mutable group geometry/permutations. Screen sequence pairs
  before full abstraction. Aggregate solver endpoint pairs after the existing
  threshold, preserving multiplicities/weights. Serialize once and intern cache
  names rather than copying/reparsing hundreds of MiB. HHDS already supplies cached
  topological traversal via `body().nodes(hhds::Node_order::forward)`; reuse it
  in hierarchy/dataflow preparation rather than computing another order.
  Preserve contextual hierarchy, shortest register hops, fanout cuts and cycles.
  TopoSearch's flat changing geometry graph stays HHDS-independent.
- [ ] Measure macro-heavy scale, memory and LP/search throughput; the 50/5,000
  synthetic LP test is correctness evidence, not a large-design performance win.
- [ ] Future: HHDS instance/leaf guidance attributes through the existing typed
  topology/constraint seam, with explicit precedence and incremental tests.
- [ ] Future: joint macro/IO timing optimization, pin-access-aware channels,
  nonrectangular/flexible groups and intra-group pin-access re-annealing.
