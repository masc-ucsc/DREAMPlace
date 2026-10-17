// BSD-3-Clause. See LICENSE.
#pragma once
#include <string>
#include <vector>
namespace DreamPlace::native {
// Centers are independent targets, not hard regions for their member cells.
struct FloorplanCircle {
  std::string name;
  double x = 0, y = 0, radius = 1, area = 1;
  int family = 0, parent = -1;
};
struct FloorplanOptions {
  double alignment_tolerance = 1; // center offset / smaller radius; strict <
  double spacing_tolerance = .25; // maximum range / mean gap in a class
  double spacing_penalty = .02;  // normalized squared gap error + penalty/class
  double symmetry_tolerance = 1; // midpoint offset / smaller radius
  double max_displacement = 1;   // per-axis local movement / own radius
};
struct FloorplanRelation {
  std::string kind;       // alignment, spacing, symmetry
  int axis = 0;           // 0: x coordinate, 1: y coordinate
  std::vector<int> nodes; // spacing: consecutive (a,b) edge endpoints
  double value = 0;       // fitted line, gap, or symmetry-axis coordinate
};
struct FloorplanSnap {
  std::vector<FloorplanCircle> circles;
  std::vector<FloorplanRelation> relations;
  int attempted = 0, rejected = 0, infeasible_families = 0;
  double seconds = 0, maximum_overlap = 0, maximum_displacement = 0;
  bool feasible = true;
};
// Greedy local relation selection and weighted least-squares fitting. Circle
// separation uses conservative halfspaces along the input pair directions.
// Alignment, spacing and symmetry relations require the same family AND parent.
// Root circles (parent -1) share the implicit top scope. Hard same-family
// separation still applies across parents; parent/child translations are
// inherited. Failed families retain their input centers and are explicitly
// reported. Core bounds apply to centers; circles are interaction ranges, not
// fences.
FloorplanSnap snap_floorplan(const std::vector<FloorplanCircle> &circles,
                             double xl, double yl, double xh, double yh,
                             const FloorplanOptions &options = {});
} // namespace DreamPlace::native
