// BSD-3-Clause. See LICENSE.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "macro_legalize/src/sp_legalize.h"

namespace {
using DreamPlace::SequencePairNode;
bool overlap(const SequencePairNode& a, double ax, double ay, const SequencePairNode& b, double bx, double by) {
  return std::min(ax + a.w, bx + b.w) > std::max(ax, bx) && std::min(ay + a.h, by + b.h) > std::max(ay, by);
}
int check(const char* what, const std::vector<SequencePairNode>& nodes, const DreamPlace::SequencePairResult& r,
          double xh, double yh) {
  if (r.overflow != 0) {
    std::cerr << what << ": no legal packing found (overflow " << r.overflow << ")\n";
    return 1;
  }
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (r.x[i] < 0 || r.y[i] < 0 || r.x[i] + nodes[i].w > xh || r.y[i] + nodes[i].h > yh) {
      std::cerr << what << ": node " << i << " leaves the core\n";
      return 1;
    }
    if (nodes[i].pinned && (r.x[i] != nodes[i].tx || r.y[i] != nodes[i].ty)) {
      std::cerr << what << ": pinned node " << i << " moved\n";
      return 1;
    }
    for (size_t j = i + 1; j < nodes.size(); ++j)
      if (overlap(nodes[i], r.x[i], r.y[i], nodes[j], r.x[j], r.y[j])) {
        std::cerr.precision(17);
        std::cerr << what << ": nodes " << i << " and " << j << " overlap: (" << r.x[i] << ", " << r.y[i] << ", "
                  << nodes[i].w << ", " << nodes[i].h << ") (" << r.x[j] << ", " << r.y[j] << ", " << nodes[j].w
                  << ", " << nodes[j].h << ")\n";
        return 1;
      }
  }
  return 0;
}
}  // namespace

int main() {
  // Four 5x5 macros all aimed at the center of a 10x10 core: the only legal
  // pose fills the core completely (a greedy one-at-a-time search can fail).
  std::vector<SequencePairNode> full(4);
  for (auto& n : full) {
    n.w = n.h = 5;
    n.tx = n.ty = 2.5;
  }
  if (check("full 2x2", full, DreamPlace::sequencePairAnneal(full, 0, 0, 10, 10, 1, 1), 10, 10)) return 1;
  // A pinned obstacle in the lower-left corner; three 5x5 macros take the rest.
  std::vector<SequencePairNode> pinned(4);
  for (auto& n : pinned) {
    n.w = n.h = 5;
    n.tx = n.ty = 0;
  }
  pinned[0].pinned = true;
  if (check("pinned", pinned, DreamPlace::sequencePairAnneal(pinned, 0, 0, 10, 10, 1, 1), 10, 10)) return 1;
  // Mixed sizes at 84% fill (one 6x10, two 4x4, one 4x2 in 10x10), targets scrambled.
  std::vector<SequencePairNode> mixed(4);
  mixed[0].w = 6, mixed[0].h = 10, mixed[0].tx = 4, mixed[0].ty = 0;
  mixed[1].w = 4, mixed[1].h = 4, mixed[1].tx = 0, mixed[1].ty = 0;
  mixed[2].w = 4, mixed[2].h = 4, mixed[2].tx = 1, mixed[2].ty = 1;
  mixed[3].w = 4, mixed[3].h = 2, mixed[3].tx = 5, mixed[3].ty = 5;
  if (check("mixed", mixed, DreamPlace::sequencePairAnneal(mixed, 0, 0, 10, 10, 1, 1), 10, 10)) return 1;
  // HighTide litepci on asap7 at 75% utilization: 20 FakeRAMs + pcie_us, each padded by a 5 um halo per side
  // (padded macros fill 83% of the 806 um core), all aimed at the core center as a collapsed global placement.
  {
    struct Kind {
      double w, h;
      int n;
    };
    const Kind kinds[] = {{238.680, 293.761, 2}, {119.340, 73.441, 2}, {238.680, 146.881, 4},
                          {105.571, 146.881, 8}, {84.456, 146.881, 4},  {80.000, 80.000, 1}};
    // in DBU (1000/um), as the legalizer sees them: integral sites (54) and rows (270)
    const double core = 806000, site = 54, row = 270;
    std::vector<SequencePairNode> dense;
    for (const auto& k : kinds)
      for (int i = 0; i < k.n; ++i) {
        SequencePairNode m;
        m.w = std::ceil((k.w + 10) * 1000 / site) * site;
        m.h = std::ceil((k.h + 10) * 1000 / row) * row;
        m.tx = std::round((core - m.w) / 2) + 3000 * i;
        m.ty = std::round((core - m.h) / 2) - 2000 * i;
        dense.push_back(m);
      }
    if (check("litepci 75%", dense, DreamPlace::sequencePairAnneal(dense, 0, 0, core, core, site, row), core, core))
      return 1;
  }
  // Five 5x5 macros cannot fit a 10x10 core: overflow must be reported.
  std::vector<SequencePairNode> over(5);
  for (auto& n : over) n.w = n.h = 5;
  if (DreamPlace::sequencePairAnneal(over, 0, 0, 10, 10, 1, 1, 0.2).overflow <= 0) {
    std::cerr << "infeasible packing reported as legal\n";
    return 1;
  }
  return 0;
}
