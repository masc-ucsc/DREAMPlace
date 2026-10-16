// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/seed.hpp"
#include "dreamplace/native/macro/toposearch.hpp"
namespace DreamPlace::native::macro {
struct Iteration {
  std::map<std::string, double> phases;
  double gp_overflow = 0;
  bool gp_converged = false;
  std::vector<DifferenceLP::Slack> slack;
  std::vector<std::string> relaxed;
  int screened = 0;
  int round = 0, gp_iterations = 0, lp_solves = 0, evaluations = 0,
      accepted_moves = 0;
  double hpwl = 0, rudy_overflow = 0, merit = 0, anchor_deviation = 0,
         order_violation = 0, seconds = 0;
  bool legal = false, selected = false;
};
struct ConstraintViolation {
  std::string id, kind;
  double amount = 0;
};
struct Result {
  DreamPlace::native::Result placement;
  std::vector<ConstraintViolation> violations;
  Topology topology;
  Parameters resolved;
  std::vector<Iteration> iterations;
  double hpwl = 0, rudy_overflow = 0, merit = 0;
  bool legal = false;
};
// Strong exception guarantee, fixed inputs preserved. No file I/O, JSON, HHDS,
// or OpenROAD dependency. Host writes output/cache artifacts separately.
// cells.observer receives one expanded physical global_placement snapshot per
// outer round, before group legalization; other internal solver frames stay private.
Result place_h2d(Design &, std::span<const NodeInfo>,
                 const Parameters &parameters = {}, const Topology &input = {},
                 const Options &cells = {});
} // namespace DreamPlace::native::macro
