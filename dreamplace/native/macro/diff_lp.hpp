// BSD-3-Clause. See LICENSE.
#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
namespace DreamPlace::native::macro {
struct Infeasible : std::runtime_error {
  std::vector<std::string> cycle;
  explicit Infeasible(std::vector<std::string> ids);
};
// Integer pitch-unit primal with weighted linear/L1 terms and difference
// constraints. A positive penalty makes an arc soft; zero means hard.
// u or v == -1 denotes the fixed zero origin. No generic LP or GLPK backend.
struct DifferenceLP {
  using Int = std::int64_t;
  struct Variable {
    Int lower, upper;
    double coefficient = 0;
    std::string id;
  };
  struct Arc {
    int u, v;
    Int bound;
    double penalty = 0;
    std::string id;
  };
  struct L1 {
    int variable;
    Int anchor;
    double weight;
  };
  struct Slack {
    std::string id;
    Int amount;
    double weight;
  };
  struct Solution {
    std::vector<Int> x;
    std::vector<Slack> slack;
    double objective = 0, slack_cost = 0;
  };
  std::vector<Variable> variables;
  std::vector<Arc> arcs;
  std::vector<L1> displacement;
  int variable(Int lo, Int hi, double coefficient = 0, std::string id = {});
  void constrain(int u, int v, Int bound, double penalty = 0,
                 std::string id = {});
  void anchor(int variable, Int position, double weight);
  Solution solve() const;
};
} // namespace DreamPlace::native::macro
