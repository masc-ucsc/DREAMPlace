// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/uroute.hpp"
#include <memory>
#include <string>

namespace DreamPlace::native {
// The mapping is read-only and remains valid across atomic manifest/weight
// replacement. Schema and all floating-point values are checked once on load.
class UrouteWeights {
public:
  explicit UrouteWeights(const std::filesystem::path &);
  std::span<const float> tensor(const std::string &,
                                std::span<const int> shape) const;

private:
  struct Data;
  std::shared_ptr<const Data> data_;
};
// Network-only forward pass, useful for exporter golden tests. Query has four
// normalized channels, support pairs have five (features + normalized label).
// Returns one nonnegative prediction per query tile. No weight updates.
std::vector<float> uroute_infer(const UrouteWeights &, const UrouteMaps &query,
                                const std::vector<UrouteMaps> &support);
// Immutable, query-independent support encoder pyramid. Construct once for a
// fixed metric/support set; infer keeps all query scratch local to the call.
class UrouteNetwork {
public:
  UrouteNetwork(const UrouteWeights &, const std::vector<UrouteMaps> &support);
  std::vector<float> infer(const UrouteMaps &query) const;

private:
  struct Data;
  std::shared_ptr<const Data> data_;
};
struct UrouteSupport {
  std::string id;
  UrouteMaps features, labels; // labels: raw utilization, DRC density, detour
};
struct UroutePrediction {
  int bins_x = 0, bins_y = 0;
  std::vector<double> utilization, drc_density;
};
// Immutable stage snapshot. The host owns manifest/version/hash verification
// and deterministic support selection. Constructor validates all sample maps.
class UrouteContext {
public:
  UrouteContext(std::shared_ptr<const UrouteWeights>,
                std::vector<UrouteSupport>, double drc_level = 0);
  UroutePrediction predict(const UrouteMaps &) const;
  const std::vector<std::string> &sample_ids() const { return ids_; }
  int bins_x() const { return bins_x_; }
  int bins_y() const { return bins_y_; }
  double drc_level() const { return drc_level_; }

private:
  std::shared_ptr<const UrouteWeights> weights_;
  std::shared_ptr<const UrouteNetwork> utilization_, drc_;
  std::vector<std::string> ids_;
  std::array<double, 2> utilization_range_, drc_range_;
  int bins_x_ = 0, bins_y_ = 0;
  double drc_level_ = 0;
};
} // namespace DreamPlace::native
