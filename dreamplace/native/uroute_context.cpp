// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/rudy.hpp"
#include "dreamplace/native/uroute_weights.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace DreamPlace::native {
namespace {
void validate_map(const UrouteMaps &m, int channels, bool normalized, int x,
                  int y) {
  if (m.channels != channels || m.normalized != normalized || m.bins_x != x ||
      m.bins_y != y || m.values.size() != size_t(channels) * x * y ||
      m.ranges.size() != size_t(channels))
    throw std::invalid_argument("URoute support shape/normalization mismatch");
  const auto area = size_t(x) * y;
  for (int c = 0; c < channels; ++c) {
    const auto [lo, hi] = m.ranges[c];
    if (!std::isfinite(lo) || !std::isfinite(hi) || hi < lo || lo < 0 ||
        !std::isfinite(hi - lo))
      throw std::invalid_argument("invalid URoute support calibration");
    for (size_t i = c * area; i < (c + 1) * area; ++i) {
      const auto v = m.values[i];
      if (!std::isfinite(v) || v < 0 ||
          (normalized && (v > 1 || (lo == hi && v != 0))) ||
          (!normalized &&
           (v < static_cast<float>(lo) || v > static_cast<float>(hi))))
        throw std::invalid_argument("invalid URoute support map value");
    }
  }
}
double median(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  const auto mid = values.size() / 2;
  return values.size() % 2 ? values[mid]
                           : values[mid - 1] / 2 + values[mid] / 2;
}
} // namespace
UrouteContext::UrouteContext(std::shared_ptr<const UrouteWeights> weights,
                             std::vector<UrouteSupport> samples, double level)
    : weights_(std::move(weights)), drc_level_(level) {
  if (!weights_ || samples.empty() || samples.size() > 16 ||
      !std::isfinite(level) || level < 0)
    throw std::invalid_argument("URoute requires weights, 1..16 samples and a "
                                "finite nonnegative DRC level");
  bins_x_ = samples[0].features.bins_x;
  bins_y_ = samples[0].features.bins_y;
  if (bins_x_ < 8 || bins_y_ < 8 || bins_x_ > 512 || bins_y_ > 512)
    throw std::invalid_argument("invalid URoute support grid");
  std::array<std::vector<double>, 4> ranges;
  std::vector<UrouteMaps> utilization, drc;
  const auto area = size_t(bins_x_) * bins_y_;
  for (const auto &sample : samples) {
    if (sample.id.empty() ||
        std::find(ids_.begin(), ids_.end(), sample.id) != ids_.end())
      throw std::invalid_argument("empty/duplicate URoute support id");
    ids_.push_back(sample.id);
    validate_map(sample.features, 4, true, bins_x_, bins_y_);
    validate_map(sample.labels, 3, false, bins_x_, bins_y_);
    for (int metric = 0; metric < 2; ++metric) {
      auto pair = sample.features;
      pair.channels = 5;
      const auto [lo, hi] = sample.labels.ranges[metric];
      ranges[2 * metric].push_back(lo);
      ranges[2 * metric + 1].push_back(hi);
      pair.ranges.push_back({lo, hi});
      for (size_t i = 0; i < area; ++i) {
        const double raw = sample.labels.values[metric * area + i];
        pair.values.push_back(hi == lo ? 0.f
                                       : static_cast<float>(std::clamp(
                                             (raw - lo) / (hi - lo), 0., 1.)));
      }
      (metric == 0 ? utilization : drc).push_back(std::move(pair));
    }
  }
  utilization_range_ = {median(ranges[0]), median(ranges[1])};
  drc_range_ = {median(ranges[2]), median(ranges[3])};
  // Constant calibration needs no network. In particular, clean DRC maps must
  // not pay support-encoding time or retain an unused encoder pyramid.
  if (utilization_range_[0] != utilization_range_[1])
    utilization_ =
        std::make_shared<const UrouteNetwork>(*weights_, utilization);
  if (drc_range_[0] != drc_range_[1])
    drc_ = std::make_shared<const UrouteNetwork>(*weights_, drc);
}
UroutePrediction UrouteContext::predict(const UrouteMaps &query) const {
  validate_map(query, 4, true, bins_x_, bins_y_);
  UroutePrediction result{bins_x_, bins_y_, {}, {}};
  auto predict = [&](const auto &network, const auto &range) {
    // Constant support calibration fixes the physical answer independently of
    // normalized prediction. Clean DRC support maps commonly take this path.
    if (range[0] == range[1])
      return std::vector<double>(size_t(query.bins_x) * query.bins_y, range[0]);
    std::vector<double> values;
    const auto normalized = network->infer(query);
    for (float v : normalized) {
      const double value = range[0] + double(v) * (range[1] - range[0]);
      if (!std::isfinite(value))
        throw std::runtime_error("URoute support calibration overflow");
      values.push_back(value);
    }
    return values;
  };
  result.utilization = predict(utilization_, utilization_range_);
  result.drc_density = predict(drc_, drc_range_);
  return result;
}

RudyResult uroute_rudy(const Design &design, const RudyOptions &options) {
  auto base_options = options;
  base_options.model = RudyModel::rsmt;
  auto result = rudy(design, base_options);
  if (!options.uroute) {
    result.uroute_skipped = "no usable URoute context";
    return result;
  }
  const auto &context = *options.uroute;
  const auto query = uroute_features(design, options.bins_x, options.bins_y,
                                     context.bins_x(), context.bins_y());
  const auto prediction = context.predict(query);
  const auto utilization =
      uroute_resample(prediction.utilization, prediction.bins_x,
                      prediction.bins_y, options.bins_x, options.bins_y);
  result.drc_density =
      uroute_resample(prediction.drc_density, prediction.bins_x,
                      prediction.bins_y, options.bins_x, options.bins_y);
  const auto count = utilization.size();
  result.overflow =
      double(result.blocked_pins + result.unroutable_nets) / count;
  result.max_utilization =
      result.blocked_pins || result.unroutable_nets ? 2 : 0;
  for (size_t i = 0; i < count; ++i) {
    const double h = result.horizontal[i], v = result.vertical[i],
                 peak = std::max(h, v);
    // Preserve the directional ratio AND max(h,v)==predicted utilization.
    // Dividing by h+v would incorrectly halve equal-direction predictions.
    result.horizontal[i] = utilization[i] * (peak > 0 ? h / peak : 1);
    result.vertical[i] = utilization[i] * (peak > 0 ? v / peak : 1);
    const double u = std::max(utilization[i], result.pin_utilization[i]);
    result.max_utilization = std::max(result.max_utilization, u);
    result.overflow += std::max(0., u - 1) / count;
    if (context.drc_level() > 0)
      result.overflow +=
          std::max(0., result.drc_density[i] - context.drc_level()) / count;
  }
  if (!std::isfinite(result.overflow) || !std::isfinite(result.max_utilization))
    throw std::runtime_error("URoute calibrated overflow is nonfinite");
  return result;
}
} // namespace DreamPlace::native
