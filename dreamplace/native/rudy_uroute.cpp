// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/uroute_weights.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>

#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#endif

namespace DreamPlace::native {
namespace {
struct Tensor {
  int c, x, y;
  std::vector<float> data;
  Tensor(int channels, int nx, int ny)
      : c(channels), x(nx), y(ny), data(size_t(c) * x * y) {}
  size_t area() const { return size_t(x) * y; }
};
#if defined(__aarch64__) && defined(__ARM_NEON)
constexpr bool scalar_gemm_order = false;
#else
constexpr bool scalar_gemm_order = true;
#endif
Tensor conv_spatial(const UrouteWeights &w, const std::string &name,
                    const Tensor &input, int outputs, int kernel) {
  const std::array<int, 4> shape{outputs, input.c, kernel, kernel};
  const std::array<int, 1> bias_shape{outputs};
  const auto weights = w.tensor(name + ".weight", shape);
  const auto bias = w.tensor(name + ".bias", bias_shape);
  Tensor result(outputs, input.x, input.y);
  const auto area = input.area();
  const int pad = kernel / 2;
  // The scalar backend matches the f32 GEMM boundary: explicit FMA products,
  // then bias. Compiler-dependent contraction drifted at cancellation-sensitive
  // outputs. Keep the existing NEON backend's scalar-head rounding unchanged.
  // The contiguous loop needs no parallel reduction or external BLAS.
  const auto accumulate = [](float weight, float value, float sum) {
    if constexpr (scalar_gemm_order)
      return std::fma(weight, value, sum);
    else
      return sum + weight * value;
  };
  int oc = 0;
  // Reuse each input value across four independent output planes. The
  // ic/kx/ky accumulation order within every output is unchanged; there is
  // no parallel reduction or approximation. Single-channel heads use the
  // scalar remainder below.
  for (; oc + 3 < outputs; oc += 4) {
    auto *out0 = result.data.data() + size_t(oc) * area;
    auto *out1 = out0 + area;
    auto *out2 = out1 + area;
    auto *out3 = out2 + area;
    std::fill(out0, out0 + area, scalar_gemm_order ? 0.f : bias[oc]);
    std::fill(out1, out1 + area, scalar_gemm_order ? 0.f : bias[oc + 1]);
    std::fill(out2, out2 + area, scalar_gemm_order ? 0.f : bias[oc + 2]);
    std::fill(out3, out3 + area, scalar_gemm_order ? 0.f : bias[oc + 3]);
    const auto stride = size_t(input.c) * kernel * kernel;
    // Keep four small output tiles resident while accumulating channels.
    // Every output still visits ic/kx/ky in the original order.
    for (int xb = 0; xb < input.x; xb += 8) {
      for (int ic = 0; ic < input.c; ++ic) {
        const auto *in = input.data.data() + size_t(ic) * area;
        for (int kx = 0; kx < kernel; ++kx) {
          const int sx = kx - pad, x0 = std::max(xb, -sx),
                    x1 = std::min({input.x, input.x - sx, xb + 8});
          for (int ky = 0; ky < kernel; ++ky) {
            const int sy = ky - pad, y0 = std::max(0, -sy),
                      y1 = std::min(input.y, input.y - sy);
            const auto index =
                size_t(oc) * stride + (size_t(ic) * kernel + kx) * kernel + ky;
            const float w0 = weights[index], w1 = weights[index + stride],
                        w2 = weights[index + 2 * stride],
                        w3 = weights[index + 3 * stride];
            for (int x = x0; x < x1; ++x) {
              const auto *row = in + size_t(x + sx) * input.y;
              auto *d0 = out0 + size_t(x) * input.y;
              auto *d1 = out1 + size_t(x) * input.y;
              auto *d2 = out2 + size_t(x) * input.y;
              auto *d3 = out3 + size_t(x) * input.y;
              for (int y = y0; y < y1; ++y) {
                const float value = row[y + sy];
                d0[y] = accumulate(w0, value, d0[y]);
                d1[y] = accumulate(w1, value, d1[y]);
                d2[y] = accumulate(w2, value, d2[y]);
                d3[y] = accumulate(w3, value, d3[y]);
              }
            }
          }
        }
      }
    }
  }
  for (; oc < outputs; ++oc) {
    auto *out = result.data.data() + oc * area;
    std::fill(out, out + area, scalar_gemm_order ? 0.f : bias[oc]);
    for (int ic = 0; ic < input.c; ++ic) {
      const auto *in = input.data.data() + ic * area;
      for (int kx = 0; kx < kernel; ++kx) {
        const int sx = kx - pad, x0 = std::max(0, -sx),
                  x1 = std::min(input.x, input.x - sx);
        for (int ky = 0; ky < kernel; ++ky) {
          const int sy = ky - pad, y0 = std::max(0, -sy),
                    y1 = std::min(input.y, input.y - sy);
          const float weight =
              weights[((oc * input.c + ic) * kernel + kx) * kernel + ky];
          for (int x = x0; x < x1; ++x) {
            const auto *row = in + size_t(x + sx) * input.y;
            auto *dest = out + size_t(x) * input.y;
            for (int y = y0; y < y1; ++y)
              dest[y] = accumulate(weight, row[y + sy], dest[y]);
          }
        }
      }
    }
  }
  if constexpr (scalar_gemm_order)
    for (int c = 0; c < outputs; ++c)
      for (size_t i = 0; i < area; ++i)
        result.data[size_t(c) * area + i] += bias[c];
  for (float v : result.data)
    if (!std::isfinite(v))
      throw std::runtime_error("URoute convolution overflow");
  return result;
}
// A bounded K-by-16 input tile is reused across all output planes. Each
// SIMD lane is an independent pixel and visits ic/kx/ky in order; there is
// no parallel reduction. Zero padding matches the network's convolution
// contract. ReLU/sigmoid also remove any intermediate zero-sign distinction.
// The largest v1 tile is 72 KiB, rather than a full-grid im2col allocation.
Tensor conv(const UrouteWeights &w, const std::string &name,
            const Tensor &input, int outputs, int kernel,
            const Tensor *seed = nullptr, int weight_channels = 0,
            int channel_offset = 0) {
#if defined(__aarch64__) && defined(__ARM_NEON)

  if (outputs < 4 || outputs % 4) {
    // The spatial fallback has no seeded/partial-channel form.
    if (seed || weight_channels || channel_offset)
      throw std::invalid_argument(
          "seeded URoute convolution needs 4-way outputs");
    return conv_spatial(w, name, input, outputs, kernel);
  }
  using Pack = float32x4_t;
  const std::array<int, 4> shape{
      outputs, weight_channels ? weight_channels : input.c, kernel, kernel};
  const std::array<int, 1> bs{outputs};
  const auto weights = w.tensor(name + ".weight", shape);
  const auto bias = w.tensor(name + ".bias", bs);
  Tensor result(outputs, input.x, input.y);
  const auto area = input.area();
  const int pad = kernel / 2;
  constexpr int tile = 16;
  const int stride = input.c * kernel * kernel;
  const int full_stride =
      (weight_channels ? weight_channels : input.c) * kernel * kernel;
  std::vector<float> packed(size_t(stride) * tile);
  for (int x = 0; x < input.x; ++x)
    for (int y = 0; y < input.y; y += tile) {
      const int count = std::min(tile, input.y - y);
      for (int ic = 0; ic < input.c; ++ic)
        for (int kx = 0; kx < kernel; ++kx)
          for (int ky = 0; ky < kernel; ++ky) {
            const int k = (ic * kernel + kx) * kernel + ky, ix = x + kx - pad,
                      iy = y + ky - pad;
            auto *dest = packed.data() + size_t(k) * tile;
            // Interior tiles overwrite every lane: avoid zero filling and
            // use a fixed-size copy. Boundary tiles retain exact zero padding.
            if (ix >= 0 && ix < input.x && count == tile && iy >= 0 &&
                iy + tile <= input.y) {
              const auto *src = input.data.data() + size_t(ic) * area +
                                size_t(ix) * input.y + size_t(iy);
              std::memcpy(dest, src, size_t(tile) * sizeof(float));
            } else {
              std::fill_n(dest, tile, 0.f);
              const int lo = std::max(0, -iy),
                        hi = std::min(count, input.y - iy);
              if (ix >= 0 && ix < input.x && hi > lo) {
                const auto *src = input.data.data() + size_t(ic) * area +
                                  size_t(ix) * input.y + size_t(iy + lo);
                std::memcpy(dest + lo, src, size_t(hi - lo) * sizeof(float));
              }
            }
          }
      for (int oc = 0; oc < outputs; oc += 4) {
        auto initial = [&](int output, int spatial) {
          Pack value = vdupq_n_f32(seed ? 0.f : bias[output]);
          if (seed && count > spatial) {
            const auto offset =
                size_t(output) * area + size_t(x) * input.y + y + spatial;
            std::memcpy(&value, seed->data.data() + offset,
                        size_t(std::min(4, count - spatial)) * sizeof(float));
          }
          return value;
        };
        Pack a00 = initial(oc + 0, 0);
        Pack a01 = initial(oc + 0, 4);
        Pack a02 = initial(oc + 0, 8);
        Pack a03 = initial(oc + 0, 12);
        Pack a10 = initial(oc + 1, 0);
        Pack a11 = initial(oc + 1, 4);
        Pack a12 = initial(oc + 1, 8);
        Pack a13 = initial(oc + 1, 12);
        Pack a20 = initial(oc + 2, 0);
        Pack a21 = initial(oc + 2, 4);
        Pack a22 = initial(oc + 2, 8);
        Pack a23 = initial(oc + 2, 12);
        Pack a30 = initial(oc + 3, 0);
        Pack a31 = initial(oc + 3, 4);
        Pack a32 = initial(oc + 3, 8);
        Pack a33 = initial(oc + 3, 12);
        for (int k = 0; k < stride; ++k) {
          const size_t index =
              size_t(oc) * full_stride + channel_offset * kernel * kernel + k;
          const float w0 = weights[index], w1 = weights[index + full_stride],
                      w2 = weights[index + 2 * full_stride],
                      w3 = weights[index + 3 * full_stride];

          Pack v0;
          std::memcpy(&v0, packed.data() + size_t(k) * tile + 0, sizeof(Pack));
          a00 = vfmaq_n_f32(a00, v0, w0);
          a10 = vfmaq_n_f32(a10, v0, w1);
          a20 = vfmaq_n_f32(a20, v0, w2);
          a30 = vfmaq_n_f32(a30, v0, w3);
          Pack v1;
          std::memcpy(&v1, packed.data() + size_t(k) * tile + 4, sizeof(Pack));
          a01 = vfmaq_n_f32(a01, v1, w0);
          a11 = vfmaq_n_f32(a11, v1, w1);
          a21 = vfmaq_n_f32(a21, v1, w2);
          a31 = vfmaq_n_f32(a31, v1, w3);
          Pack v2;
          std::memcpy(&v2, packed.data() + size_t(k) * tile + 8, sizeof(Pack));
          a02 = vfmaq_n_f32(a02, v2, w0);
          a12 = vfmaq_n_f32(a12, v2, w1);
          a22 = vfmaq_n_f32(a22, v2, w2);
          a32 = vfmaq_n_f32(a32, v2, w3);
          Pack v3;
          std::memcpy(&v3, packed.data() + size_t(k) * tile + 12, sizeof(Pack));
          a03 = vfmaq_n_f32(a03, v3, w0);
          a13 = vfmaq_n_f32(a13, v3, w1);
          a23 = vfmaq_n_f32(a23, v3, w2);
          a33 = vfmaq_n_f32(a33, v3, w3);
        }
        const size_t base = size_t(oc) * area + size_t(x) * input.y + y;
        if (count == tile) {
          std::memcpy(result.data.data() + base + 0 * area + 0, &a00,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 0 * area + 4, &a01,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 0 * area + 8, &a02,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 0 * area + 12, &a03,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 1 * area + 0, &a10,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 1 * area + 4, &a11,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 1 * area + 8, &a12,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 1 * area + 12, &a13,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 2 * area + 0, &a20,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 2 * area + 4, &a21,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 2 * area + 8, &a22,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 2 * area + 12, &a23,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 3 * area + 0, &a30,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 3 * area + 4, &a31,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 3 * area + 8, &a32,
                      sizeof(Pack));
          std::memcpy(result.data.data() + base + 3 * area + 12, &a33,
                      sizeof(Pack));
        } else {
          if (count > 0)
            std::memcpy(result.data.data() + base + 0 * area + 0, &a00,
                        size_t(std::min(4, count - 0)) * sizeof(float));
          if (count > 4)
            std::memcpy(result.data.data() + base + 0 * area + 4, &a01,
                        size_t(std::min(4, count - 4)) * sizeof(float));
          if (count > 8)
            std::memcpy(result.data.data() + base + 0 * area + 8, &a02,
                        size_t(std::min(4, count - 8)) * sizeof(float));
          if (count > 12)
            std::memcpy(result.data.data() + base + 0 * area + 12, &a03,
                        size_t(std::min(4, count - 12)) * sizeof(float));
          if (count > 0)
            std::memcpy(result.data.data() + base + 1 * area + 0, &a10,
                        size_t(std::min(4, count - 0)) * sizeof(float));
          if (count > 4)
            std::memcpy(result.data.data() + base + 1 * area + 4, &a11,
                        size_t(std::min(4, count - 4)) * sizeof(float));
          if (count > 8)
            std::memcpy(result.data.data() + base + 1 * area + 8, &a12,
                        size_t(std::min(4, count - 8)) * sizeof(float));
          if (count > 12)
            std::memcpy(result.data.data() + base + 1 * area + 12, &a13,
                        size_t(std::min(4, count - 12)) * sizeof(float));
          if (count > 0)
            std::memcpy(result.data.data() + base + 2 * area + 0, &a20,
                        size_t(std::min(4, count - 0)) * sizeof(float));
          if (count > 4)
            std::memcpy(result.data.data() + base + 2 * area + 4, &a21,
                        size_t(std::min(4, count - 4)) * sizeof(float));
          if (count > 8)
            std::memcpy(result.data.data() + base + 2 * area + 8, &a22,
                        size_t(std::min(4, count - 8)) * sizeof(float));
          if (count > 12)
            std::memcpy(result.data.data() + base + 2 * area + 12, &a23,
                        size_t(std::min(4, count - 12)) * sizeof(float));
          if (count > 0)
            std::memcpy(result.data.data() + base + 3 * area + 0, &a30,
                        size_t(std::min(4, count - 0)) * sizeof(float));
          if (count > 4)
            std::memcpy(result.data.data() + base + 3 * area + 4, &a31,
                        size_t(std::min(4, count - 4)) * sizeof(float));
          if (count > 8)
            std::memcpy(result.data.data() + base + 3 * area + 8, &a32,
                        size_t(std::min(4, count - 8)) * sizeof(float));
          if (count > 12)
            std::memcpy(result.data.data() + base + 3 * area + 12, &a33,
                        size_t(std::min(4, count - 12)) * sizeof(float));
        }
      }
    }
  for (float value : result.data)
    if (!std::isfinite(value))
      throw std::runtime_error("URoute convolution overflow");
  return result;
#else
  if (seed || weight_channels || channel_offset)
    throw std::invalid_argument("seeded URoute convolution requires NEON");
  return conv_spatial(w, name, input, outputs, kernel);
#endif
}
void relu(Tensor &t) {
  for (auto &v : t.data)
    v = std::max(0.f, v);
}
Tensor residual(const UrouteWeights &w, const std::string &name,
                const Tensor &input) {
  auto first = conv(w, name + ".conv1", input, 64, 3);
  relu(first);
  auto output = conv(w, name + ".conv2", first, 64, 3);
  auto skip = input.c == 64 ? input : conv(w, name + ".skip", input, 64, 1);
  for (size_t i = 0; i < output.data.size(); ++i)
    output.data[i] += skip.data[i];
  relu(output);
  return output;
}
Tensor resize(const Tensor &in, int nx, int ny) {
  Tensor out(in.c, nx, ny);
  const float scale_x = float(in.x) / nx, scale_y = float(in.y) / ny;
  for (int c = 0; c < in.c; ++c) {
    for (int x = 0; x < nx; ++x) {
      const float u = std::max(0.f, scale_x * (x + .5f) - .5f);
      const int x0 = int(u), x1 = std::min(x0 + 1, in.x - 1);
      const float wx = u - x0;
      for (int y = 0; y < ny; ++y) {
        const float v = std::max(0.f, scale_y * (y + .5f) - .5f);
        const int y0 = int(v), y1 = std::min(y0 + 1, in.y - 1);
        const float wy = v - y0;
        const auto base = size_t(c) * in.area();
        const float a = in.data[base + size_t(x0) * in.y + y0];
        const float b = in.data[base + size_t(x0) * in.y + y1];
        const float d = in.data[base + size_t(x1) * in.y + y0];
        const float e = in.data[base + size_t(x1) * in.y + y1];
        float value;
        // Match PyTorch CPU bilinear evaluation: small outputs use four
        // weighted samples; larger outputs interpolate each axis in turn.
        // Rounding order matters when the final head nearly cancels.
        if (nx + ny <= 128) {
          const float wa = (1 - wx) * (1 - wy), wb = (1 - wx) * wy;
          const float wd = wx * (1 - wy), we = wx * wy;
          const float pa = wa * a, pb = wb * b, pd = wd * d, pe = we * e;
          const float tail = pd + pe;
          const float middle = pb + tail;
          value = pa + middle;
        } else {
          const float first = std::fma(1 - wy, a, wy * b);
          const float second = std::fma(1 - wy, d, wy * e);
          value = std::fma(1 - wx, first, wx * second);
        }
        out.data[size_t(c) * out.area() + size_t(x) * ny + y] = value;
      }
    }
  }
  return out;
}
Tensor concat(const Tensor &a, const Tensor &b) {
  if (a.x != b.x || a.y != b.y)
    throw std::invalid_argument("URoute fusion shape mismatch");
  Tensor result(a.c + b.c, a.x, a.y);
  std::copy(a.data.begin(), a.data.end(), result.data.begin());
  std::copy(b.data.begin(), b.data.end(), result.data.begin() + a.data.size());
  return result;
}
// Scalar f32 exponential for attention, matching the CPU Torch reference's
// SLEEF polynomial and explicit FMA order. System expf can differ by an ulp,
// which is amplified by the cancellation-sensitive decoder/head.
// Adapted from SLEEF xexpf at 5a1d179df9cf652951b59010a2d2075372d67f68,
// src/libm/sleefsimdsp.c (the SLEEF revision pinned by reference Torch).
// Copyright Naoki Shibata and contributors 2010 - 2024.
/*
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/
float attention_exp(float d) {
  if (d < -104.f)
    return 0.f;
  if (d > 100.f)
    return std::numeric_limits<float>::infinity();
  const int q =
      int(std::nearbyint(d * 1.442695040888963407359924681001892137f));
  float s = std::fma(float(q), -0.693145751953125f, d);
  s = std::fma(float(q), -1.428606765330187045e-6f, s);
  float u = 0.000198527617612853646278381f;
  u = std::fma(u, s, 0.00139304355252534151077271f);
  u = std::fma(u, s, 0.00833336077630519866943359f);
  u = std::fma(u, s, 0.0416664853692054748535156f);
  u = std::fma(u, s, 0.166666671633720397949219f);
  u = std::fma(u, s, 0.5f);
  u = 1.f + std::fma(s * s, u, s);
  return std::ldexp(u, q);
}
Tensor matching(const UrouteWeights &w, const std::string &name,
                const Tensor &query, const std::vector<Tensor> &support) {
  Tensor total(64, query.x, query.y);
  const auto area = query.area();
#if defined(__aarch64__) && defined(__ARM_NEON)
  // Concatenation visits all query channels before any support channel. Reuse
  // that exact f32 prefix, then continue each support's original FMA sequence.
  // Do not split off and add an independently reduced support convolution:
  // that would change rounding at cancellation-sensitive matching gates.
  // The fixed .z weight shape already forces one support channel count.
  std::optional<Tensor> prefix;
  const int joined = support.empty() ? 0 : query.c + support.front().c;
  if (support.size() > 1)
    prefix.emplace(conv(w, name + ".z", query, 64, 3, nullptr, joined, 0));
#endif
  for (const auto &s : support) {
#if defined(__aarch64__) && defined(__ARM_NEON)
    auto z = prefix ? conv(w, name + ".z", s, 64, 3, &*prefix, joined, query.c)
                    : conv(w, name + ".z", concat(query, s), 64, 3);
#else
    auto z = conv(w, name + ".z", concat(query, s), 64, 3);
#endif
    relu(z);
    Tensor mean(1, z.x, z.y);
    // CPU channel means sum contiguous blocks of 16 before combining them.
    // A single 64-channel running sum drifts at cancellation-sensitive heads.
    for (size_t i = 0; i < area; ++i) {
      float sum = 0;
      for (int base = 0; base < 64; base += 16) {
        float chunk = 0;
        for (int c = base; c < base + 16; ++c)
          chunk += z.data[c * area + i];
        sum += chunk;
      }
      mean.data[i] = sum;
    }
    for (auto &v : mean.data)
      v /= 64;
    auto gate = conv(w, name + ".gate", mean, 1, 7);
    for (auto &v : gate.data)
      v = 1.f / (1.f + attention_exp(-v));
    for (int c = 0; c < 64; ++c)
      for (size_t i = 0; i < area; ++i) {
        // PyTorch materializes weight * sample before adding each support.
        // Keep that f32 rounding boundary; a fused multiply-add accumulates
        // differently across shots and can fail cancellation-sensitive heads.
        // Volatile enforces the boundary even when contraction across source
        // statements is enabled; it does not change convolution accumulation.
        const volatile float contribution = gate.data[i] * s.data[c * area + i];
        total.data[c * area + i] += contribution;
      }
  }
  auto out = conv(w, name + ".out", total, 64, 1);
  relu(out);
  return out;
}
Tensor input(const UrouteMaps &m, int channels) {
  if (m.channels != channels || m.bins_x < 8 || m.bins_y < 8 ||
      m.bins_x > 512 || m.bins_y > 512 ||
      m.values.size() != size_t(channels) * m.bins_x * m.bins_y ||
      !m.normalized)
    throw std::invalid_argument(
        "invalid URoute inference input shape/normalization");
  for (float v : m.values)
    if (!std::isfinite(v) || v < 0 || v > 1)
      throw std::invalid_argument("invalid URoute inference input value");
  Tensor t(channels, m.bins_x, m.bins_y);
  t.data = m.values;
  return t;
}
} // namespace
struct UrouteNetwork::Data {
  explicit Data(const UrouteWeights &w) : weights(w) {}
  UrouteWeights weights;
  int bins_x = 0, bins_y = 0;
  std::array<std::vector<Tensor>, 4> levels;
};
UrouteNetwork::UrouteNetwork(const UrouteWeights &w,
                             const std::vector<UrouteMaps> &pairs) {
  if (pairs.empty() || pairs.size() > 16)
    throw std::invalid_argument("URoute needs 1..16 support pairs");
  auto data = std::make_shared<Data>(w);
  data->bins_x = pairs.front().bins_x;
  data->bins_y = pairs.front().bins_y;
  std::vector<Tensor> support;
  for (const auto &pair : pairs) {
    if (pair.bins_x != data->bins_x || pair.bins_y != data->bins_y)
      throw std::invalid_argument("URoute support grid mismatch");
    auto s = conv(w, "s_embed", input(pair, 5), 64, 3);
    relu(s);
    support.push_back(std::move(s));
  }
  for (int level = 0; level < 4; ++level) {
    for (auto &s : support)
      s = residual(w, "support." + std::to_string(level), s);
    data->levels[level] = std::move(support);
    // A moved-from vector is valid but not guaranteed empty. Do not depend on
    // its implementation-specific state before creating the next level.
    support.clear();
    if (level != 3)
      for (const auto &s : data->levels[level])
        support.push_back(resize(s, s.x / 2, s.y / 2));
  }
  data_ = std::move(data);
}
std::vector<float> UrouteNetwork::infer(const UrouteMaps &query) const {
  if (query.bins_x != data_->bins_x || query.bins_y != data_->bins_y)
    throw std::invalid_argument("URoute support/query grid mismatch");
  const auto &w = data_->weights;
  auto q = conv(w, "q_embed", input(query, 4), 64, 3);
  relu(q);
  std::vector<Tensor> matched;
  for (int level = 0; level < 4; ++level) {
    const auto index = "." + std::to_string(level);
    q = residual(w, "query" + index, q);
    matched.push_back(matching(w, "matches" + index, q, data_->levels[level]));
    if (level != 3)
      q = resize(q, q.x / 2, q.y / 2);
  }
  auto out = std::move(q);
  for (int level = 3; level >= 0; --level) {
    const auto &m = matched[level];
    if (out.x != m.x || out.y != m.y)
      out = resize(out, m.x, m.y);
    out = residual(w, "decoder." + std::to_string(level), concat(m, out));
  }
  out = conv(w, "output", out, 1, 1);
  relu(out);
  for (float v : out.data)
    if (!std::isfinite(v))
      throw std::runtime_error("URoute inference overflow");
  return out.data;
}
std::vector<float> uroute_infer(const UrouteWeights &w, const UrouteMaps &query,
                                const std::vector<UrouteMaps> &pairs) {
  return UrouteNetwork(w, pairs).infer(query);
}
} // namespace DreamPlace::native
