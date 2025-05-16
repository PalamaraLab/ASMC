//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include "Simd.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>
#include <iostream>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "Simd.cpp"

#include <hwy/foreach_target.h>
#include <hwy/highway.h>

HWY_BEFORE_NAMESPACE();

namespace asmc::HWY_NAMESPACE
{

namespace hn = hwy::HWY_NAMESPACE;

int getNumSimdLanes_hwy()
{
  return static_cast<int>(hn::Lanes(hn::ScalableTag<float>()));
}

void printRuntimeSimdInfo_hwy()
{
  fmt::print("Detected {} float lanes and targeting {}\n", hn::Lanes(hn::ScalableTag<float>()),
             hwy::TargetName(HWY_TARGET));
}

void validateBatchSize_hwy(const int batchSize)
{
  const hn::ScalableTag<float> d;
  const size_t lanes = hn::Lanes(d); // Determine the SIMD lane count

  // Check if batchSize is divisible by the lane count
  if (batchSize <= 0 || batchSize % lanes != 0) {
    throw std::runtime_error(fmt::format(
        "Error: the batch size ({}) must be a positive multiple of the SIMD lane width ({})", batchSize, lanes));
  }
}

void warnIfSimdMismatch_hwy(const std::string& expected)
{
  const std::string actual = hwy::TargetName(HWY_TARGET);

  if (expected != actual) {
    fmt::print(std::cerr,
               "[Warning] SIMD backend mismatch: this test was designed to run with {}, but is running with {}. "
               "Numerical differences may be possible.\n",
               expected, actual);
    std::cerr.flush();
  }
}

void calculateScalingBatch_hwy(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings,
                               Eigen::Ref<Eigen::ArrayXf> sums, const int batchSize, const int numStates)
{
  const hn::ScalableTag<float> d;
  const int lanes = static_cast<int>(hn::Lanes(d));

  const auto zero = hn::Zero(d);
  for (int i = 0; i < batchSize; i += lanes) {
    hn::Store(zero, d, &sums(i));
  }

  for (int stateIdx = 0; stateIdx < numStates; ++stateIdx) {
    for (int batchItem = 0; batchItem < batchSize; batchItem += lanes) {
      auto idx = stateIdx * batchSize + batchItem;
      auto sum_vec = hn::Load(d, &sums(batchItem));
      auto vec_vec = hn::Load(d, &vec(idx));
      sum_vec = hn::Add(sum_vec, vec_vec);
      hn::Store(sum_vec, d, &sums(batchItem));
    }
  }

  for (int batchItem = 0; batchItem < batchSize; batchItem += lanes) {
    auto sum_vec = hn::Load(d, &sums(batchItem));
    auto one = hn::Set(d, 1.0f);
    auto scale_vec = hn::Div(one, sum_vec);
    hn::Store(scale_vec, d, &scalings(batchItem));
  }
}

void applyScalingBatch_hwy(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings, const int batchSize,
                           const int numStates)
{
  const hn::ScalableTag<float> d;
  const int lanes = static_cast<int>(hn::Lanes(d));

  for (int stateIdx = 0; stateIdx < numStates; ++stateIdx) {
    for (int batchItem = 0; batchItem < batchSize; batchItem += lanes) {
      auto vec_ptr = &vec(stateIdx * batchSize + batchItem);
      auto scale_ptr = &scalings(batchItem);

      auto vec_lanes = hn::Load(d, vec_ptr);
      auto scale_lanes = hn::Load(d, scale_ptr);
      auto scaled = hn::Mul(vec_lanes, scale_lanes);

      hn::Store(scaled, d, vec_ptr);
    }
  }
}

void normalizeAlphaWithBeta_hwy(Eigen::Ref<Eigen::ArrayXf> alpha, Eigen::Ref<Eigen::ArrayXf> beta,
                                Eigen::Ref<Eigen::ArrayXf> scale, int batchSize, int numStates, int from, int to)
{
  const hn::ScalableTag<float> d;
  const int lanes = static_cast<int>(hn::Lanes(d));

  // Zero out scale
  for (int pos = from; pos < to; ++pos) {
    for (int v = 0; v < batchSize; v += lanes) {
      hn::Store(hn::Zero(d), d, &scale(pos * batchSize + v));
    }
  }

  // Multiply alpha * beta, store in alpha, accumulate into scale
  for (int pos = from; pos < to; ++pos) {
    for (int k = 0; k < numStates; ++k) {
      for (int v = 0; v < batchSize; v += lanes) {
        const int ind = (pos * numStates + k) * batchSize + v;
        auto alpha_vec = hn::Load(d, &alpha(ind));
        auto beta_vec = hn::Load(d, &beta(ind));
        auto prod = hn::Mul(alpha_vec, beta_vec);
        hn::Store(prod, d, &alpha(ind));

        auto scale_vec = hn::Load(d, &scale(pos * batchSize + v));
        hn::Store(hn::Add(scale_vec, prod), d, &scale(pos * batchSize + v));
      }
    }
  }

  // scale = 1.0 / scale
  for (int pos = from; pos < to; ++pos) {
    for (int v = 0; v < batchSize; v += lanes) {
      const int ind = pos * batchSize + v;
      auto scale_vec = hn::Load(d, &scale(ind));
      hn::Store(hn::Div(hn::Set(d, 1.0f), scale_vec), d, &scale(ind));
    }
  }

  // alpha *= scale
  for (int pos = from; pos < to; ++pos) {
    for (int k = 0; k < numStates; ++k) {
      for (int v = 0; v < batchSize; v += lanes) {
        const int ind = (pos * numStates + k) * batchSize + v;
        auto alpha_vec = hn::Load(d, &alpha(ind));
        auto scale_vec = hn::Load(d, &scale(pos * batchSize + v));
        hn::Store(hn::Mul(alpha_vec, scale_vec), d, &alpha(ind));
      }
    }
  }
}

void updateAlphaColumn_hwy(Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> previousAlpha, int batchSize,
                           int k)
{
  const hn::ScalableTag<float> d;
  const int lanes = static_cast<int>(hn::Lanes(d));

  const int offset_k = k * batchSize;
  const int offset_kplus = (k + 1) * batchSize;

  for (int v = 0; v < batchSize; v += lanes) {
    auto next_alpha = hn::Load(d, &alphaC[offset_kplus + v]);
    auto prev_alpha = hn::Load(d, &previousAlpha[offset_k + v]);
    hn::Store(hn::Add(next_alpha, prev_alpha), d, &alphaC[offset_k + v]);
  }
}

} // namespace asmc::HWY_NAMESPACE

HWY_AFTER_NAMESPACE();

#if HWY_ONCE

namespace asmc
{

// Export functions with the correct signature
HWY_EXPORT(getNumSimdLanes_hwy);
HWY_EXPORT(validateBatchSize_hwy);
HWY_EXPORT(printRuntimeSimdInfo_hwy);
HWY_EXPORT(warnIfSimdMismatch_hwy);
HWY_EXPORT(calculateScalingBatch_hwy);
HWY_EXPORT(applyScalingBatch_hwy);
HWY_EXPORT(normalizeAlphaWithBeta_hwy);
HWY_EXPORT(updateAlphaColumn_hwy);

int getNumSimdLanes()
{
  return HWY_DYNAMIC_DISPATCH(getNumSimdLanes_hwy)();
}

void printRuntimeSimdInfo()
{
  return HWY_DYNAMIC_DISPATCH(printRuntimeSimdInfo_hwy)();
}

void validateBatchSize(int batchSize)
{
  return HWY_DYNAMIC_DISPATCH(validateBatchSize_hwy)(batchSize);
}

void warnIfSimdMismatch(const std::string& expected)
{
  return HWY_DYNAMIC_DISPATCH(warnIfSimdMismatch_hwy)(expected);
}

void calculateScalingBatch(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings,
                           Eigen::Ref<Eigen::ArrayXf> sums, int batchSize, int numStates)
{
  return HWY_DYNAMIC_DISPATCH(calculateScalingBatch_hwy)(vec, scalings, sums, batchSize, numStates);
}

void applyScalingBatch(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings, const int batchSize,
                       const int numStates)
{
  HWY_DYNAMIC_DISPATCH(applyScalingBatch_hwy)(vec, scalings, batchSize, numStates);
}

void normalizeAlphaWithBeta(Eigen::Ref<Eigen::ArrayXf> alpha, Eigen::Ref<Eigen::ArrayXf> beta,
                            Eigen::Ref<Eigen::ArrayXf> scale, int batchSize, int numStates, int from, int to)
{
  HWY_DYNAMIC_DISPATCH(normalizeAlphaWithBeta_hwy)(alpha, beta, scale, batchSize, numStates, from, to);
}

void updateAlphaColumn(Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> previousAlpha, int batchSize,
                       int k)
{
  HWY_DYNAMIC_DISPATCH(updateAlphaColumn_hwy)(alphaC, previousAlpha, batchSize, k);
}

} // namespace asmc

#endif // HWY_ONCE
