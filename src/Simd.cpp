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
  fflush(stdout);
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

void updateAlphaColumn_hwy(Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> previousAlpha, int batchSize, int numStates)
{
  for (int k = numStates - 2; k >= 0; k--) {

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
}

void updateAlphaForwardStep_hwy(Eigen::Ref<Eigen::ArrayXf> nextAlpha, Eigen::Ref<Eigen::ArrayXf> previousAlpha,
                            Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> AU, const float* B,
                            const float* U, const float* D, const std::vector<float>& columnRatios,
                            const std::vector<float>& emission1AtSite, const std::vector<float>& emission0minus1AtSite,
                            const std::vector<float>& emission2minus0AtSite, Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch,
                            Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch, int batchSize, int numStates, int pos)
{
  const hn::ScalableTag<float> d;
  const int lanes = static_cast<int>(hn::Lanes(d));

  for (int k = 0; k < numStates; k++) {
    const int offset_k = k * batchSize;
    const int offset_kplus1 = (k + 1) * batchSize;
    const int offset_kminus1 = (k - 1) * batchSize;
    const int obs_offset = pos * batchSize;

    const auto D_k = hn::Set(d, D[k]);
    const auto B_k = (k < numStates - 1) ? hn::Set(d, B[k]) : hn::Zero(d);
    const auto em1 = hn::Set(d, emission1AtSite[k]);
    const auto em0m1 = hn::Set(d, emission0minus1AtSite[k]);
    const auto em2m0 = hn::Set(d, emission2minus0AtSite[k]);

    auto U_km1 = k > 0 ? hn::Set(d, U[k - 1]) : hn::Zero(d);
    auto col_km1 = k > 0 ? hn::Set(d, columnRatios[k - 1]) : hn::Zero(d);

    for (int v = 0; v < batchSize; v += lanes) {

      auto AU_v = hn::Load(d, &AU[v]);

      if (k > 0) {
        auto prev_km1 = hn::Load(d, &previousAlpha[offset_kminus1 + v]);
        auto term1 = hn::Mul(U_km1, prev_km1);
        auto term2 = hn::Mul(col_km1, AU_v);
        AU_v = hn::Add(term1, term2);
        hn::Store(AU_v, d, &AU[v]);
      }

      auto prev_k = hn::Load(d, &previousAlpha[offset_k + v]);
      auto term = hn::Add(AU_v, hn::Mul(D_k, prev_k));

      if (k < numStates - 1) {
        auto ac_kplus1 = hn::Load(d, &alphaC[offset_kplus1 + v]);
        term = hn::Add(term, hn::Mul(B_k, ac_kplus1));
      }

      auto obs0 = hn::Load(d, &obsIsZeroBatch[obs_offset + v]);
      auto obs2 = hn::Load(d, &obsIsTwoBatch[obs_offset + v]);

      auto emission_term = hn::Add(em1, hn::Mul(em0m1, obs0));
      emission_term = hn::Add(emission_term, hn::Mul(em2m0, obs2));

      hn::Store(hn::Mul(emission_term, term), d, &nextAlpha[offset_k + v]);
    }
  }
}

void computeBetaEmissionProduct_hwy(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> lastComputedBeta,
                                    Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                                    const std::vector<float>& emission1AtSite,
                                    const std::vector<float>& emission0minus1AtSite,
                                    const std::vector<float>& emission2minus0AtSite, int batchSize, int numStates,
                                    int pos)
{
  namespace hn = hwy::HWY_NAMESPACE;
  const hn::ScalableTag<float> d;
  const int lanes = hn::Lanes(d);
  const int obsOffset = (pos + 1) * batchSize;

  for (int k = 0; k < numStates; ++k) {
    const auto em1 = hn::Set(d, emission1AtSite[k]);
    const auto em0m1 = hn::Set(d, emission0minus1AtSite[k]);
    const auto em2m0 = hn::Set(d, emission2minus0AtSite[k]);

    const int offset = k * batchSize;

    for (int v = 0; v < batchSize; v += lanes) {
      auto obs0 = hn::Load(d, &obsIsZeroBatch[obsOffset + v]);
      auto obs2 = hn::Load(d, &obsIsTwoBatch[obsOffset + v]);

      auto emission_k = hn::Add(em1, hn::Mul(em0m1, obs0));
      emission_k = hn::Add(emission_k, hn::Mul(em2m0, obs2));

      auto beta = hn::Load(d, &lastComputedBeta[offset + v]);
      auto product = hn::Mul(emission_k, beta);

      hn::Store(product, d, &vec[offset + v]);
    }
  }
}

void computeBetaUpwardSweep_hwy(Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> vec,
                                const std::vector<float>& U, const std::vector<float>& RR, int batchSize, int numStates)
{
  namespace hn = hwy::HWY_NAMESPACE;
  const hn::ScalableTag<float> d;
  const int lanes = hn::Lanes(d);

  for (int k = numStates - 2; k >= 0; --k) {
    const auto U_k = hn::Set(d, U[k]);
    const auto RR_k = hn::Set(d, RR[k]);

    const int offset_k = k * batchSize;
    const int offset_kplus1 = (k + 1) * batchSize;

    for (int v = 0; v < batchSize; v += lanes) {
      const auto vec_next = hn::Load(d, &vec[offset_kplus1 + v]);
      const auto bu_next = hn::Load(d, &BU[offset_kplus1 + v]);

      const auto term1 = hn::Mul(U_k, vec_next);
      const auto term2 = hn::Mul(RR_k, bu_next);

      hn::Store(hn::Add(term1, term2), d, &BU[offset_k + v]);
    }
  }
}

void computeBetaFinalCombine_hwy(Eigen::Ref<Eigen::ArrayXf> currentBeta, Eigen::Ref<Eigen::ArrayXf> BL,
                                 Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> vec,
                                 const std::vector<float>& B, const std::vector<float>& D, int batchSize, int numStates)
{
  namespace hn = hwy::HWY_NAMESPACE;
  const hn::ScalableTag<float> d;
  const int lanes = hn::Lanes(d);

  for (int k = 0; k < numStates; ++k) {
    const auto D_k = hn::Set(d, D[k]);
    hn::Vec<decltype(d)> B_km1;

    if (k > 0) {
      B_km1 = hn::Set(d, B[k - 1]);
    }

    const int offset_k = k * batchSize;
    const int offset_kminus1 = (k - 1) * batchSize;

    for (int v = 0; v < batchSize; v += lanes) {
      auto BL_v = hn::Load(d, &BL[v]);

      if (k > 0) {
        auto vec_km1 = hn::Load(d, &vec[offset_kminus1 + v]);
        BL_v = hn::Add(BL_v, hn::Mul(B_km1, vec_km1));
        hn::Store(BL_v, d, &BL[v]);
      }

      const auto vec_k = hn::Load(d, &vec[offset_k + v]);
      const auto BU_k = hn::Load(d, &BU[offset_k + v]);
      const auto D_term = hn::Mul(D_k, vec_k);

      auto total = hn::Add(BL_v, hn::Add(D_term, BU_k));
      hn::Store(total, d, &currentBeta[offset_k + v]);
    }
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
HWY_EXPORT(updateAlphaForwardStep_hwy);
HWY_EXPORT(computeBetaEmissionProduct_hwy);
HWY_EXPORT(computeBetaUpwardSweep_hwy);
HWY_EXPORT(computeBetaFinalCombine_hwy);

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

void updateAlphaColumn(Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> previousAlpha, int batchSize, int numStates)
{
  HWY_DYNAMIC_DISPATCH(updateAlphaColumn_hwy)(alphaC, previousAlpha, batchSize, numStates);
}

void updateAlphaForwardStep(Eigen::Ref<Eigen::ArrayXf> nextAlpha, Eigen::Ref<Eigen::ArrayXf> previousAlpha,
                            Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> AU, const float* B,
                            const float* U, const float* D, const std::vector<float>& columnRatios,
                            const std::vector<float>& emission1AtSite, const std::vector<float>& emission0minus1AtSite,
                            const std::vector<float>& emission2minus0AtSite, Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch,
                            Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch, int batchSize, int numStates, int pos)
{
  HWY_DYNAMIC_DISPATCH(updateAlphaForwardStep_hwy)
  (nextAlpha, previousAlpha, alphaC, AU, B, U, D, columnRatios, emission1AtSite, emission0minus1AtSite,
   emission2minus0AtSite, obsIsZeroBatch, obsIsTwoBatch, batchSize, numStates, pos);
}

void computeBetaEmissionProduct(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> lastComputedBeta,
                                Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                                const std::vector<float>& emission1AtSite,
                                const std::vector<float>& emission0minus1AtSite,
                                const std::vector<float>& emission2minus0AtSite, int batchSize, int numStates, int pos)
{
  HWY_DYNAMIC_DISPATCH(computeBetaEmissionProduct_hwy)
  (vec, lastComputedBeta, obsIsZeroBatch, obsIsTwoBatch, emission1AtSite, emission0minus1AtSite, emission2minus0AtSite,
   batchSize, numStates, pos);
}

void computeBetaUpwardSweep(Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> vec, const std::vector<float>& U,
                            const std::vector<float>& RR, int batchSize, int numStates)
{
  HWY_DYNAMIC_DISPATCH(computeBetaUpwardSweep_hwy)(BU, vec, U, RR, batchSize, numStates);
}

void computeBetaFinalCombine(Eigen::Ref<Eigen::ArrayXf> currentBeta, Eigen::Ref<Eigen::ArrayXf> BL,
                             Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> vec, const std::vector<float>& B,
                             const std::vector<float>& D, int batchSize, int numStates)
{
  HWY_DYNAMIC_DISPATCH(computeBetaFinalCombine_hwy)(currentBeta, BL, BU, vec, B, D, batchSize, numStates);
}

} // namespace asmc

#endif // HWY_ONCE
