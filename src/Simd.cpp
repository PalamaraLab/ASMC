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

// Implementation of the validation function with batchSize as a parameter
void validateBatchSize_hwy(const int batchSize)
{
  const hn::ScalableTag<float> d;
  const size_t lanes = hn::Lanes(d); // Determine the SIMD lane count

  // Check if batchSize is divisible by the lane count
  if (batchSize % lanes != 0) {
    throw std::invalid_argument(
        fmt::format("Error: the batch size ({}) must be a multiple of the SIMD lane width ({})", batchSize, lanes));
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

} // namespace asmc

#endif // HWY_ONCE
