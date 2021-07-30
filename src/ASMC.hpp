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

#ifndef ASMC_HPP
#define ASMC_HPP

#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"

#include <string>
#include <vector>

namespace ASMC
{

class ASMC
{

private:
  DecodingParams mParams;
  Data mData;
  HMM mHmm;

public:
  /**
   * ASMC constructor with full control over parameters, by manually specifying a DecodingParams object.
   *
   * @param params the decoding parameters
   */
  explicit ASMC(DecodingParams params);

  /**
   * ASMC constructor that will set sensible defaults. If you wish to fine-tune parameters, use the constructor that
   * takes a DecodingParams object, which you can configure manually.
   *
   * @param inFileRoot the input file root
   * @param decodingQuantFile the decoding quantities file
   * @param outFileRoot the output file root, default to the input file root
   */
  ASMC(const std::string& inFileRoot, const std::string& decodingQuantFile, const std::string& outFileRoot = "");

  DecodingReturnValues decodeAllInJob();

  void decodePairs(const std::vector<unsigned long>& hapIndicesA, const std::vector<unsigned long>& hapIndicesB,
                   bool perPairPosteriors = false, bool sumOfPosteriors = false, bool perPairPosteriorMeans = false,
                   bool perPairMAPs = false);

  void decodePairs(const std::vector<std::string>& hapIdsA, const std::vector<std::string>& hapIdsB,
                   bool perPairPosteriors = false, bool sumOfPosteriors = false, bool perPairPosteriorMeans = false,
                   bool perPairMAPs = false);

  DecodePairsReturnStruct getCopyOfResults();

  const DecodePairsReturnStruct& getRefOfResults();
};
} // namespace ASMC

#endif
