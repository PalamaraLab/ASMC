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

#include "ASMC.hpp"
#include "HmmUtils.hpp"

#include <fmt/format.h>

#include <exception>
#include <iostream>

ASMC::ASMC::ASMC(DecodingParams params) : mParams{std::move(params)}, mData{mParams}, mHmm{mData, mParams}
{
}

ASMC::ASMC::ASMC(const std::string& inFileRoot, const std::string& decodingQuantFile, const std::string& outFileRoot)
    : mParams{inFileRoot,
              decodingQuantFile,
              outFileRoot.empty() ? inFileRoot : outFileRoot,
              1,
              1,
              "array",
              false,
              true,
              false,
              false,
              0.f,
              false,
              true,
              false,
              "",
              false,
              true,
              true},
      mData{mParams}, mHmm{mData, mParams}
{
}

DecodingReturnValues ASMC::ASMC::decodeAllInJob()
{
  std::cout << "Decoding job " << mParams.jobInd << " of " << mParams.jobs << "\n\n";

  std::cout << "Will decode " << mParams.decodingModeString << " data." << std::endl;
  std::cout << "Output will have prefix: " << mParams.outFileRoot << std::endl;

  if (mParams.compress) {
    std::cout << "Will use classic emission model (no CSFS)." << std::endl;
  } else {
    std::cout << "Minimum marker distance to use CSFS is set to " << mParams.skipCSFSdistance << "." << std::endl;
  }

  if (mParams.useAncestral) {
    std::cout << "Assuming ancestral alleles are correctly encoded." << std::endl;
  }

  if (mParams.doPosteriorSums) {
    std::cout << "Will output sum of posterior tables for all pairs." << std::endl;
  }

  if (mParams.doMajorMinorPosteriorSums) {
    std::cout << "Will output sum of posterior tables for all pairs, partitioned by major/minor alleles." << std::endl;
  }

  mHmm.decodeAll(mParams.jobs, mParams.jobInd);
  return mHmm.getDecodingReturnValues();
}

void ASMC::ASMC::decodePairs(const std::vector<unsigned long>& hapIndicesA,
                             const std::vector<unsigned long>& hapIndicesB, bool perPairPosteriors,
                             bool sumOfPosteriors, bool perPairPosteriorMeans, bool perPairMAPs)
{
  if (hapIndicesA.empty() || hapIndicesA.size() != hapIndicesB.size()) {
    throw std::runtime_error(
        fmt::format("Vector of A indices ({}) must be the same size as vector of B indices ({}).\n", hapIndicesA.size(),
                    hapIndicesB.size()));
  }

  mHmm.getDecodePairsReturnStruct().initialise(hapIndicesA, hapIndicesB, mData.sites,
                                               mHmm.getDecodingQuantities().states, perPairPosteriors, sumOfPosteriors,
                                               perPairPosteriorMeans, perPairMAPs);
  mHmm.setStorePerPairPosteriorMean(perPairPosteriorMeans);
  mHmm.setStorePerPairMap(perPairMAPs);
  mHmm.setStorePerPairPosterior(perPairPosteriors);
  mHmm.setStoreSumOfPosterior(sumOfPosteriors);
  mHmm.decodeHapPairs(hapIndicesA, hapIndicesB);
  mHmm.finishDecoding();
  mHmm.getDecodePairsReturnStruct().finaliseCalculations();
}

void ASMC::ASMC::decodePairs(const std::vector<std::string>& hapIdsA, const std::vector<std::string>& hapIdsB,
                             bool perPairPosteriors, bool sumOfPosteriors, bool perPairPosteriorMeans, bool perPairMAPs)
{
  if (hapIdsA.size() != hapIdsB.size()) {
    throw std::runtime_error(fmt::format("Vector of A IDs ({}) must be the same size as vector of B IDs ({}).\n",
                                         hapIdsA.size(), hapIdsB.size()));
  }

  std::vector<unsigned long> hapsIndicesA;
  std::vector<unsigned long> hapsIndicesB;

  hapsIndicesA.resize(hapIdsA.size());
  hapsIndicesB.resize(hapIdsB.size());

  for (auto i = 0ul; i < hapIdsA.size(); ++i) {
    auto [strIdA, hapA] = asmc::combinedIdToIndPlusHap(hapIdsA.at(i));
    auto [strIdB, hapB] = asmc::combinedIdToIndPlusHap(hapIdsB.at(i));

    unsigned long indIdA = asmc::getIndIdxFromIdString(mData.IIDList, strIdA);
    unsigned long indIdB = asmc::getIndIdxFromIdString(mData.IIDList, strIdB);

    hapsIndicesA.at(i) = asmc::dipToHapId(indIdA, hapA);
    hapsIndicesB.at(i) = asmc::dipToHapId(indIdB, hapB);
  }

  decodePairs(hapsIndicesA, hapsIndicesB, perPairPosteriors, sumOfPosteriors, perPairPosteriorMeans, perPairMAPs);
}

DecodePairsReturnStruct ASMC::ASMC::getCopyOfResults()
{
  return mHmm.getDecodePairsReturnStruct();
}

const DecodePairsReturnStruct& ASMC::ASMC::getRefOfResults()
{
  return mHmm.getDecodePairsReturnStruct();
}
