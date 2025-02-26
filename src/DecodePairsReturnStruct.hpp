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

#ifndef FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
#define FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP

#include <Eigen/Core>

#include <string>
#include <tuple>
#include <vector>

struct DecodePairsReturnStruct {

private:

  bool m_storeFullPosteriors = false;
  bool m_storeSumOfPosteriors = false;
  bool m_storePerPairPosteriors = false;
  bool m_storePerPairMAPs = false;

  std::size_t numWritten = 0ul;

public:
  void initialise(const std::vector<unsigned long>& individualsA, const std::vector<unsigned long>& individualsB,
                  long int numSites, long int numStates, bool _fullPosteriors = false, bool _sumOfPosteriors = false,
                  bool _perPairPosteriors = false, bool _perPairMAPs = false)
  {
    numWritten = 0ul;
    Eigen::Index numPairsToDecode = individualsA.size();

    m_storeFullPosteriors = _fullPosteriors;
    m_storeSumOfPosteriors = _sumOfPosteriors;
    m_storePerPairPosteriors = _perPairPosteriors;
    m_storePerPairMAPs = _perPairMAPs;

    perPairIndices.resize(numPairsToDecode);

    if (m_storeFullPosteriors) {
      perPairPosteriors.resize(numPairsToDecode);
      for (auto& arr : perPairPosteriors) {
        arr.resize(numStates, numSites);
      }
    }

    if (m_storeSumOfPosteriors) {
      sumOfPosteriors.resize(numStates, numSites);
      sumOfPosteriors.setZero();
    }

    if (m_storePerPairPosteriors) {
      perPairPosteriorMeans.resize(numPairsToDecode, numSites);
      minPosteriorMeans.resize(numSites);
      argminPosteriorMeans.resize(numSites);
    }

    if (m_storePerPairMAPs) {
      perPairMAPs.resize(numPairsToDecode, numSites);
      minMAPs.resize(numSites);
      argminMAPs.resize(numSites);
    }
  }

  /// iHapIdx, iHapId, jHapIdx, jHapId
  std::vector<std::tuple<unsigned long, std::string, unsigned long, std::string>> perPairIndices;

  /// The full set of posteriors: for each pair this is a (states * numSites) matrix
  std::vector<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> perPairPosteriors;

  /// The sum of all posteriors in perPairPosteriors: a (states * numSites) matrix
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sumOfPosteriors;

  /// Posterior means: each row is an array of length numSites
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> perPairPosteriorMeans;

  Eigen::Array<float, 1, Eigen::Dynamic, Eigen::RowMajor> minPosteriorMeans;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> argminPosteriorMeans;

  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> perPairMAPs;

  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> minMAPs;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> argminMAPs;

  void incrementNumWritten()
  {
    numWritten += 1;
  }

  void finaliseCalculations()
  {
    for (Eigen::Index siteIdx = 0ll; siteIdx < perPairPosteriorMeans.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      minPosteriorMeans(siteIdx) = perPairPosteriorMeans.col(siteIdx).minCoeff(&argmin);
      argminPosteriorMeans(siteIdx) = static_cast<int>(argmin);
    }

    for (Eigen::Index siteIdx = 0ll; siteIdx < perPairMAPs.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      minMAPs(siteIdx) = perPairMAPs.col(siteIdx).minCoeff(&argmin);
      argminMAPs(siteIdx) = static_cast<int>(argmin);
    }
  }

  [[nodiscard]] std::size_t getNumWritten() const
  {
    return numWritten;
  }
};

#endif // FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
