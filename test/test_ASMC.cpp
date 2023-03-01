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

#include "catch.hpp"

#include <utility>
#include <vector>

#include "ASMC.hpp"

#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/ostream.h>

TEST_CASE("test ASMC decodePairsArray", "[ASMC]")
{
  ASMC::ASMC asmc(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                  ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  asmc.setStorePerPairMap();
  asmc.setStorePerPairPosterior();
  asmc.setStorePerPairPosteriorMean();
  asmc.setStorePerPairMap();

  std::vector<unsigned long> indA = {1, 2, 3};
  std::vector<unsigned long> indB = {2, 3, 4};
  asmc.decodePairs(indA, indB);
  auto result = asmc.getRefOfResults();

  SECTION("test decode pair summarize")
  {
    REQUIRE(result.perPairIndices.size() == 3ul);

    // 0.1% margin in this test as the results can vary between pure and avx/sse
    REQUIRE(result.perPairPosteriorMeans(0, 0) == Approx(15968.91016f).margin(15968.91016f * 0.001f));
    REQUIRE(result.perPairPosteriorMeans(1, 8) == Approx(27963.49805f).margin(27963.49805f * 0.001f));
    REQUIRE(result.perPairPosteriorMeans(2, 29) == Approx(48573.32812f).margin(48573.32812f * 0.001f));

    REQUIRE(result.perPairMAPs(0, 0) == 29);
    REQUIRE(result.perPairMAPs(1, 1234) == 65);
    REQUIRE(result.perPairMAPs(2, 7) == 33);

    // Check that the posteriors actually sum to one
    for (Eigen::Index idx = 0ll; idx < result.perPairPosteriors.size(); ++idx) {
      REQUIRE(result.perPairPosteriors.at(idx).colwise().sum().isOnes(1e-2));
    }
  }
}

TEST_CASE("test ASMC decodePairsSequence", "[ASMC]")
{
  ASMC::ASMC asmc(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300",
                  ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz", "", "sequence");

  asmc.setStorePerPairMap();
  asmc.setStorePerPairPosterior();
  asmc.setStorePerPairPosteriorMean();
  asmc.setStorePerPairMap();

  std::vector<unsigned long> indA = {5, 6};
  std::vector<unsigned long> indB = {7, 8};
  asmc.decodePairs(indA, indB);
  auto result = asmc.getRefOfResults();

  SECTION("test decode pair summarize")
  {
    REQUIRE(result.perPairIndices.size() == 2ul);

    // 0.1% margin in this test as the results can vary between pure and avx/sse
    REQUIRE(result.perPairPosteriorMeans(0, 0) == Approx(801.06647f).margin(801.06647f * 0.001f));
    REQUIRE(result.perPairPosteriorMeans(1, 8) == Approx(17953.60938f).margin(17953.60938f * 0.001f));

    REQUIRE(result.perPairMAPs(0, 0) == 16);
    REQUIRE(result.perPairMAPs(1, 1234) == 61);

    // Check that the posteriors actually sum to one
    for (Eigen::Index idx = 0ll; idx < result.perPairPosteriors.size(); ++idx) {
      REQUIRE(result.perPairPosteriors.at(idx).colwise().sum().isOnes(1e-2));
    }
  }
}

TEST_CASE("test other get methods", "[ASMC]")
{
  ASMC::ASMC asmc(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                  ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  const std::vector<float>& expectedTimes = asmc.getExpectedTimes();
  CHECK(expectedTimes.at(0) == Approx(14.999777896567f).margin(1e-5));
  CHECK(expectedTimes.at(4) == Approx(135.698150766900f).margin(1e-5));
}

TEST_CASE("test from and to", "[ASMC]")
{
  ASMC::ASMC asmc_full(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                       ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  ASMC::ASMC asmc_part(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                       ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  std::vector<unsigned long> indA = {1, 2, 3, 4, 5};
  std::vector<unsigned long> indB = {6, 7, 8, 9, 10};

  asmc_full.setStorePerPairMap();
  asmc_full.setStorePerPairPosterior();
  asmc_full.setStorePerPairPosteriorMean();
  asmc_full.setStoreSumOfPosterior();

  asmc_part.setStorePerPairMap();
  asmc_part.setStorePerPairPosterior();
  asmc_part.setStorePerPairPosteriorMean();
  asmc_part.setStoreSumOfPosterior();

  asmc_full.decodePairs(indA, indB);
  auto result_full = asmc_full.getRefOfResults();
  asmc_part.setStorePerPairMap();

  const unsigned lo = 1000;
  const unsigned hi = 1100;
  const unsigned long windowSize = static_cast<unsigned long>(hi - lo);
  asmc_part.decodePairs(indA, indB, lo, hi, 0.5f);
  auto result_part = asmc_part.getRefOfResults();

  SECTION("test part sizes are correct")
  {
    REQUIRE(result_part.perPairPosteriors.front().rows() == result_full.perPairPosteriors.front().rows());
    REQUIRE(result_part.perPairPosteriors.front().cols() == windowSize);

    REQUIRE(result_part.sumOfPosteriors.rows() == result_full.sumOfPosteriors.rows());
    REQUIRE(result_part.sumOfPosteriors.cols() == windowSize);

    REQUIRE(result_part.perPairPosteriorMeans.rows() == result_full.perPairPosteriorMeans.rows());
    REQUIRE(result_part.perPairPosteriorMeans.cols() == windowSize);

    REQUIRE(result_part.minPosteriorMeans.cols() == windowSize);
    REQUIRE(result_part.argminPosteriorMeans.cols() == windowSize);

    REQUIRE(result_part.perPairMAPs.rows() == result_full.perPairMAPs.rows());
    REQUIRE(result_part.perPairMAPs.cols() == windowSize);

    REQUIRE(result_part.minMAPs.cols() == windowSize);
    REQUIRE(result_part.argminMAPs.cols() == windowSize);
  }

  SECTION("test parts match full analysis")
  {
    for(auto i = 0ul; i < indA.size(); ++i) {
        REQUIRE(
          (result_full.perPairPosteriors.at(i).middleCols(static_cast<unsigned long>(lo), windowSize)
            - result_part.perPairPosteriors.at(i)).abs().maxCoeff() < 1e-6
        );
    }
  }
}
