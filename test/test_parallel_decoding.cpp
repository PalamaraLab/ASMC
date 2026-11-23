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

// ASMC::decodePairsParallel against ASMC::decodePairs. Both run the same batched recursions on the same batches, so
// the results must be identical, whatever the number of threads.

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <string>
#include <vector>

#include "ASMC.hpp"
#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"
#include "HmmUtils.hpp"

namespace
{

/// Parameters for one of the example datasets. A batch size of 16 is valid under every SIMD target, and keeps the
/// alpha and beta tables that each thread allocates for the 97,117-site sequence data under 1 GB.
DecodingParams makeParams(const std::string& inputRoot, const std::string& decodingMode)
{
  DecodingParams params(inputRoot, ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz", "", 1, 1,
                        decodingMode);
  params.batchSize = 16;
  return params;
}

/// Haplotype pairs from the 300 haplotypes of the example data. 37 pairs make two full batches of 16 and a final
/// batch of 5, which is padded up to a multiple of the SIMD lane count.
struct TestPairs {
  std::vector<unsigned long> hapsA;
  std::vector<unsigned long> hapsB;

  explicit TestPairs(const unsigned long numPairs)
  {
    for (auto pair = 0ul; pair < numPairs; ++pair) {
      hapsA.push_back(2ul * pair);
      hapsB.push_back(299ul - 3ul * pair);
    }
  }
};

/// Check that every output decodePairsParallel produces is identical to decodePairs's, for 1 to 4 threads.
/// With 1 thread, decodePairsParallel falls back to decodePairs's code.
void checkParallelMatchesSequential(ASMC::ASMC& asmc, const TestPairs& pairs, const unsigned from, const unsigned to)
{
  asmc.decodePairs(pairs.hapsA, pairs.hapsB, from, to);
  const DecodePairsReturnStruct expected = asmc.getCopyOfResults();

  for (const int numThreads : {1, 2, 3, 4}) {
    INFO(numThreads << " threads, " << pairs.hapsA.size() << " pairs, sites [" << from << ", " << to << ")");
    asmc.decodePairsParallel(pairs.hapsA, pairs.hapsB, from, to, 0.5f, numThreads);
    const DecodePairsReturnStruct actual = asmc.getCopyOfResults();

    CHECK(actual.perPairIndices == expected.perPairIndices);
    REQUIRE(actual.perPairPosteriorMeans.rows() == expected.perPairPosteriorMeans.rows());
    REQUIRE(actual.perPairPosteriorMeans.cols() == expected.perPairPosteriorMeans.cols());
    CHECK((actual.perPairPosteriorMeans == expected.perPairPosteriorMeans).all());
    CHECK((actual.minPosteriorMeans == expected.minPosteriorMeans).all());
    CHECK((actual.argminPosteriorMeans == expected.argminPosteriorMeans).all());
    REQUIRE(actual.perPairMAPs.rows() == expected.perPairMAPs.rows());
    REQUIRE(actual.perPairMAPs.cols() == expected.perPairMAPs.cols());
    CHECK((actual.perPairMAPs == expected.perPairMAPs).all());
    CHECK((actual.minMAPs == expected.minMAPs).all());
    CHECK((actual.argminMAPs == expected.argminMAPs).all());
  }
}

/**
 * A window [from, to), starting from site 1000 onwards, whose decoded sites, burn-in of 0.5 cM included, number
 * `numSegments` checkpoint segments plus `extraSites`. The burn-in moves in steps as the window grows, so both ends
 * are searched. Returns {0, 0} if there is none.
 */
std::pair<unsigned, unsigned> findWindow(const Data& data, const unsigned numSegments, const unsigned extraSites)
{
  const auto numSites = static_cast<unsigned>(data.sites);
  for (unsigned from = 1000u; from < 2000u; ++from) {
    const unsigned decodedFrom = asmc::getFromPosition(data.geneticPositions, from, 0.5f);
    for (unsigned to = from + 1u; to <= numSites; ++to) {
      const unsigned decodedSites = asmc::getToPosition(data.geneticPositions, to, 0.5f) - decodedFrom;
      if (decodedSites == numSegments * HMM::checkpointSegmentLength(decodedSites) + extraSites) {
        return {from, to};
      }
      if (decodedSites > numSegments * HMM::checkpointSegmentLength(decodedSites) + extraSites) {
        break;
      }
    }
  }
  return {0u, 0u};
}

} // namespace

TEST_CASE("parallel decoding matches decodePairs at checkpoint segment boundaries", "[ASMC][parallel]")
{
  const DecodingParams params = makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array", "array");
  const Data data(params);
  ASMC::ASMC asmc(params);
  asmc.setStorePerPairMap();

  // {segments, extra sites}: whole segments only, one site into a new segment, a partial last segment, and fewer
  // sites than one segment
  const std::vector<std::pair<unsigned, unsigned>> layouts = {{2u, 0u}, {2u, 1u}, {1u, 255u}, {0u, 100u}};
  for (const auto& [numSegments, extraSites] : layouts) {
    INFO(numSegments << " segments and " << extraSites << " sites");
    const auto [from, to] = findWindow(data, numSegments, extraSites);
    REQUIRE(to > 0u);
    checkParallelMatchesSequential(asmc, TestPairs(37ul), from, to);
  }
}

TEST_CASE("parallel decoding matches decodePairs in array mode", "[ASMC][parallel]")
{
  ASMC::ASMC asmc(makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array", "array"));
  asmc.setStorePerPairMap();

  SECTION("all sites")
  {
    checkParallelMatchesSequential(asmc, TestPairs(37ul), 0u, 0u);
  }
  SECTION("a window")
  {
    checkParallelMatchesSequential(asmc, TestPairs(37ul), 1000u, 3000u);
  }
  SECTION("fewer pairs than SIMD lanes")
  {
    checkParallelMatchesSequential(asmc, TestPairs(1ul), 0u, 0u);
  }
}

TEST_CASE("parallel decoding matches decodePairs in sequence mode", "[ASMC][parallel]")
{
  ASMC::ASMC asmc(makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300", "sequence"));
  asmc.setStorePerPairMap();
  checkParallelMatchesSequential(asmc, TestPairs(37ul), 5000u, 7000u);
}

TEST_CASE("parallel decoding rejects outputs it does not produce", "[ASMC][parallel]")
{
  ASMC::ASMC asmc(makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array", "array"));
  const TestPairs pairs(3ul);

  // Checked for every thread count, so the behaviour does not depend on whether OpenMP is available
  SECTION("per-pair posteriors")
  {
    asmc.setStorePerPairPosterior();
    for (const int numThreads : {1, 2}) {
      CHECK_THROWS_AS(asmc.decodePairsParallel(pairs.hapsA, pairs.hapsB, 0u, 0u, 0.5f, numThreads),
                      std::invalid_argument);
    }
  }
  SECTION("the sum of posteriors")
  {
    asmc.setStoreSumOfPosterior();
    for (const int numThreads : {1, 2}) {
      CHECK_THROWS_AS(asmc.decodePairsParallel(pairs.hapsA, pairs.hapsB, 0u, 0u, 0.5f, numThreads),
                      std::invalid_argument);
    }
  }
}
