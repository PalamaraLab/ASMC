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

// ASMC regression tests: each decodes a few pairs, or one small job, under one decoding configuration and compares the
// output with snapshots. See regression_snapshots.hpp for how to regenerate the snapshots.

#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "ASMC.hpp"
#include "DecodingParams.hpp"
#include "regression_snapshots.hpp"

namespace
{

using asmc::test::checkMatrixMatchesSnapshot;
using asmc::test::exactMatch;
using asmc::test::floatOutput;

constexpr const char* arrayRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
constexpr const char* sequenceRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300";
constexpr const char* decodingQuantities = ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";

/// First haplotype of each pair decoded by the regression tests: two pairs within an individual, two across.
std::vector<unsigned long> pairHapsA()
{
  return {0ul, 2ul, 17ul, 150ul};
}

/// Second haplotype of each pair; see pairHapsA.
std::vector<unsigned long> pairHapsB()
{
  return {1ul, 5ul, 230ul, 151ul};
}

/**
 * ASMC parameters for the example data at `inputRoot`.
 *
 * @param decodingMode "array" or "sequence"
 * @param compress use the classic emission model (no CSFS)
 * @param useAncestral treat allele 1 as derived instead of folding to the minor allele
 */
DecodingParams asmcParams(const std::string& inputRoot, const std::string& decodingMode, const bool compress = false,
                          const bool useAncestral = false)
{
  return DecodingParams(inputRoot, decodingQuantities, "", 1, 1, decodingMode, false, true, compress, useAncestral,
                        std::numeric_limits<float>::quiet_NaN());
}

/**
 * Decode the regression pairs over sites [from, to) and compare their posterior means and MAPs with snapshots
 * `<name>.perPairPosteriorMeans.gz` and `<name>.perPairMAPs.gz`. `to` = 0 means the last site.
 */
void checkDecodedPairsMatchSnapshots(DecodingParams params, const std::string& name, const unsigned from = 0u,
                                     const unsigned to = 0u)
{
  ASMC::ASMC asmc(std::move(params));
  asmc.setStorePerPairPosteriorMean(true);
  asmc.setStorePerPairMap(true);
  asmc.decodePairs(pairHapsA(), pairHapsB(), from, to, 0.5f);

  const DecodePairsReturnStruct& results = asmc.getRefOfResults();
  checkMatrixMatchesSnapshot(results.perPairPosteriorMeans, name + ".perPairPosteriorMeans.gz", floatOutput);
  checkMatrixMatchesSnapshot(results.perPairMAPs, name + ".perPairMAPs.gz", exactMatch);
}

} // namespace

TEST_CASE("ASMC regression: array data", "[ASMC_regression]")
{
  checkDecodedPairsMatchSnapshots(asmcParams(arrayRoot, "array"), "asmc_array");
}

TEST_CASE("ASMC regression: sequence data in a window", "[ASMC_regression]")
{
  // The sequence data has about 97,000 sites: a window keeps this fast. A small batch keeps the alpha and beta
  // buffers (sites x states x batch size) small; the batch size does not change the results.
  DecodingParams params = asmcParams(sequenceRoot, "sequence");
  params.batchSize = 16;
  checkDecodedPairsMatchSnapshots(params, "asmc_sequence_window", 5000u, 7000u);
}

TEST_CASE("ASMC regression: compressed emissions", "[ASMC_regression]")
{
  checkDecodedPairsMatchSnapshots(asmcParams(arrayRoot, "array", /*compress=*/true), "asmc_compress");
}

TEST_CASE("ASMC regression: unfolded data with ancestral alleles", "[ASMC_regression]")
{
  checkDecodedPairsMatchSnapshots(asmcParams(arrayRoot, "array", false, /*useAncestral=*/true), "asmc_ancestral");
}

TEST_CASE("ASMC regression: posterior sums for one job of decodeAll", "[ASMC_regression]")
{
  // Job 8 of 1000 decodes about 45 of the 44,850 pairs, in one padded batch
  const DecodingParams params(arrayRoot, decodingQuantities, "", /*_jobs=*/1000, /*_jobInd=*/8, "array", false, true,
                              false, false, std::numeric_limits<float>::quiet_NaN(), false,
                              /*_doPosteriorSums=*/true, false, "", false, /*_doMajorMinorPosteriorSums=*/true);
  ASMC::ASMC asmc(params);
  // The ASMC constructors turn on per-pair posterior-mean storage, which decodeAllInJob cannot fill: it throws
  // std::out_of_range unless storage is turned off first (#57)
  asmc.setStorePerPairPosteriorMean(false);
  const DecodingReturnValues sums = asmc.decodeAllInJob();

  // The sums are sites x states; every 50th site keeps the snapshots small
  const Eigen::Index siteStride = 50;
  checkMatrixMatchesSnapshot(asmc::test::everyNthRow(sums.sumOverPairs, siteStride),
                             "asmc_decode_all_job.sumOverPairs.gz", floatOutput);
  checkMatrixMatchesSnapshot(asmc::test::everyNthRow(sums.sumOverPairs00, siteStride),
                             "asmc_decode_all_job.sumOverPairs00.gz", floatOutput);
  checkMatrixMatchesSnapshot(asmc::test::everyNthRow(sums.sumOverPairs01, siteStride),
                             "asmc_decode_all_job.sumOverPairs01.gz", floatOutput);
  checkMatrixMatchesSnapshot(asmc::test::everyNthRow(sums.sumOverPairs11, siteStride),
                             "asmc_decode_all_job.sumOverPairs11.gz", floatOutput);
}
