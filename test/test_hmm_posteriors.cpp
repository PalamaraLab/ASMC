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

// Numerical checks on HMM posteriors: the batched SIMD path against the non-batched scalar path (HMM::decode),
// under every SIMD target available at runtime, plus batching, windowing, normalisation and summary statistics.

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <hwy/targets.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"
#include "HmmUtils.hpp"
#include "test_fixtures.hpp"

namespace
{

using asmc::test::decodeBatched;
using asmc::test::decodeReference;
using asmc::test::maxAbsPosteriorDifference;
using asmc::test::ReferencePosterior;
using asmc::test::RowMajorPosterior;

// The batched path runs the same recursions as HMM::decode, vectorised across pairs, so the two should differ only by
// floating-point rounding (for example, SIMD kernels may use fused multiply-add). Posteriors lie in [0, 1].
constexpr float posteriorTolerance = 1e-5f;

// A batch size that is a multiple of the lane count of every SIMD target, so one HMM can be decoded under each.
constexpr int batchSizeForAllTargets = 16;

// Burn-in either side of a decoding window, in centimorgans. getFromPosition asserts that this is positive.
constexpr float cmBurnIn = 0.5f;

/// First haplotype of each test pair from exampleFile.n300 (300 haplotypes). With pairHapsB: two pairs within one
/// individual and three across individuals. A function rather than a global, so construction cannot throw before main.
std::vector<unsigned long> pairHapsA()
{
  return {0ul, 2ul, 17ul, 150ul, 298ul};
}

/// Second haplotype of each test pair; see pairHapsA.
std::vector<unsigned long> pairHapsB()
{
  return {1ul, 5ul, 230ul, 151ul, 299ul};
}

/// Parameters for one of the example datasets, with a batch size valid under every SIMD target.
DecodingParams makeParams(const std::string& inputRoot, const std::string& decodingMode)
{
  DecodingParams params(inputRoot, ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz", "", 1, 1,
                        decodingMode);
  params.batchSize = batchSizeForAllTargets;
  return params;
}

DecodingParams arrayParams()
{
  return makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array", "array");
}

DecodingParams sequenceParams()
{
  return makeParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300", "sequence");
}

/// Restricts Highway's dynamic dispatch to one SIMD target while this object exists, and restores the default after,
/// even if an assertion throws.
class ScopedSimdTarget
{
public:
  explicit ScopedSimdTarget(const int64_t target)
  {
    hwy::SetSupportedTargetsForTest(target);
  }

  ~ScopedSimdTarget()
  {
    hwy::SetSupportedTargetsForTest(0);
  }

  ScopedSimdTarget(const ScopedSimdTarget&) = delete;
  ScopedSimdTarget& operator=(const ScopedSimdTarget&) = delete;
  ScopedSimdTarget(ScopedSimdTarget&&) = delete;
  ScopedSimdTarget& operator=(ScopedSimdTarget&&) = delete;
};

/// Reference posteriors from HMM::decode for each of the test pairs over sites [from, to).
std::vector<ReferencePosterior> decodeReferencePairs(HMM& hmm, const unsigned from, const unsigned to)
{
  const std::vector<unsigned long> hapsA = pairHapsA();
  const std::vector<unsigned long> hapsB = pairHapsB();
  std::vector<ReferencePosterior> references;
  references.reserve(hapsA.size());
  for (auto pair = 0ul; pair < hapsA.size(); ++pair) {
    references.push_back(decodeReference(hmm, hapsA.at(pair), hapsB.at(pair), from, to));
  }
  return references;
}

/**
 * Check that each pair's batched posterior over [from, to) matches the reference decoded over the same sites plus
 * burn-in. The batched path decodes the burn-in but stores only the window, so the reference must cover the same
 * extended region for the two to agree.
 */
void checkWindowMatchesReference(HMM& hmm, const Data& data, const unsigned from, const unsigned to)
{
  const unsigned decodedFrom = asmc::getFromPosition(data.geneticPositions, from, cmBurnIn);
  const unsigned decodedTo = asmc::getToPosition(data.geneticPositions, to, cmBurnIn);
  const std::vector<ReferencePosterior> references = decodeReferencePairs(hmm, decodedFrom, decodedTo);
  const DecodePairsReturnStruct batched = decodeBatched(hmm, pairHapsA(), pairHapsB(), from, to, cmBurnIn);

  for (auto pair = 0ul; pair < references.size(); ++pair) {
    INFO("window [" << from << ", " << to << "), pair " << pair);
    REQUIRE(batched.perPairPosteriors.at(pair).cols() == static_cast<Eigen::Index>(to - from));
    CHECK(maxAbsPosteriorDifference(batched.perPairPosteriors.at(pair), references.at(pair), from) <
          posteriorTolerance);
  }
}

/// Largest distance from 1 of any site's posterior sum.
float maxNormalisationError(const RowMajorPosterior& posterior)
{
  return (posterior.colwise().sum() - 1.f).abs().maxCoeff();
}

/// Largest relative difference between the stored posterior means of one pair and sum_k posterior(k) * time(k).
float maxPosteriorMeanError(const DecodePairsReturnStruct& results, const std::size_t pair,
                            const std::vector<float>& expectedTimes)
{
  const RowMajorPosterior& posterior = results.perPairPosteriors.at(pair);
  const Eigen::Map<const Eigen::ArrayXf> times(expectedTimes.data(), static_cast<Eigen::Index>(expectedTimes.size()));
  const Eigen::ArrayXf recomputedMeans = (posterior.colwise() * times).colwise().sum().transpose();
  const Eigen::ArrayXf storedMeans = results.perPairPosteriorMeans.row(static_cast<Eigen::Index>(pair)).transpose();
  return ((storedMeans - recomputedMeans).abs() / recomputedMeans).maxCoeff();
}

/// Number of sites at which the stored MAP state of one pair is not the state with the largest posterior.
long countMapMismatches(const DecodePairsReturnStruct& results, const std::size_t pair)
{
  const RowMajorPosterior& posterior = results.perPairPosteriors.at(pair);
  long mismatches = 0;
  for (Eigen::Index site = 0; site < posterior.cols(); ++site) {
    Eigen::Index mostLikelyState = 0;
    posterior.col(site).maxCoeff(&mostLikelyState);
    if (results.perPairMAPs(static_cast<Eigen::Index>(pair), site) != mostLikelyState) {
      ++mismatches;
    }
  }
  return mismatches;
}

/// How far HMM::decodeSummarize's output is from the MAP and posterior mean recomputed from the posterior.
struct SummaryErrors {
  long mapMismatches = 0;   // sites whose MAP time is not the expected time of the most likely state
  float maxMeanError = 0.f; // largest relative error in the posterior mean
};

/**
 * Compare HMM::decodeSummarize's per-site MAP times and posterior means with values recomputed from the posterior.
 *
 * @param summary (MAP times, posterior means) from HMM::decodeSummarize
 * @param posterior the posterior from HMM::decode for the same pair, indexed [state][site]
 * @param expectedTimes the expected coalescence time of each state
 */
SummaryErrors compareSummaryWithPosterior(const std::pair<std::vector<float>, std::vector<float>>& summary,
                                          const ReferencePosterior& posterior, const std::vector<float>& expectedTimes)
{
  const auto& [mapTimes, meanTimes] = summary;
  SummaryErrors errors;
  for (auto site = 0ul; site < mapTimes.size(); ++site) {
    auto mostLikelyState = 0ul;
    float mean = 0.f;
    for (auto state = 0ul; state < posterior.size(); ++state) {
      mean += posterior.at(state).at(site) * expectedTimes.at(state);
      if (posterior.at(state).at(site) > posterior.at(mostLikelyState).at(site)) {
        mostLikelyState = state;
      }
    }
    if (mapTimes.at(site) != expectedTimes.at(mostLikelyState)) {
      ++errors.mapMismatches;
    }
    errors.maxMeanError = std::max(errors.maxMeanError, std::abs(meanTimes.at(site) - mean) / mean);
  }
  return errors;
}

} // namespace

TEST_CASE("batched posteriors match HMM::decode under every SIMD target", "[HMM][SIMD]")
{
  const DecodingParams params = arrayParams();
  const Data data(params);
  HMM hmm(data, params);
  const auto numSites = static_cast<unsigned>(data.sites);
  const std::vector<ReferencePosterior> references = decodeReferencePairs(hmm, 0u, numSites);

  for (const int64_t target : hwy::SupportedAndGeneratedTargets()) {
    const ScopedSimdTarget onlyThisTarget(target);
    const DecodePairsReturnStruct batched = decodeBatched(hmm, pairHapsA(), pairHapsB(), 0u, numSites, cmBurnIn);

    for (auto pair = 0ul; pair < references.size(); ++pair) {
      INFO("SIMD target " << hwy::TargetName(target) << ", pair " << pair);
      CHECK(maxAbsPosteriorDifference(batched.perPairPosteriors.at(pair), references.at(pair), 0u) <
            posteriorTolerance);
      CHECK(maxNormalisationError(batched.perPairPosteriors.at(pair)) < posteriorTolerance);
    }
  }
}

TEST_CASE("a pair's posterior does not depend on the rest of its batch", "[HMM]")
{
  const DecodingParams params = arrayParams();
  const Data data(params);
  HMM hmm(data, params);
  const auto numSites = static_cast<unsigned>(data.sites);

  // 19 pairs fill one batch of 16 and leave 3 for a final batch padded up to a multiple of the SIMD lane count
  const auto numPairs = 19ul;
  std::vector<unsigned long> hapsA;
  std::vector<unsigned long> hapsB;
  hapsA.reserve(numPairs);
  hapsB.reserve(numPairs);
  for (auto pair = 0ul; pair < numPairs; ++pair) {
    hapsA.push_back(pair);
    hapsB.push_back(299ul - pair);
  }
  const DecodePairsReturnStruct together = decodeBatched(hmm, hapsA, hapsB, 0u, numSites, cmBurnIn);

  // Lanes are computed independently, so a pair should give bit-identical results in any batch
  for (const unsigned long pair : {0ul, 18ul}) {
    INFO("pair " << pair);
    const DecodePairsReturnStruct alone =
        decodeBatched(hmm, {hapsA.at(pair)}, {hapsB.at(pair)}, 0u, numSites, cmBurnIn);
    CHECK((alone.perPairPosteriors.at(0) == together.perPairPosteriors.at(pair)).all());
  }
}

TEST_CASE("windows at the ends of the chromosome match HMM::decode", "[HMM]")
{
  const DecodingParams params = arrayParams();
  const Data data(params);
  HMM hmm(data, params);
  const auto numSites = static_cast<unsigned>(data.sites);

  SECTION("first site only")
  {
    checkWindowMatchesReference(hmm, data, 0u, 1u);
  }
  SECTION("starting at the first site")
  {
    checkWindowMatchesReference(hmm, data, 0u, 200u);
  }
  SECTION("ending at the last site")
  {
    checkWindowMatchesReference(hmm, data, numSites - 200u, numSites);
  }
  SECTION("last site only")
  {
    checkWindowMatchesReference(hmm, data, numSites - 1u, numSites);
  }
}

// Known bug: in sequence mode, HMM::forwardBatch and HMM::backwardBatch take two steps per site and assign the
// intermediate result with `previousAlpha = nextAlpha` and `lastComputedBeta = previousBeta`. These are Eigen::Maps
// into the alpha and beta buffers, so the assignment copies the intermediate values over the stored alpha of the
// previous site and the stored beta of the next site, which then enter that site's posterior. Batched posteriors
// differ from HMM::decode by up to about 0.5 as a result. Remove [!shouldfail] once that is fixed.
TEST_CASE("sequence-mode posteriors match HMM::decode", "[HMM][!shouldfail]")
{
  // The sequence data has about 97,000 sites, so decode a window to keep the test fast
  const DecodingParams params = sequenceParams();
  const Data data(params);
  HMM hmm(data, params);
  checkWindowMatchesReference(hmm, data, 5000u, 7000u);
}

TEST_CASE("batched posterior means and MAPs agree with the posteriors", "[HMM]")
{
  const DecodingParams params = arrayParams();
  const Data data(params);
  HMM hmm(data, params);
  const auto numSites = static_cast<unsigned>(data.sites);
  const std::vector<float>& expectedTimes = hmm.getDecodingQuantities().expectedTimes;

  const DecodePairsReturnStruct batched = decodeBatched(hmm, pairHapsA(), pairHapsB(), 0u, numSites, cmBurnIn);
  for (auto pair = 0ul; pair < batched.perPairPosteriors.size(); ++pair) {
    INFO("pair " << pair);
    CHECK(maxPosteriorMeanError(batched, pair, expectedTimes) < 1e-5f);
    CHECK(countMapMismatches(batched, pair) == 0);
  }
}

TEST_CASE("HMM::decodeSummarize agrees with HMM::decode", "[HMM]")
{
  const DecodingParams params = arrayParams();
  const Data data(params);
  HMM hmm(data, params);
  const std::vector<float>& expectedTimes = hmm.getDecodingQuantities().expectedTimes;

  const PairObservations observations = hmm.makePairObs(1, 0u, 2, 0u);
  const SummaryErrors errors =
      compareSummaryWithPosterior(hmm.decodeSummarize(observations), hmm.decode(observations), expectedTimes);
  CHECK(errors.mapMismatches == 0);
  CHECK(errors.maxMeanError < 1e-5f);
}
