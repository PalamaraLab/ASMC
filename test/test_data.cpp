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

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "Data.hpp"
#include "DecodingParams.hpp"

namespace
{

constexpr const char* asmcExampleRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
constexpr const char* fastSmcExampleRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
constexpr const char* asmcDecodingQuantities =
    ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";
constexpr const char* fastSmcDecodingQuantities =
    ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";

/// FastSMC parameters for job `jobInd` of `jobs`, with or without the hashing step.
DecodingParams fastSmcParams(const bool hashing, const int jobs, const int jobInd)
{
  DecodingParams params(fastSmcExampleRoot, fastSmcDecodingQuantities, "out", /*_fastSMC=*/true);
  params.hashing = hashing;
  params.jobs = jobs;
  params.jobInd = jobInd;
  params.validateParamsFastSMC();
  return params;
}

} // namespace

TEST_CASE("ASMC loads every sample, whatever the job", "[Data]")
{
  const Data allJobs(DecodingParams(asmcExampleRoot, asmcDecodingQuantities));
  REQUIRE(allJobs.individuals.size() == allJobs.sampleSize);

  // HMM::decodeAll divides the pairs between jobs, so every job needs every sample, with the right genotypes
  const int jobs = 4;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(DecodingParams(asmcExampleRoot, asmcDecodingQuantities, "", jobs, jobInd));
    REQUIRE(oneJob.individuals.size() == allJobs.individuals.size());
    for (auto ind = 0ul; ind < allJobs.individuals.size(); ++ind) {
      CHECK(oneJob.individuals.at(ind).genotype1 == allJobs.individuals.at(ind).genotype1);
      CHECK(oneJob.individuals.at(ind).genotype2 == allJobs.individuals.at(ind).genotype2);
    }
  }
}

TEST_CASE("FastSMC without hashing loads every sample", "[Data]")
{
  const int jobs = 4;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(fastSmcParams(/*hashing=*/false, jobs, jobInd));
    CHECK(oneJob.individuals.size() == oneJob.sampleSize);
  }
}

TEST_CASE("FastSMC with hashing loads only its job's windows, and the jobs cover every sample", "[Data]")
{
  // With 4 jobs there are two windows. Jobs 1 and 4 compare a window with itself and load only that window; jobs 2
  // and 3 compare the two windows with each other and load both.
  const int jobs = 4;
  std::vector<bool> loadedByAnyJob;
  bool someJobLoadsASubset = false;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(fastSmcParams(/*hashing=*/true, jobs, jobInd));
    loadedByAnyJob.resize(oneJob.sampleSize, false);

    auto numInJob = 0ul;
    for (auto sample = 0ul; sample < oneJob.sampleSize; ++sample) {
      if (oneJob.isSampleInJob(sample)) {
        ++numInJob;
        loadedByAnyJob.at(sample) = true;
      }
    }
    CHECK(oneJob.individuals.size() == numInJob);
    someJobLoadsASubset = someJobLoadsASubset || numInJob < oneJob.sampleSize;
  }

  CHECK(someJobLoadsASubset);
  CHECK(std::find(loadedByAnyJob.begin(), loadedByAnyJob.end(), false) == loadedByAnyJob.end());
}
