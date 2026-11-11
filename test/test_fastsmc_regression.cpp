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

// FastSMC regression tests: each runs FastSMC on the example data under one configuration and compares the IBD output
// with a snapshot. See regression_snapshots.hpp for how to regenerate the snapshots.

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include "BinaryDataReader.hpp"
#include "DecodingParams.hpp"
#include "FastSMC.hpp"
#include "regression_snapshots.hpp"
#include "test_fixtures.hpp"

namespace
{

using asmc::test::checkLinesMatchSnapshot;
using asmc::test::floatOutput;

/**
 * FastSMC parameters shared by the regression tests: the example data, hashing on, per-segment age estimates and
 * segment lengths on. Output goes to `outFileRoot`. Callers adjust these and then call validateParamsFastSMC.
 */
DecodingParams fastSmcParams(const std::string& outFileRoot)
{
  DecodingParams params;
  params.decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";
  params.inFileRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
  params.outFileRoot = outFileRoot;
  params.decodingModeString = "array";
  params.usingCSFS = true;
  params.batchSize = 32;
  params.recallThreshold = 3;
  params.min_m = 1.5;
  params.hashing = true;
  params.FastSMC = true;
  params.BIN_OUT = false;
  params.outputIbdSegmentLength = true;
  params.time = 50;
  params.noConditionalAgeEstimates = true;
  params.doPerPairMAP = true;
  params.doPerPairPosteriorMean = true;
  params.useKnownSeed = true;
  params.hashingOnly = false;
  return params;
}

/// Path of the output FastSMC writes for `params`, in text or binary form.
std::string outputPath(const DecodingParams& params)
{
  return params.outFileRoot + "." + std::to_string(params.jobInd) + "." + std::to_string(params.jobs) +
         (params.BIN_OUT ? ".FastSMC.bibd.gz" : ".FastSMC.ibd.gz");
}

/// Validate `params`, run FastSMC, and return its text output, one line per IBD segment.
std::vector<std::string> runFastSmc(DecodingParams params)
{
  params.validateParamsFastSMC();
  ASMC::FastSMC fastSMC(params);
  fastSMC.run();
  return asmc::test::readLines(outputPath(params));
}

/// FastSMC without hashing, for job 8 of 400: 113 pairs, so the last batch of 32 is padded.
DecodingParams noHashingJobParams(const std::string& outFileRoot)
{
  DecodingParams params = fastSmcParams(outFileRoot);
  params.hashing = false;
  params.jobs = 400;
  params.jobInd = 8;
  return params;
}

} // namespace

TEST_CASE("FastSMC regression: hashing", "[FastSMC_regression]")
{
  const asmc::test::TempDirectory temp;
  checkLinesMatchSnapshot(runFastSmc(fastSmcParams(temp.file("out").string())), "fastsmc_hashing.ibd.gz", floatOutput);
}

TEST_CASE("FastSMC regression: hashing only", "[FastSMC_regression]")
{
  const asmc::test::TempDirectory temp;
  DecodingParams params = fastSmcParams(temp.file("out").string());
  params.hashingOnly = true;
  params.outputIbdSegmentLength = false;
  params.doPerPairMAP = false;
  params.doPerPairPosteriorMean = false;
  checkLinesMatchSnapshot(runFastSmc(params), "fastsmc_hashing_only.ibd.gz", floatOutput);
}

TEST_CASE("FastSMC regression: no hashing, one job", "[FastSMC_regression]")
{
  const asmc::test::TempDirectory temp;
  checkLinesMatchSnapshot(runFastSmc(noHashingJobParams(temp.file("out").string())), "fastsmc_no_hashing_job.ibd.gz",
                          floatOutput);
}

TEST_CASE("FastSMC regression: hashing, one job of several, conditional age estimates", "[FastSMC_regression]")
{
  // Job 2 of 4 compares the second window of samples with the first
  const asmc::test::TempDirectory temp;
  DecodingParams params = fastSmcParams(temp.file("out").string());
  params.jobs = 4;
  params.jobInd = 2;
  params.noConditionalAgeEstimates = false;
  checkLinesMatchSnapshot(runFastSmc(params), "fastsmc_hashing_job_conditional.ibd.gz", floatOutput);
}

TEST_CASE("FastSMC regression: binary output matches text output", "[FastSMC_regression]")
{
  // Same configuration as "no hashing, one job", but writing binary output, read back with BinaryDataReader
  const asmc::test::TempDirectory temp;
  DecodingParams params = noHashingJobParams(temp.file("out").string());
  params.BIN_OUT = true;
  params.validateParamsFastSMC();
  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  std::vector<std::string> lines;
  BinaryDataReader reader(outputPath(params));
  while (reader.moreLinesInFile()) {
    lines.push_back(reader.getNextLine().toString());
  }

  // Never regenerate from binary output: the snapshot belongs to the text test
  if (!asmc::test::regeneratingSnapshots()) {
    checkLinesMatchSnapshot(lines, "fastsmc_no_hashing_job.ibd.gz", floatOutput);
  }
}
