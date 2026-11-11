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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "DecodingParams.hpp"

using Catch::Matchers::ContainsSubstring;

namespace
{

/// Run DecodingParams::processCommandLineArgsFastSMC on `options`, preceded by a program name and the two mandatory
/// options, as main() would receive them from the command line.
bool parseFastSmcCommandLine(DecodingParams& params, const std::vector<std::string>& options)
{
  std::vector<std::string> args = {"FastSMC", "--inFileRoot", "in", "--outFileRoot", "out"};
  args.insert(args.end(), options.begin(), options.end());

  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& arg : args) {
    argv.push_back(arg.data());
  }
  return params.processCommandLineArgsFastSMC(static_cast<int>(argv.size()), argv.data());
}

} // namespace

TEST_CASE("test DecodingParams", "[DecodingParams]")
{
  std::string inFileRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
  std::string decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";

  SECTION("test array folded") {
    DecodingParams params(inFileRoot, decodingQuantFile);
    REQUIRE(params.decodingMode == DecodingMode::arrayFolded);
    REQUIRE(params.compress == false);
  }

  SECTION("test sequence folded") {
    DecodingParams params(inFileRoot, decodingQuantFile,
        "", // _outFileRoot
        1, // _jobs
        1, // _jobInd
        "sequence", // _decodingModeString, override default
        false, // _decodingSequence
        true, // _usingCSFS
        true, // _compress, override default
        false, // _useAncestral
        std::numeric_limits<float>::quiet_NaN() // _skipCSFSdistance, override default
    );
    REQUIRE(params.decodingMode == DecodingMode::sequenceFolded);
    REQUIRE(params.compress == true);
  }

  SECTION("test sequence") {
    DecodingParams params(inFileRoot, decodingQuantFile,
        "", // _outFileRoot
        1, // _jobs
        1, // _jobInd
        "sequence", // _decodingModeString, override default
        false, // _decodingSequence
        true, // _usingCSFS
        false, // _compress
        true // _useAncestral, override default
    );
    REQUIRE(params.decodingMode == DecodingMode::sequence);
  }
}

TEST_CASE("FastSMC command line honours --mode", "[DecodingParams]")
{
  DecodingParams params;

  SECTION("array by default")
  {
    REQUIRE(parseFastSmcCommandLine(params, {}));
    CHECK_FALSE(params.decodingSequence);
    CHECK(params.decodingMode == DecodingMode::arrayFolded);
  }

  SECTION("sequence when requested")
  {
    REQUIRE(parseFastSmcCommandLine(params, {"--mode", "sequence"}));
    CHECK(params.decodingSequence);
    CHECK(params.decodingMode == DecodingMode::sequenceFolded);
  }
}

TEST_CASE("FastSMC options that default to on can be turned off", "[DecodingParams]")
{
  DecodingParams params;

  SECTION("on by default")
  {
    REQUIRE(parseFastSmcCommandLine(params, {}));
    CHECK(params.hashing);
    CHECK(params.outputIbdSegmentLength);
    CHECK(params.doPerPairMAP);
    CHECK(params.doPerPairPosteriorMean);
  }

  SECTION("a bare flag still turns an option on")
  {
    REQUIRE(parseFastSmcCommandLine(params, {"--hashing"}));
    CHECK(params.hashing);
  }

  SECTION("=0 turns an option off")
  {
    REQUIRE(parseFastSmcCommandLine(
        params, {"--hashing=0", "--segmentLength=0", "--perPairMAP=0", "--perPairPosteriorMeans=0"}));
    CHECK_FALSE(params.hashing);
    CHECK_FALSE(params.outputIbdSegmentLength);
    CHECK_FALSE(params.doPerPairMAP);
    CHECK_FALSE(params.doPerPairPosteriorMean);
  }
}

TEST_CASE("compress works without setting skipCSFSdistance", "[DecodingParams]")
{
  const std::string inFileRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
  const std::string decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";
  const float nan = std::numeric_limits<float>::quiet_NaN();

  SECTION("ASMC constructor with the default distance")
  {
    const DecodingParams params(inFileRoot, decodingQuantFile, "", 1, 1, "array", false, true, /*_compress=*/true);
    CHECK(std::isinf(params.skipCSFSdistance));
  }

  SECTION("ASMC constructor with an explicit distance is rejected")
  {
    CHECK_THROWS_AS(DecodingParams(inFileRoot, decodingQuantFile, "", 1, 1, "array", false, true, /*_compress=*/true,
                                   false, /*_skipCSFSdistance=*/0.1f),
                    std::invalid_argument);
  }

  SECTION("FastSMC parameters can be validated more than once")
  {
    DecodingParams params(inFileRoot, decodingQuantFile, "out", /*_fastSMC=*/true);
    params.compress = true;
    params.skipCSFSdistance = nan;
    CHECK_NOTHROW(params.validateParamsFastSMC());
    CHECK_NOTHROW(params.validateParamsFastSMC());
    CHECK(std::isinf(params.skipCSFSdistance));
  }
}

TEST_CASE("invalid parameters throw instead of exiting", "[DecodingParams]")
{
  const std::string inFileRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
  const std::string decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";

  CHECK_THROWS_WITH(DecodingParams(inFileRoot, decodingQuantFile, "", 1, 1, "nonsense"),
                    ContainsSubstring("Unknown decoding mode \"nonsense\""));
  CHECK_THROWS_AS(DecodingParams(inFileRoot, decodingQuantFile, "out", /*_fastSMC=*/false), std::invalid_argument);

  DecodingParams params;
  CHECK_FALSE(parseFastSmcCommandLine(params, {"--mode", "nonsense"}));
}
