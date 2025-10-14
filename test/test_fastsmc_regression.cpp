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

#include <string>
#include <vector>

#include "utils/StringUtils.hpp"

#include "DecodingParams.hpp"
#include "FastSMC.hpp"
#include "FileUtils.hpp"
#include "Simd.hpp"
#include "Timer.hpp"

TEST_CASE("test FastSMC with hashing regression test", "[FastSMC_regression]")
{
  // This test is designed to work with AVX2. It may work with other SIMD instruction sets, but it is possible that
  // the test could fail with minor numerical differences.
  asmc::warnIfSimdMismatch("AVX2");

  DecodingParams params;
  params.decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";
  params.inFileRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
  params.outFileRoot = "/tmp/FastSMC";
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

  params.validateParamsFastSMC();

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  SECTION("regression test")
  {
    const auto expectedNumLines = 1587ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_DATA_DIR "/testing/fastsmc/regression/regression_output.ibd.gz");
      for (std::string f_line; getline(fin_regression, f_line);) {
        regressionLines.emplace_back(f_line);
      }
      fin_regression.close();
    }

    // Read lines from generated test output into a vector of strings
    std::vector<std::string> generatedLines;
    generatedLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_generated;
      fin_generated.openOrExit(params.outFileRoot + ".1.1.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == expectedNumLines);
    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      auto reg_split = asmc::splitTextByDelimiter(regressionLines.at(lineNum), "\t");
      auto gen_split = asmc::splitTextByDelimiter(generatedLines.at(lineNum), "\t");

      CHECK(reg_split.size() == 13ul);
      CHECK(gen_split.size() == 13ul);

      // Check the string fields
      CHECK(reg_split.at(0) == gen_split.at(0));
      CHECK(reg_split.at(1) == gen_split.at(1));
      CHECK(reg_split.at(3) == gen_split.at(3));
      CHECK(reg_split.at(4) == gen_split.at(4));

      // Check the integer fields
      CHECK(std::stoul(reg_split.at(2)) == std::stoul(gen_split.at(2)));
      CHECK(std::stoul(reg_split.at(5)) == std::stoul(gen_split.at(5)));
      CHECK(std::stoul(reg_split.at(6)) == std::stoul(gen_split.at(6)));
      CHECK(std::stoul(reg_split.at(7)) == std::stoul(gen_split.at(7)));
      CHECK(std::stoul(reg_split.at(8)) == std::stoul(gen_split.at(8)));

      // Check the float fields
      CHECK(std::stold(reg_split.at(9)) == Approx(std::stold(gen_split.at(9))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(10)) == Approx(std::stold(gen_split.at(10))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(11)) == Approx(std::stold(gen_split.at(11))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(12)) == Approx(std::stold(gen_split.at(12))).epsilon(0.001));
    }
  }
}

TEST_CASE("test FastSMC with hashing only (GERMLINE2)", "[FastSMC_regression]")
{
  // This test is designed to work with AVX2. It may work with other SIMD instruction sets, but it is possible that
  // the test could fail with minor numerical differences.
  asmc::warnIfSimdMismatch("AVX2");

  DecodingParams params;
  params.decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";
  params.inFileRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
  params.outFileRoot = "/tmp/FastSMCHashing";
  params.decodingModeString = "array";
  params.usingCSFS = true;
  params.batchSize = 32;
  params.recallThreshold = 3;
  params.min_m = 1.5;
  params.hashing = true;
  params.FastSMC = true;
  params.BIN_OUT = false;
  params.outputIbdSegmentLength = false;
  params.time = 50;
  params.noConditionalAgeEstimates = true;
  params.doPerPairMAP = false;
  params.doPerPairPosteriorMean = false;
  params.useKnownSeed = true;
  params.hashingOnly = true;

  params.validateParamsFastSMC();

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  SECTION("regression test")
  {
    const auto expectedNumLines = 495ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_DATA_DIR "/testing/fastsmc/regression/regression_output_hashing.ibd.gz");
      for (std::string f_line; getline(fin_regression, f_line);) {
        regressionLines.emplace_back(f_line);
      }
      fin_regression.close();
    }

    // Read lines from generated test output into a vector of strings
    std::vector<std::string> generatedLines;
    generatedLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_generated;
      fin_generated.openOrExit(params.outFileRoot + ".1.1.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == expectedNumLines);
    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      auto reg_split = asmc::splitTextByDelimiter(regressionLines.at(lineNum), "\t");
      auto gen_split = asmc::splitTextByDelimiter(generatedLines.at(lineNum), "\t");

      CHECK(reg_split.size() == 9ul);
      CHECK(gen_split.size() == 9ul);

      // Check the string fields
      CHECK(reg_split.at(0) == gen_split.at(0));
      CHECK(reg_split.at(1) == gen_split.at(1));
      CHECK(reg_split.at(3) == gen_split.at(3));
      CHECK(reg_split.at(4) == gen_split.at(4));

      // Check the integer fields
      CHECK(std::stoul(reg_split.at(2)) == std::stoul(gen_split.at(2)));
      CHECK(std::stoul(reg_split.at(5)) == std::stoul(gen_split.at(5)));
      CHECK(std::stoul(reg_split.at(6)) == std::stoul(gen_split.at(6)));
      CHECK(std::stoul(reg_split.at(7)) == std::stoul(gen_split.at(7)));
      CHECK(std::stoul(reg_split.at(8)) == std::stoul(gen_split.at(8)));
    }
  }
}


TEST_CASE("test FastSMC without hashing regression test", "[FastSMC_regression]")
{
  // This test is designed to work with AVX2. It may work with other SIMD instruction sets, but it is possible that
  // the test could fail with minor numerical differences.
  asmc::warnIfSimdMismatch("AVX2");

  DecodingParams params;
  params.decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";
  params.inFileRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
  params.outFileRoot = "/tmp/FastSMCresults";
  params.decodingModeString = "array";
  params.usingCSFS = true;
  params.batchSize = 32;
  params.recallThreshold = 3;
  params.min_m = 1.5;
  params.hashing = false;
  params.FastSMC = true;
  params.BIN_OUT = false;
  params.outputIbdSegmentLength = true;
  params.time = 50;
  params.noConditionalAgeEstimates = true;
  params.doPerPairMAP = true;
  params.doPerPairPosteriorMean = true;
  params.jobInd = 7;
  params.jobs = 25;
  params.useKnownSeed = true;

  params.validateParamsFastSMC();

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  SECTION("regression test")
  {
    const auto expectedNumLines = 498ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_DATA_DIR "/testing/fastsmc/regression/regression_output_no_hashing.ibd.gz");
      for (std::string f_line; getline(fin_regression, f_line);) {
        regressionLines.emplace_back(f_line);
      }
      fin_regression.close();
    }

    // Read lines from generated test output into a vector of strings
    std::vector<std::string> generatedLines;
    generatedLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_generated;
      fin_generated.openOrExit(params.outFileRoot + ".7.25.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == expectedNumLines);
    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      auto reg_split = asmc::splitTextByDelimiter(regressionLines.at(lineNum), "\t");
      auto gen_split = asmc::splitTextByDelimiter(generatedLines.at(lineNum), "\t");

      CHECK(reg_split.size() == 13ul);
      CHECK(gen_split.size() == 13ul);

      // Check the string fields
      CHECK(reg_split.at(0) == gen_split.at(0));
      CHECK(reg_split.at(1) == gen_split.at(1));
      CHECK(reg_split.at(3) == gen_split.at(3));
      CHECK(reg_split.at(4) == gen_split.at(4));

      // Check the integer fields
      CHECK(std::stoul(reg_split.at(2)) == std::stoul(gen_split.at(2)));
      CHECK(std::stoul(reg_split.at(5)) == std::stoul(gen_split.at(5)));
      CHECK(std::stoul(reg_split.at(6)) == std::stoul(gen_split.at(6)));
      CHECK(std::stoul(reg_split.at(7)) == std::stoul(gen_split.at(7)));
      CHECK(std::stoul(reg_split.at(8)) == std::stoul(gen_split.at(8)));

      // Check the float fields
      CHECK(std::stold(reg_split.at(9)) == Approx(std::stold(gen_split.at(9))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(10)) == Approx(std::stold(gen_split.at(10))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(11)) == Approx(std::stold(gen_split.at(11))).epsilon(0.001));
      CHECK(std::stold(reg_split.at(12)) == Approx(std::stold(gen_split.at(12))).epsilon(0.001));
    }
  }
}
