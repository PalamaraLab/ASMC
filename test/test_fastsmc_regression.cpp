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

#include <iostream>
#include <string>
#include <vector>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "FileUtils.hpp"
#include "HMM.hpp"
#include "Timer.hpp"



TEST_CASE("test FastSMC + hashing regression test", "[FastSMC_regression]")
{
  DecodingParams params;
  params.decodingQuantFile = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/example.decodingQuantities.gz";
  params.inFileRoot = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/example";
  params.outFileRoot = "/tmp/FastSMCresults";
  params.decodingModeString = "array";
  params.foldData = true;
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

  assert(params.validateParamsFastSMC());

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  SECTION("regression test")
  {
    const auto expectedNumLines = 1524ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_FILE_DIR "/FASTSMC_EXAMPLE/regression_output.ibd.gz");
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
      REQUIRE(regressionLines.at(lineNum) == generatedLines.at(lineNum));
    }
  }
}


TEST_CASE("test FastSMC without hashing regression test", "[FastSMC_regression]")
{
  DecodingParams params;
  params.decodingQuantFile = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/example.decodingQuantities.gz";
  params.inFileRoot = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/example";
  params.outFileRoot = "/tmp/FastSMCresults";
  params.decodingModeString = "array";
  params.foldData = true;
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
  params.jobs = 9;
  params.useKnownSeed = true;

  assert(params.validateParamsFastSMC());

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  SECTION("regression test")
  {
    const auto expectedNumLines = 2986ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_FILE_DIR "/FASTSMC_EXAMPLE/regression_output_no_hashing.ibd.gz");
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
      fin_generated.openOrExit(params.outFileRoot + ".7.9.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == expectedNumLines);
    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      REQUIRE(regressionLines.at(lineNum) == generatedLines.at(lineNum));
    }
  }
}