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

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "ASMC.hpp"
#include "FileUtils.hpp"

#include <Eigen/Core>
#include <fmt/format.h>

#include <sstream>

TEST_CASE("test ASMC regression", "[HMM_regression]")
{
  // we only needed to set doPosteriorSums to true, but because C++ does
  // not have keyword arguments we need to go through everything
  DecodingParams params(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                        ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz",
                        "",      // _outFileRoot
                        1,       // _jobs
                        1,       // _jobInd
                        "array", // _decodingModeString
                        false,   // _decodingSequence
                        true,    // _usingCSFS
                        false,   // _compress
                        false,   // _useAncestral
                        0.f,     // _skipCSFSdistance
                        false,   // _noBatches
                        true     // _doPosteriorSums
  );

  std::vector<unsigned long> indToDecodeA = {1ul, 2ul, 3ul};
  std::vector<unsigned long> indToDecodeB = {2ul, 3ul, 4ul};

  ASMC::ASMC asmc(params);
  asmc.setStorePerPairPosteriorMean();
  asmc.setStorePerPairMap();
  asmc.decodePairs(indToDecodeA, indToDecodeB);

  auto res = asmc.getRefOfResults();

  SECTION("regression test per pair posterior means")
  {

    CHECK(res.perPairPosteriorMeans.rows() == 3ll);
    CHECK(res.perPairPosteriorMeans.cols() == 6760ll);

    std::string regressionFile = ASMC_DATA_DIR "/testing/asmc/regression/regression.perPairPosteriorMeans.gz";
    FileUtils::AutoGzIfstream fin;
    fin.openOrExit(regressionFile);

    for (auto rowIdx = 0ul; rowIdx < indToDecodeA.size(); ++rowIdx) {
      std::string line;
      getline(fin, line);
      std::istringstream iss(line);
      std::vector<float> rowAsFloats = {std::istream_iterator<float>(iss), std::istream_iterator<float>()};

      CHECK(rowAsFloats.size() == 6760ul);
      for (auto colIdx = 0ul; colIdx < rowAsFloats.size(); ++colIdx) {
        CHECK(res.perPairPosteriorMeans(rowIdx, colIdx) == Approx(rowAsFloats.at(colIdx)));
      }
    }
  }

  SECTION("regression test per pair MAP")
  {
    CHECK(res.perPairMAPs.rows() == 3ll);
    CHECK(res.perPairMAPs.cols() == 6760ll);

    std::string regressionFile = ASMC_DATA_DIR "/testing/asmc/regression/regression.perPairMAP.gz";
    FileUtils::AutoGzIfstream fin;
    fin.openOrExit(regressionFile);

    for (auto rowIdx = 0ul; rowIdx < indToDecodeA.size(); ++rowIdx) {
      std::string line;
      getline(fin, line);
      std::istringstream iss(line);
      std::vector<int> rowAsInts = {std::istream_iterator<int>(iss), std::istream_iterator<int>()};

      CHECK(rowAsInts.size() == 6760ul);
      for (auto colIdx = 0ul; colIdx < rowAsInts.size(); ++colIdx) {
        CHECK(res.perPairMAPs(rowIdx, colIdx) == rowAsInts.at(colIdx));
      }
    }
  }
}
