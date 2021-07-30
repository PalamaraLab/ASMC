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

#include "HMM.hpp"
#include "FileUtils.hpp"
#include <sstream>

TEST_CASE("test hmm with regression test", "[HMM_regression]")
{
  // we only needed to set doPosteriorSums to true, but because C++ does
  // not have keyword arguments we need to go through everything
  DecodingParams params(
      ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array",
      ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz",
      "", // _outFileRoot
      1, // _jobs
      1, // _jobInd
      "array", // _decodingModeString
      false, // _decodingSequence
      true, // _usingCSFS
      false, // _compress
      false, // _useAncestral
      0.f, // _skipCSFSdistance
      false, // _noBatches
      true // _doPosteriorSums
      );

  Data data(params);
  HMM hmm(data, params);

  REQUIRE(data.individuals.size() > 20);

  SECTION("regression test")
  {
    std::string regressionFile = ASMC_FILE_DIR "/../ASMC_SRC/TESTS/data/regression_test_original.gz";
    FileUtils::AutoGzIfstream fin;
    fin.openOrExit(regressionFile);
    hmm.decodeAll(params.jobs, params.jobInd);
    const DecodingReturnValues& decodingReturnValues = hmm.getDecodingReturnValues();
    int pos = 0;
    for( std::string line; getline(fin, line); )
    {
      std::ostringstream oss;
      for (uint k = 0; k < hmm.getDecodingQuantities().states; k++) {
        if (k) oss << "\t";
        oss << decodingReturnValues.sumOverPairs(pos,k);
      }
      REQUIRE(oss.str() == line);
      pos++;
    }
    REQUIRE(pos == decodingReturnValues.sumOverPairs.rows());
  }
}
