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

#include "DecodingParams.hpp"

TEST_CASE("test DecodingParams", "[DecodingParams]")
{
  std::string inFileRoot = ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array";
  std::string decodingQuantFile = ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz";

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
        nan("") // _skipCSFSdistance, override default
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
