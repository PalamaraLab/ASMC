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

#include <utility>
#include <vector>

#include "HMM.hpp"

TEST_CASE("test hmm functions", "[HMM]")
{
  DecodingParams params(
      ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
      ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  Data data(params);
  HMM hmm(data, params);

  REQUIRE(data.individuals.size() > 20);

  SECTION("test decode pair summarize")
  {
    PairObservations pairObs = hmm.makePairObs(1, 0, 2, 0);
    std::vector<std::vector<float>> decodeResult = hmm.decode(pairObs);
    std::pair<std::vector<float>, std::vector<float>> decodeSummary = hmm.decodeSummarize(pairObs);
    // check that the MAP and posterior mean are the same length
    REQUIRE(decodeSummary.first.size() == decodeSummary.second.size());
    REQUIRE(decodeSummary.first.size() == decodeResult[0].size());
  }

  SECTION("test decode pair")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePair(0, 9);
    REQUIRE(hmm.getBatchBuffer().size() == 4);
    hmm.decodePair(1, 1);
    REQUIRE(hmm.getBatchBuffer().size() == 5);
  }

  SECTION("test decode pairs")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePairs({ 0, 1 }, { 9, 1 });
    REQUIRE(hmm.getBatchBuffer().size() == 5);
  }

  SECTION("test finishDecoding")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePair(0, 9);
    REQUIRE(hmm.getBatchBuffer().size() == 4);
    hmm.finishDecoding();
    REQUIRE(hmm.getBatchBuffer().size() == 0);
  }

  SECTION("test fill up buffer")
  {
    // default batch size is 64
    for (int i = 1; i <= 64 / 4; ++i) {
      hmm.decodePair(0, i);
    }

    // buffer should be empty now
    REQUIRE(hmm.getBatchBuffer().size() == 0);
  }
}
