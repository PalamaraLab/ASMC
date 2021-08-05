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

#include "ASMC.hpp"

#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/ostream.h>

TEST_CASE("test ASMC decodeAllInJob", "[ASMC]")
{
  DecodingParams params(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                        ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  ASMC::ASMC asmc(params);

  //  auto result = asmc.decodeAllInJob();
  //
  //  SECTION("test decode pair summarize")
  //  {
  //    REQUIRE(result.sumOverPairs.size() == 466440ul);
  //  }
}

TEST_CASE("test ASMC decodePairs", "[ASMC]")
{
  ASMC::ASMC asmc(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                  ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");

  asmc.setStorePerPairMap();
  asmc.setStorePerPairPosterior();
  asmc.setStorePerPairPosteriorMean();
  asmc.setStorePerPairMap();

  std::vector<unsigned long> indA = {1, 2, 3};
  std::vector<unsigned long> indB = {2, 3, 4};
  asmc.decodePairs(indA, indB);
  auto result = asmc.getRefOfResults();

  SECTION("test decode pair summarize")
  {
    REQUIRE(result.perPairIndices.size() == 3ul);

    REQUIRE(result.perPairPosteriorMeans(0, 0) == Approx(15968.91016f).margin(1.f));
    REQUIRE(result.perPairPosteriorMeans(1, 8) == Approx(27963.49805f).margin(1.f));
    REQUIRE(result.perPairPosteriorMeans(2, 29) == Approx(48573.32812f).margin(1.f));

    REQUIRE(result.perPairMAPs(0, 0) == 29);
    REQUIRE(result.perPairMAPs(1, 1234) == 65);
    REQUIRE(result.perPairMAPs(2, 7) == 33);

    for (Eigen::Index idx = 0ll; idx < result.perPairPosteriors.size(); ++idx) {
      REQUIRE((result.perPairPosteriors.at(idx).colwise().sum() - result.perPairPosteriorMeans.row(idx)).isZero(1e-2));
    }
  }
}
