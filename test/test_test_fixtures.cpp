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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <filesystem>

#include "test_fixtures.hpp"

using Catch::Approx;

TEST_CASE("test fixture helpers", "[fixtures]")
{
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);

  CHECK(std::filesystem::is_regular_file(input.sampleFile));
  CHECK(std::filesystem::is_regular_file(input.hapFile));
  CHECK(std::filesystem::is_regular_file(input.mapFile));

  const auto posterior = asmc::test::makeColumnNormalisedPosterior(3u, 4u);
  REQUIRE(posterior.rows() == 3);
  REQUIRE(posterior.cols() == 4);
  CHECK(posterior.col(0).sum() == Approx(1.f));
  CHECK(posterior.col(3).sum() == Approx(1.f));

  const auto buffer = asmc::test::makeBatchedPosteriorBuffer({posterior}, 2);
  REQUIRE(buffer.size() == 3 * 4 * 2);
  CHECK(buffer(0) == posterior(0, 0));
  CHECK(buffer(1) == 0.f);

  const auto obs = asmc::test::makePairObservations({true, false}, {false, true});
  CHECK(obs.obsBits.size() == 2ul);
  CHECK(obs.homMinorBits.size() == 2ul);
}
