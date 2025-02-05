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

#include "Simd.hpp"

TEST_CASE("Print runtime SIMD information", "[SIMD]")
{
  asmc::printRuntimeSimdInfo();
}

TEST_CASE("Check batchsize validation", "[SIMD]")
{
  const int numSimdLanes = asmc::getNumSimdLanes();
  REQUIRE(numSimdLanes > 1);

  const int compositeBatchSize = 2 * numSimdLanes;
  const int nonCompositeBatchSize = compositeBatchSize + 1;

  REQUIRE_NOTHROW(asmc::validateBatchSize(compositeBatchSize));
  REQUIRE_THROWS_AS(asmc::validateBatchSize(nonCompositeBatchSize), std::runtime_error);
}
