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

#include <string>

#include "DecodingQuantities.hpp"

using Catch::Matchers::ContainsSubstring;

TEST_CASE("test validate decoding quantities file", "[DecodingQuantities]")
{
  std::string nonExistentDecodingQuantitiesFile = ASMC_DATA_DIR "/random_nonexistent_file.txt";
  std::string goodDecodingQuantitiesFile = ASMC_DATA_DIR "/testing/asmc/decoding_quantities_good.txt";
  std::string badDecodingQuantitiesFile = ASMC_DATA_DIR "/testing/asmc/decoding_quantities_bad.txt";

  SECTION("test nonexistent file")
  {
    CHECK_THROWS_WITH(DecodingQuantities{nonExistentDecodingQuantitiesFile},
                      ContainsSubstring("random_nonexistent_file.txt does not exist"));
  }

  SECTION("test good file")
  {
    CHECK_NOTHROW(DecodingQuantities{goodDecodingQuantitiesFile});
  }

  SECTION("test bad file")
  {
    CHECK_THROWS_WITH(DecodingQuantities{badDecodingQuantitiesFile},
                      ContainsSubstring("decoding_quantities_bad.txt does not seem to contain the correct information") &&
                          ContainsSubstring("but instead found \"this file does not start with"));
  }
}
