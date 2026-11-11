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
#include "test_fixtures.hpp"

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

TEST_CASE("malformed decoding quantities throw a readable error", "[DecodingQuantities]")
{
  const asmc::test::TempDirectory temp;
  const auto file = temp.file("decoding_quantities.txt");

  SECTION("a row with the wrong number of entries")
  {
    asmc::test::writeTextFile(file, "TransitionType\nSMC\nStates\n3\nExpectedTimes\n1.0 2.0\n");
    CHECK_THROWS_WITH(DecodingQuantities{file.string()},
                      ContainsSubstring("Parsed 2 ExpectedTimes entries, but expected 3"));
  }

  SECTION("a section before the number of states")
  {
    asmc::test::writeTextFile(file, "TransitionType\nSMC\nExpectedTimes\n1.0 2.0 3.0\n");
    CHECK_THROWS_WITH(DecodingQuantities{file.string()},
                      ContainsSubstring("Parsed ExpectedTimes before parsing states"));
  }
}
