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


// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <limits>
#include <type_traits>

#include "StringUtils.hpp"

TEST_CASE("test string conversions", "[StringUtils]") {

  std::string str_1_0 = "1.0";
  std::string str_1e0 = "1.0E0";
  std::string str_minus1_0 = "-1.0";

  // Not representable as float or double
  std::string too_small_str = "3.20676899524985E-310";
  long double too_small_ld = 3.20676899524985E-310L;

  // Acceptable non-numbers
  std::string nan = "NAN";
  std::string inf = "INF";

  // Not representable at all
  std::string way_too_small = "1.23E-1000000000";

  // Not convertible
  std::string greeting = "hello";

  // Check stof
  CHECK(StringUtils::stof(str_1_0) == 1.f);
  CHECK(StringUtils::stof(str_1e0) == 1.f);
  CHECK(StringUtils::stof(str_minus1_0) == -1.f);
  CHECK(StringUtils::stof(too_small_str) == static_cast<float>(too_small_ld));
  CHECK(StringUtils::stof(way_too_small) == 0.f);
  CHECK(std::isnan(StringUtils::stof(nan)));
  CHECK(std::isinf(StringUtils::stof(inf)));
  CHECK_THROWS_AS(StringUtils::stof(greeting), std::invalid_argument);

  // Check stod
  CHECK(StringUtils::stod(str_1_0) == 1.0);
  CHECK(StringUtils::stod(str_1e0) == 1.0);
  CHECK(StringUtils::stod(str_minus1_0) == -1.0);
  CHECK(StringUtils::stod(too_small_str) == static_cast<double>(too_small_ld));
  CHECK(std::isnan(StringUtils::stod(nan)));
  CHECK(std::isinf(StringUtils::stod(inf)));
  CHECK(StringUtils::stod(way_too_small) == 0.0);
  CHECK_THROWS_AS(StringUtils::stod(greeting), std::invalid_argument);
}

TEST_CASE("test SOME coverage", "[coverage]") {

    auto s = StringUtils::findDelimiters("this;string;has;five;semi;colons", ";");

    CHECK(s == ";;;;;");
}


