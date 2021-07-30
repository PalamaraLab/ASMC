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

#include "hashing/ExtendHash.hpp"
#include "hashing/Individuals.hpp"
#include "hashing/Match.hpp"
#include "hashing/SeedHash.hpp"
#include "hashing/Utils.hpp"

TEST_CASE("ExtendHash", "[hashing]")
{
  ExtendHash e(4ul, 2ul, true);
  REQUIRE(e.size() == 0ul);
  REQUIRE(e.getWordSize() == 4ul);

  //todo: test ExtendHash
}

TEST_CASE("individuals", "[hashing]")
{
  Individuals ind(8ul, 3ul, 5u);

  REQUIRE(ind.getIdNum() == 5u);
  REQUIRE(ind.getWordSize() == 8ul);
  REQUIRE(ind.getNumReadAhead() == 3ul);

  // Check up to 10 - but internally we're just going 0-1-2-0-1-2-0-1-2-0
  for (auto i = 0; i < 10; ++i) {
    REQUIRE(ind.getWordHash(i) == 0ul);
    REQUIRE(ind.getWordString(i) == "00000000");
  }

  ind.setMarker(0, 0);
  ind.setMarker(1, 2);

  ind.setMarker(2, 2);
  ind.setMarker(2, 3);

  REQUIRE(ind.getWordHash(0) == 1ul);
  REQUIRE(ind.getWordString(0) == "00000001");

  REQUIRE(ind.getWordHash(1) == 4ul);
  REQUIRE(ind.getWordString(1) == "00000100");

  REQUIRE(ind.getWordHash(2) == 12ul);
  REQUIRE(ind.getWordString(2) == "00001100");

  // Clear 2
  ind.clear(2);
  REQUIRE(ind.getWordHash(2) == 0ul);
  REQUIRE(ind.getWordString(2) == "00000000");

  // Clear 1 by clearing 4
  ind.clear(4);
  REQUIRE(ind.getWordHash(1) == 0ul);
  REQUIRE(ind.getWordString(1) == "00000000");
}

TEST_CASE("match", "[hashing]")
{
  SECTION("default construction")
  {
    Match m(4ul);
    REQUIRE(m.getWordSize() == 4ul);
    REQUIRE(m.getGaps() == 0u);
    REQUIRE(m.getInterval()[0] == 0);
    REQUIRE(m.getInterval()[1] == 0);

    m.addGap();
    m.addGap();
    REQUIRE(m.getGaps() == 2u);

    m.extend(5);
    REQUIRE(m.getInterval()[1] == 5);
  }

  SECTION("explicit constructor")
  {
    Match m(4, 7);
    REQUIRE(m.getWordSize() == 4ul);
    REQUIRE(m.getGaps() == 0u);
    REQUIRE(m.getInterval()[0] == 7);
    REQUIRE(m.getInterval()[1] == 7);

    m.extend(5);
    REQUIRE(m.getInterval()[1] == 7);

    m.extend(8);
    REQUIRE(m.getInterval()[1] == 8);
  }

  SECTION("print method")
  {
    //TODO: this method is harder to test because it requires access to an HMM instance
  }
}

TEST_CASE("SeedHash", "[hashing]")
{
  SeedHash s;
  REQUIRE(s.size() == 0ul);

  //todo: test SeedHash
}

TEST_CASE("utils", "[hashing]")
{
  SECTION("cmBetween")
  {
    std::vector<float> genPos = {0.00402186f, 0.0388124f, 0.0567817f, 0.0668489f, 0.0915063f, 0.12783f,  0.198618f,
                                 0.199045f,   0.250093f,  0.259338f,  0.293267f,  0.294899f,  0.316173f, 0.353332f,
                                 0.354553f,   0.357123f,  0.359118f,  0.395468f,  0.41749f,   0.421739f, 0.453347f,
                                 0.471302f,   0.535031f,  0.548733f,  0.574022f,  0.604538f,  0.620419f};

    SECTION("both words are inside vector")
    {
      const int wordSize = 4;
      const int w1 = 0;
      const int w2 = 3;
      const int w3 = 5;

      REQUIRE(asmc::cmBetween(w1, w2, genPos, wordSize) == 100.0 * (genPos.at(15) - genPos.at(0)));
      REQUIRE(asmc::cmBetween(w1, w3, genPos, wordSize) == 100.0 * (genPos.at(23) - genPos.at(0)));
      REQUIRE(asmc::cmBetween(w2, w3, genPos, wordSize) == 100.0 * (genPos.at(23) - genPos.at(12)));
    }

    SECTION("second word overflows vector")
    {
      const int wordSize = 4;
      const int w1 = 0;
      const int w2 = 1;
      const int w3 = 10;

      REQUIRE(asmc::cmBetween(w1, w3, genPos, wordSize) == 100.0 * (genPos.back() - genPos.at(0)));
      REQUIRE(asmc::cmBetween(w2, w3, genPos, wordSize) == 100.0 * (genPos.back() - genPos.at(4)));
    }
  }
}