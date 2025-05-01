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

#include <sstream>
#include <string>
#include <vector>

#include "HmmUtils.hpp"
#include "Simd.hpp"

#include "Eigen/Core"

using Catch::Matchers::StartsWith;

TEST_CASE("test HMM utility free functions", "[HmmUtils]")
{

  SECTION("test subsetXorVec")
  {
    std::vector<bool> v1 = {false, false, true, true, false, false};
    std::vector<bool> v2 = {false, false, true, true, true, false};

    std::vector<bool> fullXor = {false, false, false, false, true, false};
    REQUIRE(asmc::subsetXorVec(v1, v2) == fullXor);
    REQUIRE(asmc::subsetXorVec(v1, v2, 0, 6) == fullXor);

    std::vector<bool> partialXor = {false, true};
    REQUIRE(asmc::subsetXorVec(v1, v2, 3, 5) == partialXor);
  }

  SECTION("test subsetAndVec")
  {
    std::vector<bool> v1 = {false, false, true, true, false, false};
    std::vector<bool> v2 = {false, false, true, true, true, false};

    std::vector<bool> fullAnd = {false, false, true, true, false, false};
    REQUIRE(asmc::subsetAndVec(v1, v2) == fullAnd);
    REQUIRE(asmc::subsetAndVec(v1, v2, 0, 6) == fullAnd);

    std::vector<bool> partialAnd = {true, false};
    REQUIRE(asmc::subsetAndVec(v1, v2, 3, 5) == partialAnd);
  }

  SECTION("test printVector")
  {
    std::vector<float> v = {1.23f, 2.34f, -3.45f};
    std::string answer = "\t1.23\t2.34\t-3.45\n";

    std::stringstream ss;
    asmc::printVector<float>(v, ss);

    REQUIRE(answer == ss.str());
  }

  SECTION("test printPctTime")
  {
    std::stringstream ss1;
    std::stringstream ss2;
    std::stringstream ss3;

    asmc::printPctTime("outputPerPair", 0.23678, ss1);
    asmc::printPctTime("combine", 1.0, ss2);
    asmc::printPctTime("other", 0.0, ss3);

    std::string answer1 = "Time in outputPerPair  :  23.7%\n";
    std::string answer2 = "Time in combine        : 100.0%\n";
    std::string answer3 = "Time in other          :   0.0%\n";

    REQUIRE(answer1 == ss1.str());
    REQUIRE(answer2 == ss2.str());
    REQUIRE(answer3 == ss3.str());
  }

  SECTION("test getSumOfVector")
  {
    std::vector<float> v1 = {1.23f, 2.34f, -3.45f};
    float answerF = 1.23f + 2.34f - 3.45f;
    REQUIRE(answerF == Approx(asmc::getSumOfVector(v1)));

    std::vector<long> v2 = {123l, 234l, -345l};
    long answerL = 123l + 234l - 345l;
    REQUIRE(answerL == asmc::getSumOfVector(v2));
  }

  SECTION("test elementWiseMultVectorScalar")
  {
    std::vector<float> v = {1.23f, 2.34f, -3.45f};
    float val = 4.56f;

    std::vector<float> answer = {1.23f * 4.56f, 2.34f * 4.56f, -3.45f * 4.56f};
    std::vector<float> calculated = asmc::elementWiseMultVectorScalar(v, val);
    REQUIRE(answer.size() == calculated.size());

    for (auto i = 0ul; i < v.size(); ++i) {
      REQUIRE(answer[i] == Approx(calculated[i]));
    }
  }

  SECTION("test elementWiseMultVectorVector")
  {
    std::vector<float> v1 = {1.23f, 2.34f, -3.45f};
    std::vector<float> v2 = {4.56f, -5.67f, 6.78f};

    std::vector<float> answer = {1.23f * 4.56f, 2.34f * -5.67f, -3.45f * 6.78f};
    std::vector<float> calculated = asmc::elementWiseMultVectorVector(v1, v2);
    REQUIRE(answer.size() == calculated.size());

    for (auto i = 0ul; i < v1.size(); ++i) {
      REQUIRE(answer[i] == Approx(calculated[i]));
    }
  }

  SECTION("test elementWiseMultMatrixMatrix")
  {
    std::vector<std::vector<float>> m1 = {{1.2f, 2.3f, -3.4f}, {4.5f, 5.6f, 6.7f}};
    std::vector<std::vector<float>> m2 = {{9.8f, -8.7f, 7.6f}, {6.5f, 5.4f, 4.3f}};

    std::vector<std::vector<float>> answer = {{1.2f * 9.8f, 2.3f * -8.7f, -3.4f * 7.6f},
                                              {4.5f * 6.5f, 5.6f * 5.4f, 6.7f * 4.3f}};
    std::vector<std::vector<float>> calculated = asmc::elementWiseMultMatrixMatrix(m1, m2);

    REQUIRE(answer.size() == calculated.size());
    REQUIRE(answer.at(0).size() == calculated.at(0).size());

    for (auto i = 0ul; i < m1.size(); ++i) {
      for (auto j = 0ul; j < m1[0].size(); ++j) {
        REQUIRE(answer[i][j] == Approx(calculated[i][j]));
      }
    }
  }

  SECTION("test normalizeMatrixColumns")
  {
    std::vector<std::vector<float>> mat = {{1.2f, 2.3f, -3.4f}, {4.5f, 5.6f, 6.7f}};
    std::vector<float> sums = {1.2f + 4.5f, 2.3f + 5.6f, -3.4f + 6.7f};

    std::vector<std::vector<float>> answer = {{1.2f / sums[0], 2.3f / sums[1], -3.4f / sums[2]},
                                              {4.5f / sums[0], 5.6f / sums[1], 6.7f / sums[2]}};
    std::vector<std::vector<float>> calculated = asmc::normalizeMatrixColumns(mat);

    REQUIRE(answer.size() == calculated.size());
    REQUIRE(answer.at(0).size() == calculated.at(0).size());

    for (auto i = 0ul; i < mat.size(); ++i) {
      for (auto j = 0ul; j < mat[0].size(); ++j) {
        REQUIRE(answer[i][j] == Approx(calculated[i][j]));
      }
    }
  }

  SECTION("test fillMatrixColumn")
  {
    std::vector<std::vector<float>> mat = {{0.f, 0.f, 0.f}, {0.f, 0.f, 0.f}};
    std::vector<float> newCol = {1.2f, -2.3f};

    std::vector<std::vector<float>> answer = {{0.f, 1.2f, 0.f}, {0.f, -2.3f, 0.f}};
    asmc::fillMatrixColumn(mat, newCol, 1ul);

    for (auto i = 0ul; i < mat.size(); ++i) {
      for (auto j = 0ul; j < mat[0].size(); ++j) {
        REQUIRE(answer[i][j] == Approx(mat[i][j]));
      }
    }
  }

  SECTION("test roundMorgans")
  {
    // Check smaller than minGenetic rounds to minGenetic
    {
      const float minGenetic = 0.5f;
      REQUIRE(asmc::roundMorgans(0.4f, 0, minGenetic) == minGenetic);
      REQUIRE(asmc::roundMorgans(0.4f, 1, minGenetic) == minGenetic);
      REQUIRE(asmc::roundMorgans(0.4f, 2, minGenetic) == minGenetic);
      REQUIRE(asmc::roundMorgans(0.f, 5, minGenetic) == minGenetic);
      REQUIRE(asmc::roundMorgans(-1.f, 7, minGenetic) == minGenetic);
    }

    // Check rounding as expected
    {
      const float minGenetic = 1e-10f;
      const float a = 0.123456f;
      REQUIRE(asmc::roundMorgans(a, 0, minGenetic) == 0.1f);
      REQUIRE(asmc::roundMorgans(a, 1, minGenetic) == 0.12f);
      REQUIRE(asmc::roundMorgans(a, 2, minGenetic) == 0.123f);
      REQUIRE(asmc::roundMorgans(a, 3, minGenetic) == 0.1235f);
      REQUIRE(asmc::roundMorgans(a, 4, minGenetic) == 0.12346f);
    }
  }

  SECTION("test roundPhysical")
  {
    // Check small values round to 1
    {
      REQUIRE(asmc::roundPhysical(-1, 0) == 1);
      REQUIRE(asmc::roundPhysical(-1, 1) == 1);
      REQUIRE(asmc::roundPhysical(-1, 2) == 1);
      REQUIRE(asmc::roundPhysical(0, 0) == 1);
      REQUIRE(asmc::roundPhysical(0, 1) == 1);
      REQUIRE(asmc::roundPhysical(0, 2) == 1);
      REQUIRE(asmc::roundPhysical(1, 0) == 1);
      REQUIRE(asmc::roundPhysical(1, 1) == 1);
      REQUIRE(asmc::roundPhysical(1, 2) == 1);
    }

    // Check rounding as expected
    {
      const int a = 123456;
      REQUIRE(asmc::roundPhysical(a, 0) == 100000);
      REQUIRE(asmc::roundPhysical(a, 1) == 120000);
      REQUIRE(asmc::roundPhysical(a, 2) == 123000);
      REQUIRE(asmc::roundPhysical(a, 3) == 123500);
      REQUIRE(asmc::roundPhysical(a, 4) == 123460);
      REQUIRE(asmc::roundPhysical(a, 5) == 123456);
    }
  }

  SECTION("test calculateScalingBatch")
  {
    // We need to make sure that this test will work independent of the number of SIMD lanes being used
    const int batchSize = std::max<int>(2, 2 * asmc::getNumSimdLanes());
    const int numStates = 2;

    Eigen::ArrayXf data(batchSize * numStates);
    Eigen::ArrayXf scalings(batchSize);
    Eigen::ArrayXf sums(batchSize);

    // set up data like 1, 2, 3, ..., 1, 2, 3, ..., ... so that scalings should be 1/(1+1), 1/(2+2), 1/(3+3) etc
    for (int stateIdx = 0; stateIdx < numStates; ++stateIdx) {
      for (int batchItem = 0; batchItem < batchSize; ++batchItem) {
        data[stateIdx * batchSize + batchItem] = 1.f + static_cast<float>(batchItem);
      }
    }

    // Calculate explicit answers by hand
    std::vector<float> answers(batchSize);
    for (auto i = 0; i < batchSize; ++i) {
      answers.at(i) = 1.f / (static_cast<float>(1ul + i) + static_cast<float>(1ul + i));
    }

    asmc::calculateScalingBatch(data, scalings, sums, batchSize, numStates);

    for (auto i = 0; i < batchSize; ++i) {
      REQUIRE(scalings[i] == Approx(answers.at(i)));
    }
  }

  SECTION("test applyScalingBatch")
  {
    // We need to make sure that this test will work independent of the number of SIMD lanes being used
    const int batchSize = std::max<int>(2, 2 * asmc::getNumSimdLanes());
    const int numStates = 2;

    Eigen::ArrayXf data(batchSize * numStates);
    Eigen::ArrayXf scalings(batchSize);

    // set up data like 1, 2, 3, ..., 1, 2, 3, ..., ...
    for (int stateIdx = 0; stateIdx < numStates; ++stateIdx) {
      for (int batchItem = 0; batchItem < batchSize; ++batchItem) {
        data[stateIdx * batchSize + batchItem] = 1.f + static_cast<float>(batchItem);
      }
    }

    // set up scalings to be 5, 6, 7, ...
    for (int batchItem = 0; batchItem < batchSize; ++batchItem) {
      scalings[batchItem] = static_cast<float>(batchItem) + 5.f;
    }

    // Calculate explicit answers by hand
    std::vector<float> answers(batchSize * numStates);
    for (int stateIdx = 0; stateIdx < numStates; ++stateIdx) {
      for (int batchItem = 0; batchItem < batchSize; ++batchItem) {
        answers.at(stateIdx * batchSize + batchItem) =
            (1.f + static_cast<float>(batchItem)) * (static_cast<float>(batchItem) + 5.f);
      }
    }

    asmc::applyScalingBatch(data, scalings, batchSize, numStates);

    for (auto i = 0; i < batchSize * numStates; ++i) {
      REQUIRE(data[i] == Approx(answers.at(i)));
    }
  }

  SECTION("test getFromPosition")
  {
    std::vector<float> geneticPositions = {0.12f, 0.23f, 0.34f, 0.45f, 0.56f, 0.67f};
    const unsigned from = 4u;  // the 0.56f element

    REQUIRE(asmc::getFromPosition(geneticPositions, from, 1.f) == 3u);  // 56 - 1 = 55, so 0.45f will be picked
    REQUIRE(asmc::getFromPosition(geneticPositions, from, 21.f) == 2u);  // 56 - 21 = 35, so 0.34f will be picked
    REQUIRE(asmc::getFromPosition(geneticPositions, from, 23.f) == 1u);  // 56 - 23 = 33, so 0.23f will be picked
    REQUIRE(asmc::getFromPosition(geneticPositions, from, 30.f) == 1u);  // 56 - 30 = 26, so 0.23f will be picked
    REQUIRE(asmc::getFromPosition(geneticPositions, from, 45.f) == 0u);  // 56 - 45 = 11, so 0.12f will be picked
    REQUIRE(asmc::getFromPosition(geneticPositions, from, 60.f) == 0u);  // anything bigger should return zero index

    // For regular ASMC, from will always be 0, and 0 should be returned no matter the distance
    REQUIRE(asmc::getFromPosition(geneticPositions, 0u, 1e-6f) == 0u);
    REQUIRE(asmc::getFromPosition(geneticPositions, 0u, 1.f) == 0u);
    REQUIRE(asmc::getFromPosition(geneticPositions, 0u, 10.f) == 0u);
  }

  SECTION("test getToPosition")
  {
    std::vector<float> geneticPositions = {0.12f, 0.23f, 0.34f, 0.45f, 0.56f, 0.67f};
    const unsigned to = 1u;  // the 0.23f element

    REQUIRE(asmc::getToPosition(geneticPositions, to, 1.f) == 3u);  // 23 + 1 = 24, so 0.34f will be picked (+1)
    REQUIRE(asmc::getToPosition(geneticPositions, to, 10.f) == 3u);  // 23 + 10 = 33, so 0.34f will be picked (+1)
    REQUIRE(asmc::getToPosition(geneticPositions, to, 12.f) == 4u);  // 23 + 12 = 35, so 0.45f will be picked (+1)
    REQUIRE(asmc::getToPosition(geneticPositions, to, 30.f) == 5u);  // 23 + 30 = 53, so 0.56f will be picked (+1)
    REQUIRE(asmc::getToPosition(geneticPositions, to, 40.f) == 6u);  // 23 + 40 = 63, so 0.67f will be picked (+1)
    REQUIRE(asmc::getToPosition(geneticPositions, to, 60.f) == 6u);  // anything bigger should return the vector size

    // For regular ASMC, `to` will always be the sequence length, and that should be returned no matter the distance
    REQUIRE(asmc::getToPosition(geneticPositions, 6u, 1e-6f) == 6u);
    REQUIRE(asmc::getToPosition(geneticPositions, 6u, 1.f) == 6u);
    REQUIRE(asmc::getToPosition(geneticPositions, 6u, 10.f) == 6u);
  }

  SECTION("test combinedIdToIndPlusHap")
  {
    CHECK_THROWS_WITH(asmc::combinedIdToIndPlusHap(""), StartsWith("Expected combined ID in form <id>_1"));
    CHECK_THROWS_WITH(asmc::combinedIdToIndPlusHap("abcde"), StartsWith("Expected combined ID in form <id>_1"));
    CHECK_THROWS_WITH(asmc::combinedIdToIndPlusHap("abcde_3"), StartsWith("Expected combined ID in form <id>_1"));
    CHECK_THROWS_WITH(asmc::combinedIdToIndPlusHap("_1"), StartsWith("Expected combined ID in form <id>_1"));

    CHECK(asmc::combinedIdToIndPlusHap("1_24_1") == std::make_pair<std::string, unsigned long>("1_24", 1ul));
    CHECK(asmc::combinedIdToIndPlusHap("1_2") == std::make_pair<std::string, unsigned long>("1", 2ul));
  }

  SECTION("test indPlusHapToCombinedId")
  {
    CHECK_THROWS_WITH(asmc::indPlusHapToCombinedId("", 1), StartsWith("Expected an individual ID and"));
    CHECK_THROWS_WITH(asmc::indPlusHapToCombinedId("123", 3), StartsWith("Expected an individual ID and"));
    CHECK_THROWS_WITH(asmc::indPlusHapToCombinedId("123", 0), StartsWith("Expected an individual ID and"));

    CHECK(asmc::indPlusHapToCombinedId("1_24", 1) == "1_24_1");
    CHECK(asmc::indPlusHapToCombinedId("1", 2) == "1_2");
  }
}
