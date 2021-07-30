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

#include "BinaryDataReader.hpp"

TEST_CASE("IbdPairDataLine default member test", "[BinaryDataReader]")
{
  IbdPairDataLine line;

  REQUIRE(line.ind1FamId == "0_00");
  REQUIRE(line.ind1Id == "0_00");
  REQUIRE(line.ind1Hap == -1);
  REQUIRE(line.ind2FamId == "0_00");
  REQUIRE(line.ind2Id == "0_00");
  REQUIRE(line.ind2Hap == -1);
  REQUIRE(line.chromosome == -1);
  REQUIRE(line.ibdStart == -1);
  REQUIRE(line.ibdEnd == -1);
  REQUIRE(line.lengthInCentimorgans == -1.f);
  REQUIRE(line.ibdScore == -1.f);
  REQUIRE(line.postEst == -1.f);
  REQUIRE(line.mapEst == -1.f);

  REQUIRE(line.toString() == "0_00\t0_00\t-1\t0_00\t0_00\t-1\t-1\t-1\t-1\t-1");

  line.lengthInCentimorgans = 1.2;
  line.postEst = 2.3;
  line.mapEst = 3.4;

  REQUIRE(line.toString() == "0_00\t0_00\t-1\t0_00\t0_00\t-1\t-1\t-1\t-1\t1.2\t-1\t2.3\t3.4");
}

TEST_CASE("BinaryDataReader real data test", "[BinaryDataReader]")
{
  BinaryDataReader dataReader(ASMC_TEST_DIR "/data/binary_output.bibd.gz");

  IbdPairDataLine line1 = dataReader.getNextLine();
  REQUIRE(line1.ind1FamId == "1_94");
  REQUIRE(line1.ind1Id == "1_94");
  REQUIRE(line1.ind1Hap == 1);
  REQUIRE(line1.ind2FamId == "1_104");
  REQUIRE(line1.ind2Id == "1_104");
  REQUIRE(line1.ind2Hap == 1);
  REQUIRE(line1.chromosome == 1);
  REQUIRE(line1.ibdStart == 8740);
  REQUIRE(line1.ibdEnd == 1660011);
  REQUIRE(line1.lengthInCentimorgans == Approx(1.86962f).epsilon(1e-5));
  REQUIRE(line1.ibdScore == Approx(0.403475f).epsilon(1e-5));
  REQUIRE(line1.postEst == Approx(146.203f).epsilon(1e-5));
  REQUIRE(line1.mapEst == Approx(24.9999f).epsilon(1e-5));

  IbdPairDataLine line2 = dataReader.getNextLine();
  REQUIRE(line2.ind1FamId == "1_94");
  REQUIRE(line2.ind1Id == "1_94");
  REQUIRE(line2.ind1Hap == 1);
  REQUIRE(line2.ind2FamId == "1_104");
  REQUIRE(line2.ind2Id == "1_104");
  REQUIRE(line2.ind2Hap == 1);
  REQUIRE(line2.chromosome == 1);
  REQUIRE(line2.ibdStart == 1679626);
  REQUIRE(line2.ibdEnd == 1679626);
  REQUIRE(line2.lengthInCentimorgans == Approx(0.f).epsilon(1e-5));
  REQUIRE(line2.ibdScore == Approx(0.0175673f).epsilon(1e-5));
  REQUIRE(line2.postEst == Approx(18029.8f).epsilon(1e-5));
  REQUIRE(line2.mapEst == Approx(24.9999f).epsilon(1e-5));

  int numLinesRead = 2;
  while(dataReader.moreLinesInFile()) {
    IbdPairDataLine line = dataReader.getNextLine();
    numLinesRead++;
  }

  REQUIRE(numLinesRead == 1520);
}
