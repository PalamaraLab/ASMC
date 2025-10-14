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

#ifndef ASMC_IBDPAIRDATALINE_HPP
#define ASMC_IBDPAIRDATALINE_HPP

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

struct IbdPairDataLine {

  std::string ind1FamId = "0_00";
  std::string ind1Id = "0_00";
  int ind1Hap = -1;

  std::string ind2FamId = "0_00";
  std::string ind2Id = "0_00";
  int ind2Hap = -1;

  int chromosome = -1;

  int ibdStart = -1;
  int ibdEnd = -1;

  float lengthInCentimorgans = -1.f;
  float ibdScore = -1.f;
  float postEst = -1.f;
  float mapEst = -1.f;

  [[nodiscard]] std::string toString() const
  {
    std::stringstream line;
    line << std::setprecision(std::numeric_limits<float>::digits10 + 1);

    line << ind1FamId << '\t' << ind1Id << '\t' << ind1Hap << '\t' << ind2FamId << '\t' << ind2Id << '\t' << ind2Hap
         << '\t' << chromosome << '\t' << ibdStart << '\t' << ibdEnd;

    if (lengthInCentimorgans != -1.f) {
      line << '\t' << lengthInCentimorgans;
    }

    if (ibdScore != -1.f) {
      line << '\t' << ibdScore;
    }

    if (postEst != -1.f) {
      line << '\t' << postEst;
    }

    if (mapEst != -1.f) {
      line << '\t' << mapEst;
    }

    return line.str();
  }
};

#endif // ASMC_IBDPAIRDATALINE_HPP
