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

#ifndef ASMC_BINARYDATAREADER_HPP
#define ASMC_BINARYDATAREADER_HPP

#include "IbdPairDataLine.hpp"

#include <zlib.h>

#include <string>
#include <vector>

class BinaryDataReader
{

private:
  /**
   * Handle to the binary zipped file, opened in the constructor and closed in the destructor
   */
  gzFile mGzBinaryFileHandle;

  bool mContainsIbdSegmentLengths = false;
  bool mContainsIbdScore = false;
  bool mContainsPosteriorAgeEstimates = false;
  bool mContainsMapAgeEstimates = false;

  int mChromosomeNumber = -1;
  unsigned mNumIds = 0u;

  std::vector<std::string> mFamIds;
  std::vector<std::string> mIIds;

  unsigned mPreReadStartOfNextLine = {};

  bool mMoreLinesInFile = true;

  void ReadHeader();

  void CheckIfNextLineExists();

public:
  explicit BinaryDataReader(const std::string& binaryFile);

  IbdPairDataLine getNextLine();

  [[nodiscard]] bool moreLinesInFile() const;

  ~BinaryDataReader();
};

#endif // ASMC_BINARYDATAREADER_HPP
