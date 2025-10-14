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

#include "BinaryDataReader.hpp"

#include <fmt/format.h>

#include <filesystem>

namespace fs = std::filesystem;

BinaryDataReader::BinaryDataReader(const std::string& binaryFile)
{
  if (!fs::is_regular_file(binaryFile)) {
    throw std::runtime_error(fmt::format("Provided path to binary file {} is not a file\n", binaryFile));
  }

  mGzBinaryFileHandle = gzopen(binaryFile.c_str(), "rb");
  ReadHeader();
  CheckIfNextLineExists();
}

BinaryDataReader::~BinaryDataReader()
{
  gzclose(mGzBinaryFileHandle);
}

void BinaryDataReader::ReadHeader()
{
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mContainsIbdSegmentLengths), sizeof(bool));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mContainsIbdScore), sizeof(bool));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mContainsPosteriorAgeEstimates), sizeof(bool));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mContainsMapAgeEstimates), sizeof(bool));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mChromosomeNumber), sizeof(int));

  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mNumIds), sizeof(unsigned));
  mFamIds.reserve(mNumIds);
  mIIds.reserve(mNumIds);

  for (unsigned i = 0; i < mNumIds; i++) {

    unsigned lengthFamId = {};
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&lengthFamId), sizeof(unsigned));
    mFamIds.emplace_back(lengthFamId, 'z');
    gzread(mGzBinaryFileHandle, &mFamIds.at(i).at(0), lengthFamId);

    unsigned lengthIId = {};
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&lengthIId), sizeof(unsigned));
    mIIds.emplace_back(lengthIId, 'z');
    gzread(mGzBinaryFileHandle, &mIIds.at(i).at(0), lengthIId);
  }
}

void BinaryDataReader::CheckIfNextLineExists()
{
  if (gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mPreReadStartOfNextLine), sizeof(unsigned)) <
      sizeof(unsigned)) {
    mMoreLinesInFile = false;
  }
}

IbdPairDataLine BinaryDataReader::getNextLine()
{
  IbdPairDataLine line;

  // We have already read the first number from this line with a call to CheckIfNextLineExists()
  unsigned ind1 = mPreReadStartOfNextLine;
  unsigned ind2 = -1;

  std::uint_least8_t hap1;
  std::uint_least8_t hap2;

  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&hap1), sizeof(std::uint_least8_t));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&ind2), sizeof(unsigned));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&hap2), sizeof(std::uint_least8_t));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.ibdStart), sizeof(int));
  gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.ibdEnd), sizeof(int));

  if (mContainsIbdSegmentLengths) {
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.lengthInCentimorgans), sizeof(float));
  }

  if (mContainsIbdScore) {
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.ibdScore), sizeof(float));
  }

  if (mContainsPosteriorAgeEstimates) {
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.postEst), sizeof(float));
  }

  if (mContainsMapAgeEstimates) {
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.mapEst), sizeof(float));
  }

  line.ind1Hap = static_cast<int>(hap1);
  line.ind2Hap = static_cast<int>(hap2);

  line.chromosome = mChromosomeNumber;

  line.ind1FamId = mFamIds.at(ind1);
  line.ind1Id = mIIds.at(ind1);

  line.ind2FamId = mFamIds.at(ind2);
  line.ind2Id = mIIds.at(ind2);

  // Pre-read first number from next line to check whether the line exists
  CheckIfNextLineExists();

  return line;
}

bool BinaryDataReader::moreLinesInFile() const
{
  return mMoreLinesInFile;
}
