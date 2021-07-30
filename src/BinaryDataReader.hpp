//
// Created by fergus on 28/08/2020.
//

#ifndef ASMC_BINARYDATAREADER_HPP
#define ASMC_BINARYDATAREADER_HPP

#include <zlib.h>

#include <array>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

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

  std::string toString()
  {
    std::stringstream line;
    line << std::setprecision(std::numeric_limits<float>::digits10 + 1);

    line << ind1FamId << '\t' << ind1Id << '\t' << ind1Hap << '\t' << ind2FamId << '\t' << ind2Id << '\t' << ind2Hap
         << '\t' << chromosome << '\t' << ibdStart << '\t' << ibdEnd;

    if (lengthInCentimorgans != -1.f) {
      line << '\t' << lengthInCentimorgans;
    }

    line << '\t' << ibdScore;

    if (postEst != -1.f) {
      line << '\t' << postEst;
    }

    if (mapEst != -1.f) {
      line << '\t' << mapEst;
    }

    return line.str();
  }
};

class BinaryDataReader
{

private:
  /**
   * Handle to the binary zipped file, opened in the constructor and closed in the destructor
   */
  gzFile mGzBinaryFileHandle;

  bool mContainsIbdSegmentLengths = false;
  bool mContainsPosteriorAgeEstimates = false;
  bool mContainsMapAgeEstimates = false;

  int mChromosomeNumber = -1;
  unsigned mNumIds = 0u;

  std::vector<std::string> mFamIds;
  std::vector<std::string> mIIds;

  unsigned mPreReadStartOfNextLine = {};

  bool mMoreLinesInFile = true;

  void ReadHeader()
  {
    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mContainsIbdSegmentLengths), sizeof(bool));
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

  void CheckIfNextLineExists()
  {
    if (gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&mPreReadStartOfNextLine), sizeof(unsigned)) <
        sizeof(unsigned)) {
      mMoreLinesInFile = false;
    }
  }

public:
  explicit BinaryDataReader(const std::string& binaryFile)
  {
    mGzBinaryFileHandle = gzopen(binaryFile.c_str(), "rb");
    ReadHeader();
    CheckIfNextLineExists();
  }

  IbdPairDataLine getNextLine()
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

    gzread(mGzBinaryFileHandle, reinterpret_cast<char*>(&line.ibdScore), sizeof(float));

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

  bool moreLinesInFile() const
  {
    return mMoreLinesInFile;
  }

  ~BinaryDataReader()
  {
    gzclose(mGzBinaryFileHandle);
  }
};

#endif // ASMC_BINARYDATAREADER_HPP
