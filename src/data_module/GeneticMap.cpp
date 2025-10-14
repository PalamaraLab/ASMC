// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "GeneticMap.hpp"

#include "utils/FileUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/VectorUtils.hpp"

#include <exception>
#include <iostream>

#include <fmt/core.h>
#include <fmt/ostream.h>

namespace asmc {

GeneticMap::GeneticMap(std::string_view mapFile) : mInputFile{mapFile} {
  validateFile();
  readFile();
  validateMap();
}

void GeneticMap::validateFile() {

  // Check file exists
  if (!fs::is_regular_file(mInputFile)) {
    throw std::runtime_error(fmt::format("Error: genetic map file {} does not exist\n", mInputFile.string()));
  }

  // Read (at most) two lines from the file
  std::vector<std::string> firstLines;
  {
    auto gzFile = gzopen(mInputFile.string().c_str(), "r");

    unsigned numLinesRead = 0u;
    while (!gzeof(gzFile) && numLinesRead < 2) {
      firstLines.emplace_back(readNextLineFromGzip(gzFile));
      numLinesRead++;
    }
    gzclose(gzFile);
  }

  // Determine whether the first two lines are valid
  std::vector<bool> validLines = {false, false};
  for (auto rowId = 0ul; rowId < firstLines.size(); ++rowId) {
    validLines.at(rowId) = validDataRow(firstLines.at(rowId));
  }

  bool validFile = false;
  bool potentialHeader = !validLines.at(0) && !firstLines.at(0).empty();

  if (validLines.at(0) && validLines.at(1)) { // both lines valid
    validFile = true;
  } else if (potentialHeader && validLines.at(1)) { // potential header followed by valid row
    mHasHeader = true;
    validFile = true;
  } else if (validLines.at(0) && !validLines.at(1)) { // only one line that's valid, or valid followed by empty
    validFile = firstLines.size() == 1ul || (firstLines.size() == 2ul && firstLines.back().empty());
  }

  if (!validFile) {
    throw std::runtime_error(fmt::format("Error: genetic map file {} should contain at least one data row with at "
                                         "least 3 tab-separated columns, but contains\n{}\n",
                                         mInputFile.string(), fmt::join(firstLines.begin(), firstLines.end(), "\n")));
  }

  mNumCols = static_cast<unsigned long>(validLines.at(0) ? splitTextByDelimiter(firstLines.at(0), "\t").size()
                                                         : splitTextByDelimiter(firstLines.at(1), "\t").size());

  mNumSites = countLinesInFile(mInputFile);
  if (mHasHeader) {
    mNumSites -= 1ul;
  }
}

bool GeneticMap::validDataRow(const std::string& row) {

  // an empty row isn't valid
  if (row.empty()) {
    return false;
  }

  // a row must contain at least 3 tab-separated values
  std::vector<std::string> rowParts = splitTextByDelimiter(row, "\t");
  if (rowParts.size() < 3ul) {
    return false;
  }

  // the first column of a valid row contains an unsigned integer (physical position)
  try {
    ulFromString(rowParts.front());
  } catch (const std::runtime_error&) {
    return false;
  }

  // the third column of a valid row contains a floating point value (genetic position)
  try {
    dblFromString(rowParts.at(2ul));
  } catch (const std::runtime_error&) {
    return false;
  }

  return true;
}

void GeneticMap::readFile() {

  mGeneticPositions.reserve(mNumSites);
  mPhysicalPositions.reserve(mNumSites);

  auto gzFile = gzopen(mInputFile.string().c_str(), "r");

  // Skip header, if present
  if (mHasHeader) {
    readNextLineFromGzip(gzFile);
  }

  while (!gzeof(gzFile)) {
    std::vector<std::string> line = splitTextByDelimiter(readNextLineFromGzip(gzFile), "\t");
    if (!line.empty()) {

      if (line.size() != mNumCols) {
        gzclose(gzFile);
        throw std::runtime_error(
            fmt::format("Error: Genetic map file {} line {} contains {} columns, but the first data row contains {}\n",
                        mInputFile.string(), 1ul + mGeneticPositions.size(), line.size(), mNumCols));
      }

      try {
        mPhysicalPositions.emplace_back(ulFromString(line.at(0ul)));
        mGeneticPositions.emplace_back(dblFromString(line.at(2ul)));
      } catch (const std::runtime_error&) {
        gzclose(gzFile);
        throw std::runtime_error(fmt::format(
            "Error: Genetic map file {} line {} should contain an unsigned integer physical position in the first "
            "column and a floating point genetic position in the third column, but found {} and {}\n",
            mInputFile.string(), 1ul + mGeneticPositions.size(), line.at(0ul), line.at(2ul)));
      }
    }
  }

  gzclose(gzFile);
}

void GeneticMap::validateMap() {
  if (!isStrictlyIncreasing(mPhysicalPositions)) {
    fmt::print(std::cout, "Warning: genetic map file {} physical positions are not strictly increasing\n",
               mInputFile.string());
    for (auto i = 1ul; i < mPhysicalPositions.size(); ++i) {
      if (mPhysicalPositions[i] <= mPhysicalPositions[i - 1]) {
        fmt::print(std::cout, "indices {} and {} have consecutive values {} and {}\n", i - 1, i,
                   mPhysicalPositions[i - 1], mPhysicalPositions[i]);
      }
    }
  }
  if (!isIncreasing(mGeneticPositions)) {
    fmt::print(std::cout, "Warning: genetic map file {} genetic positions are not increasing\n", mInputFile.string());
    for (auto i = 1ul; i < mGeneticPositions.size(); ++i) {
      if (mGeneticPositions[i] < mGeneticPositions[i - 1]) {
        fmt::print(std::cout, "indices {} and {} have consecutive values {} and {}\n", i - 1, i,
                   mGeneticPositions[i - 1], mGeneticPositions[i]);
      }
    }
  }

  // Check whether genetic positions are in an appropriate range to be in Centimorgans rather than Morgans.
  // Should be very roughly 1cM <-> 10^6 base pairs (check 90% are within 0.4 - 2.5 mega base pairs per Centimorgan
  unsigned long numOutOfRange = 0ul;
  const double upper = 1e6 / 0.4;
  const double lower = 1e6 / 2.5;
  for (auto i = 0ul; i < mGeneticPositions.size(); ++i) {
    if (upper * mGeneticPositions.at(i) < static_cast<double>(mPhysicalPositions.at(i)) ||
        lower * mGeneticPositions.at(i) > static_cast<double>(mPhysicalPositions.at(i))) {
      numOutOfRange++;
    }
  }
  if (const double ratio = static_cast<double>(numOutOfRange) / static_cast<double>(mGeneticPositions.size());
      ratio > 0.1) {
    fmt::print(std::cout,
               "Warning: {:.1f}% of entries in the genetic map file {} are not in the expected range for a human "
               "genome (0.4-2.5 mega base pairs per Centimorgan). Please check that your map file is providing genetic "
               "positions in Centimorgans.\n",
               100.0 * ratio, mInputFile.string());
  }
}

unsigned long GeneticMap::getNumSites() const {
  return mNumSites;
}

unsigned long GeneticMap::getNumCols() const {
  return mNumCols;
}

const std::vector<double>& GeneticMap::getGeneticPositions() const {
  return mGeneticPositions;
}

const std::vector<unsigned long>& GeneticMap::getPhysicalPositions() const {
  return mPhysicalPositions;
}

unsigned long GeneticMap::hasHeader() const {
  return mHasHeader;
}

} // namespace asmc
