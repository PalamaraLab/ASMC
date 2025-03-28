// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "PlinkMap.hpp"

#include "utils/FileUtils.hpp"
#include "utils/StringUtils.hpp"
#include "utils/VectorUtils.hpp"

#include <exception>
#include <iostream>

#include <fmt/core.h>
#include <fmt/ostream.h>

namespace asmc {

PlinkMap::PlinkMap(std::string_view mapFile) : mInputFile{mapFile} {
  validateFile();
  readFile();
  validateMap();
}

void PlinkMap::validateFile() {

  // Check file exists
  if (!fs::is_regular_file(mInputFile)) {
    throw std::runtime_error(fmt::format("Error: PLINK map file {} does not exist\n", mInputFile.string()));
  }

  // Check that the file contains either 3 or 4 tab-separated columns
  auto gzFile = gzopen(mInputFile.string().c_str(), "r");
  std::vector<std::string> firstLine = splitTextByDelimiter(readNextLineFromGzip(gzFile), "\t");
  mNumCols = static_cast<unsigned long>(firstLine.size());
  gzclose(gzFile);

  if (!(mNumCols == 3ul || mNumCols == 4ul)) {
    throw std::runtime_error(
        fmt::format("Error: PLINK map file {} should contain either 3 or 4 tab-separated columns, but contains {}\n",
                    mInputFile.string(), mNumCols));
  }

  // Count the number of lines in the file
  mNumSites = countLinesInFile(mInputFile);
}

void PlinkMap::readFile() {
  // Reserve space in vectors
  mChrIds.reserve(mNumSites);
  mSnpIds.reserve(mNumSites);
  if (mNumCols == 4ul) {
    mGeneticPositions.reserve(mNumSites);
  }
  mPhysicalPositions.reserve(mNumSites);

  const unsigned long chrCol = 0ul;
  const unsigned long snpCol = 1ul;
  const unsigned long genCol = 2ul;
  const unsigned long physCol = mNumCols == 4ul ? 3ul : 2ul;

  auto gzFile = gzopen(mInputFile.string().c_str(), "r");

  while (!gzeof(gzFile)) {
    std::vector<std::string> line = splitTextByDelimiter(readNextLineFromGzip(gzFile), "\t");
    if (!line.empty()) {

      if (line.size() != mNumCols) {
        gzclose(gzFile);
        throw std::runtime_error(
            fmt::format("Error: PLINK map file {} line {} contains {} columns, but line 1 contains {} columns\n",
                        mInputFile.string(), 1ul + mChrIds.size(), line.size(), mNumCols));
      }

      mChrIds.emplace_back(line.at(chrCol));
      mSnpIds.emplace_back(line.at(snpCol));
      if (mNumCols == 4ul) {
        try {
          mGeneticPositions.emplace_back(dblFromString(line.at(genCol)));
        } catch (const std::runtime_error& e) {
          gzclose(gzFile);
          throw std::runtime_error(fmt::format(
              "Error: PLINK map file {} line {} column {}: expected floating point but got {}\n{}\n",
              mInputFile.string(), 1ul + mGeneticPositions.size(), 1ul + genCol, line.at(physCol), e.what()));
        }
      }

      try {
        mPhysicalPositions.emplace_back(ulFromString(line.at(physCol)));
      } catch (const std::runtime_error& e) {
        gzclose(gzFile);
        throw std::runtime_error(fmt::format(
            "Error: PLINK map file {} line {} column {}: expected unsigned integer but got {}\n{}\n",
            mInputFile.string(), 1ul + mPhysicalPositions.size(), 1ul + physCol, line.at(physCol), e.what()));
      }
    }
  }

  gzclose(gzFile);
}

void PlinkMap::validateMap() {
  if (!isStrictlyIncreasing(mPhysicalPositions)) {
    fmt::print(std::cout, "Warning: PLINK map file {} physical positions are not strictly increasing\n", mInputFile.string());
    for (auto i = 1ul; i < mPhysicalPositions.size(); ++i) {
      if (mPhysicalPositions[i] <= mPhysicalPositions[i - 1]) {
        fmt::print(std::cout, "indices {} and {} have consecutive values {} and {}\n", i - 1, i, mPhysicalPositions[i - 1],
                   mPhysicalPositions[i]);
      }
    }
  }
  if (!isIncreasing(mGeneticPositions)) {
    fmt::print(std::cout, "Warning: PLINK map file {} genetic positions are not increasing\n", mInputFile.string());
    for (auto i = 1ul; i < mGeneticPositions.size(); ++i) {
      if (mGeneticPositions[i] < mGeneticPositions[i - 1]) {
        fmt::print(std::cout, "indices {} and {} have consecutive values {} and {}\n", i - 1, i, mGeneticPositions[i - 1],
                   mGeneticPositions[i]);
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

unsigned long PlinkMap::getNumSites() const {
  return mNumSites;
}

unsigned long PlinkMap::getNumCols() const {
  return mNumCols;
}

const std::vector<std::string>& PlinkMap::getChrIds() const {
  return mChrIds;
}

const std::vector<std::string>& PlinkMap::getSnpIds() const {
  return mSnpIds;
}

const std::vector<double>& PlinkMap::getGeneticPositions() const {
  return mGeneticPositions;
}

const std::vector<unsigned long>& PlinkMap::getPhysicalPositions() const {
  return mPhysicalPositions;
}

} // namespace asmc
