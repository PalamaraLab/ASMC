// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef DATA_MODULE_GENETIC_MAP_HPP
#define DATA_MODULE_GENETIC_MAP_HPP

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>

namespace asmc {

namespace fs = std::filesystem;

/**
 * A class that reads and stores a genetic map.
 */
class GeneticMap {

private:

  /** Path to the input file */
  fs::path mInputFile;

  /** Whether the map has a header (determined automatically) */
  unsigned long mHasHeader = false;

  /** Number of sites/SNPs (#rows in the map file) */
  unsigned long mNumSites{};

  /** Number of columns in the map file: a valid map must contain at least 3 */
  unsigned long mNumCols{};

  /** The genetic positions in column three of the map. These are either in Morgans or Centimorgans */
  std::vector<double> mGeneticPositions;

  /** The physical positions in column one of the map */
  std::vector<unsigned long> mPhysicalPositions;

  /**
   * Check that:
   * - the file exists
   * - there is at least one row of data in addition to an optional header row
   * - the first (non-header) row contains at least 3 or four tab-separated columns
   * and count the number of rows in the file
   */
  void validateFile();

  /**
   * Check that a row from the map file:
   * - contains at least three tab-separated columns
   * - the first column contains an unsigned integer (physical position)
   *
   * @param row the row to validate (a string)
   * @return whether the row is valid
   */
  static bool validDataRow(const std::string& row);

  /**
   * Read each line into the relevant vectors, checking each line has the same number of columns and that the column
   * containing physical positions contains positive integer values.
   */
  void readFile();

  /**
   * Once the map has been read from file, validate that the genetic and physical positions are strictly increasing.
   * A runtime error will be thrown if this is not the case.
   */
  void validateMap();

public:
  /**
   * Read a genetic .map file: a tab-separated text file with or without a header header row, and one line per variant
   * with the following fields:
   * 1. physical position in base pairs (integer)
   * 2. \todo: what is column 2?
   * 3. genetic position in either Morgans or Centimorgans; ASMC assumes Centimorgans (float)
   * 4. other unspecified data columns
   * All lines must have the same number of columns.
   *
   * @param mapFile path to the .map file
   */
  explicit GeneticMap(std::string_view mapFile);

  [[nodiscard]] unsigned long getNumSites() const;
  [[nodiscard]] unsigned long getNumCols() const;
  [[nodiscard]] unsigned long hasHeader() const;
  [[nodiscard]] const std::vector<double>& getGeneticPositions() const;
  [[nodiscard]] const std::vector<unsigned long>& getPhysicalPositions() const;
};

} // namespace asmc

#endif // DATA_MODULE_GENETIC_MAP_HPP
