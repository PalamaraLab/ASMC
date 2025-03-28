// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef DATA_MODULE_PLINK_MAP_HPP
#define DATA_MODULE_PLINK_MAP_HPP

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>

namespace asmc {

namespace fs = std::filesystem;

/**
 * A class that reads and stores a PLINK map.
 */
class PlinkMap {

private:

  /** Path to the input file */
  fs::path mInputFile;

  /** Number of sites/SNPs (#rows in the map file) */
  unsigned long mNumSites{};

  /** Number of columns in the map file: a valid map can contain either 3 or 4 */
  unsigned long mNumCols{};

  /** The chromosome IDs in the map. These need not necessarily be numeric IDs */
  std::vector<std::string> mChrIds;

  /** The site/SNP IDs. These are string variables */
  std::vector<std::string> mSnpIds;

  /** The genetic positions in the map. These are either in Morgans or Centimorgans, but are optional in the file */
  std::vector<double> mGeneticPositions;

  /** The physical positions */
  std::vector<unsigned long> mPhysicalPositions;

  /**
   * Check that:
   * - the file exists
   * - the first row contains either 3 or four tab-separated columns
   * and count the number of rows in the file
   */
  void validateFile();

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
   * Read a PLINK .map file: a text file with no header file, and one line per variant with the following 3-4 fields:
   * 1. chromosome code / ID (string)
   * 2. variant/SNP identifier (string)
   * 3. (optional) genetic position in either Morgans or Centimorgans; ASMC assumes Centimorgans (float)
   * 4. physical position in base pairs (integer)
   * All lines must have the same number of columns.
   *
   * @param mapFile path to the .map file
   */
  explicit PlinkMap(std::string_view mapFile);

  [[nodiscard]] unsigned long getNumSites() const;
  [[nodiscard]] unsigned long getNumCols() const;
  [[nodiscard]] const std::vector<std::string>& getChrIds() const;
  [[nodiscard]] const std::vector<std::string>& getSnpIds() const;
  [[nodiscard]] const std::vector<double>& getGeneticPositions() const;
  [[nodiscard]] const std::vector<unsigned long>& getPhysicalPositions() const;
};

} // namespace asmc

#endif // DATA_MODULE_PLINK_MAP_HPP
