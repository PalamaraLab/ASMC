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

// Helpers for regression tests, which compare output with snapshots stored in ASMC_data/testing/regression.
//
// To rewrite every snapshot from the current code, for example after a deliberate change in output, run the
// regression tests with ASMC_REGENERATE_SNAPSHOTS=1:
//
//   ASMC_REGENERATE_SNAPSHOTS=1 ./ASMC_regression
//
// then review the changes in the ASMC_data submodule and commit them there.

#ifndef ASMC_REGRESSION_SNAPSHOTS_HPP
#define ASMC_REGRESSION_SNAPSHOTS_HPP

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "FileUtils.hpp"

namespace asmc::test
{

/// Path of a snapshot file, given its name.
inline std::string snapshotPath(const std::string& name)
{
  return std::string(ASMC_DATA_DIR "/testing/regression/") + name;
}

/// Whether the tests should rewrite snapshots rather than compare with them.
inline bool regeneratingSnapshots()
{
  const char* value = std::getenv("ASMC_REGENERATE_SNAPSHOTS");
  return value != nullptr && std::string(value) == "1";
}

/// How far a value may differ from its snapshot: |actual - expected| <= absolute + relative * |expected|.
struct Tolerance {
  double relative = 0.0;
  double absolute = 0.0;
};

/// For integer output, such as MAP states and positions.
constexpr Tolerance exactMatch{};

/// For floating-point output. Results differ slightly between SIMD targets (e.g. through fused multiply-add), so
/// snapshots made on one machine must still pass on another.
constexpr Tolerance floatOutput{1e-3, 1e-6};

/// Read every line of a (possibly gzipped) text file.
inline std::vector<std::string> readLines(const std::string& path)
{
  std::vector<std::string> lines;
  FileUtils::AutoGzIfstream in;
  in.openOrExit(path);
  for (std::string line; getline(in, line);) {
    lines.push_back(line);
  }
  in.close();
  return lines;
}

/// Write lines to a gzipped text file.
inline void writeLines(const std::vector<std::string>& lines, const std::string& path)
{
  FileUtils::AutoGzOfstream out;
  out.openOrExit(path);
  for (const auto& line : lines) {
    out << line << '\n';
  }
  out.close();
}

/// Split a line into whitespace-separated fields.
inline std::vector<std::string> splitFields(const std::string& line)
{
  std::vector<std::string> fields;
  std::istringstream stream(line);
  for (std::string field; stream >> field;) {
    fields.push_back(field);
  }
  return fields;
}

/// Parse a field as a number, if the whole field is one.
inline bool parseNumber(const std::string& field, double& value)
{
  std::istringstream stream(field);
  stream >> value;
  return !stream.fail() && stream.eof();
}

/// Whether a field is written as an integer, so that it must match exactly.
inline bool looksLikeInteger(const std::string& field)
{
  return field.find_first_of(".eEnNiI") == std::string::npos;
}

/**
 * Whether two fields match: exactly if either is not a number or both are integers, otherwise within `tolerance`.
 */
inline bool fieldsMatch(const std::string& expected, const std::string& actual, const Tolerance tolerance)
{
  double expectedValue = 0.0;
  double actualValue = 0.0;
  const bool bothNumbers = parseNumber(expected, expectedValue) && parseNumber(actual, actualValue);
  if (!bothNumbers || (looksLikeInteger(expected) && looksLikeInteger(actual))) {
    return expected == actual;
  }
  if (std::isnan(expectedValue) || std::isnan(actualValue)) {
    return std::isnan(expectedValue) && std::isnan(actualValue);
  }
  return std::abs(actualValue - expectedValue) <= tolerance.absolute + (tolerance.relative * std::abs(expectedValue));
}

/**
 * Compare lines of output with the lines of snapshot `snapshotName`, field by field (see fieldsMatch), and record one
 * Catch2 check for the whole comparison. When regenerating, write `lines` as the new snapshot instead.
 */
inline void checkLinesMatchSnapshot(const std::vector<std::string>& lines, const std::string& snapshotName,
                                    const Tolerance tolerance)
{
  if (regeneratingSnapshots()) {
    std::filesystem::create_directories(std::filesystem::path(snapshotPath(snapshotName)).parent_path());
    writeLines(lines, snapshotPath(snapshotName));
    WARN("Regenerated snapshot " << snapshotName << " (" << lines.size() << " lines)");
    return;
  }

  const std::vector<std::string> expected = readLines(snapshotPath(snapshotName));
  INFO("snapshot " << snapshotName);
  REQUIRE(lines.size() == expected.size());

  std::size_t mismatchedLines = 0;
  std::string firstMismatch;
  for (std::size_t lineIdx = 0; lineIdx < expected.size(); ++lineIdx) {
    const auto expectedFields = splitFields(expected.at(lineIdx));
    const auto actualFields = splitFields(lines.at(lineIdx));
    bool match = expectedFields.size() == actualFields.size();
    for (std::size_t field = 0; match && field < expectedFields.size(); ++field) {
      match = fieldsMatch(expectedFields.at(field), actualFields.at(field), tolerance);
    }
    if (!match) {
      ++mismatchedLines;
      if (firstMismatch.empty()) {
        firstMismatch = "line " + std::to_string(lineIdx + 1) + ": expected \"" + expected.at(lineIdx) +
                        "\" but got \"" + lines.at(lineIdx) + "\"";
      }
    }
  }
  INFO(firstMismatch);
  CHECK(mismatchedLines == 0);
}

/// Format each row of a matrix as one line of space-separated values, precise enough to round-trip a float.
template <typename Derived> std::vector<std::string> matrixToLines(const Eigen::DenseBase<Derived>& matrix)
{
  std::vector<std::string> lines;
  lines.reserve(static_cast<std::size_t>(matrix.rows()));
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    std::ostringstream line;
    line << std::setprecision(std::numeric_limits<float>::max_digits10);
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      line << (col == 0 ? "" : " ") << matrix(row, col);
    }
    lines.push_back(line.str());
  }
  return lines;
}

/// Compare a matrix with snapshot `snapshotName`, one line per row (see checkLinesMatchSnapshot).
template <typename Derived>
void checkMatrixMatchesSnapshot(const Eigen::DenseBase<Derived>& matrix, const std::string& snapshotName,
                                const Tolerance tolerance)
{
  checkLinesMatchSnapshot(matrixToLines(matrix), snapshotName, tolerance);
}

/// Every `stride`-th row of a matrix, starting with row 0, to keep snapshots of large outputs small.
inline Eigen::ArrayXXf everyNthRow(const Eigen::ArrayXXf& matrix, const Eigen::Index stride)
{
  const Eigen::Index numRows = (matrix.rows() + stride - 1) / stride;
  Eigen::ArrayXXf rows(numRows, matrix.cols());
  for (Eigen::Index row = 0; row < numRows; ++row) {
    rows.row(row) = matrix.row(row * stride);
  }
  return rows;
}

} // namespace asmc::test

#endif // ASMC_REGRESSION_SNAPSHOTS_HPP
