// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef DATA_MODULE_FILE_UTILS_HPP
#define DATA_MODULE_FILE_UTILS_HPP

#include <filesystem>
#include <string>
#include <string_view>

#include <zlib.h>

namespace asmc {

namespace fs = std::filesystem;

/**
 * Read the next line from a gzip file.
 *
 * @param gzFileHandle handle to a file opened with zlib's gzopen
 * @return a string containing the next line contained in the gzip file, without a trailing newline character
 */
std::string readNextLineFromGzip(gzFile& gzFileHandle);

/**
 * Count the number of non-empty lines in a file that may or may not be gzipped.
 *
 * @param filePath path to the file
 * @return the number of non-empty lines in the file
 */
unsigned long countLinesInFile(const fs::path& filePath);

} // namespace asmc

#endif // DATA_MODULE_FILE_UTILS_HPP
