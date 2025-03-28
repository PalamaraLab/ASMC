// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <array>
#include <string>

#include <zlib.h>

namespace asmc {

std::string readNextLineFromGzip(gzFile& gzFileHandle) {

  std::array<char, 512> buffer = {};
  std::string line;

  do {
    auto successful_read = gzgets(gzFileHandle, buffer.data(), static_cast<int>(buffer.size()));

    if (successful_read != Z_NULL) {
      line += buffer.data();
    }

  } while (!line.empty() && line.back() != '\n' && !gzeof(gzFileHandle));

  return stripBack(line);
}

unsigned long countLinesInFile(const fs::path& filePath) {
  auto gzFile = gzopen(filePath.string().c_str(), "r");

  unsigned long numLines = 0ul;
  while (!gzeof(gzFile)) {
    std::string line = readNextLineFromGzip(gzFile);
    if (!line.empty()) {
      numLines++;
    }
  }

  gzclose(gzFile);

  return numLines;
}

} // namespace asmc
