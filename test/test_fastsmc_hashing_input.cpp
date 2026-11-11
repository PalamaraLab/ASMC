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

// Checks that FastSMC's hashing step reads its .hap input correctly. Each test runs hashing only (no HMM decoding) on
// two inputs that should give identical matches, and compares the IBD output line by line.

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <filesystem>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "DecodingParams.hpp"
#include "FastSMC.hpp"
#include "FileUtils.hpp"
#include "test_fixtures.hpp"

namespace
{

namespace fs = std::filesystem;

constexpr const char* exampleRoot = ASMC_DATA_DIR "/examples/fastsmc/example";

// Fields of a .hap line: chromosome, SNP ID, position, two alleles, then one field per haplotype.
constexpr std::size_t firstHaplotypeField = 5;

std::vector<std::string> splitFields(const std::string& line)
{
  std::vector<std::string> fields;
  std::istringstream stream(line);
  for (std::string field; stream >> field;) {
    fields.push_back(field);
  }
  return fields;
}

std::string joinFields(const std::vector<std::string>& fields)
{
  std::string line;
  for (const auto& field : fields) {
    line += (line.empty() ? "" : " ") + field;
  }
  return line;
}

/**
 * Write FastSMC input files at `root`: copies of the example .samples and .map, and a .hap.gz made by applying
 * `transformFields` to the fields of each line of the example haps.
 */
void writeExampleInput(const fs::path& root, const std::function<void(std::vector<std::string>&)>& transformFields)
{
  fs::copy_file(std::string(exampleRoot) + ".samples", root.string() + ".samples");
  fs::copy_file(std::string(exampleRoot) + ".map", root.string() + ".map");

  FileUtils::AutoGzIfstream haps;
  haps.openOrExit(std::string(exampleRoot) + ".hap.gz");
  FileUtils::AutoGzOfstream transformedHaps;
  transformedHaps.openOrExit(root.string() + ".hap.gz");
  for (std::string line; getline(haps, line);) {
    std::vector<std::string> fields = splitFields(line);
    transformFields(fields);
    transformedHaps << joinFields(fields) << '\n';
  }
  transformedHaps.close();
  haps.close();
}

/**
 * If the frequency of allele 1 at the site in `fields` is outside [minMaf, 1 - minMaf], set every haplotype to 0, so
 * that the site can't seed or extend a match. Uses the same calculation and comparison as FastSMC::run.
 *
 * @return true if the site was blanked
 */
bool blankSiteOutsideMafRange(std::vector<std::string>& fields, const float minMaf)
{
  const auto numHaplotypes = fields.size() - firstHaplotypeField;
  std::size_t numAlleleOne = 0;
  for (auto field = firstHaplotypeField; field < fields.size(); ++field) {
    numAlleleOne += fields.at(field) == "1" ? 1 : 0;
  }

  const auto maf = static_cast<float>(static_cast<double>(numAlleleOne) / static_cast<double>(numHaplotypes));
  if (maf >= minMaf && maf <= 1 - minMaf) {
    return false;
  }
  for (auto field = firstHaplotypeField; field < fields.size(); ++field) {
    fields.at(field) = "0";
  }
  return true;
}

/// Run FastSMC hashing only on the input at `root` and return the lines of its IBD output.
std::vector<std::string> runHashingOnly(const fs::path& root, const float minMaf)
{
  DecodingParams params(root.string(), ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz",
                        root.string(), /*_fastSMC=*/true);
  params.hashingOnly = true;
  params.min_maf = minMaf;
  params.useKnownSeed = true;
  params.validateParamsFastSMC();

  ASMC::FastSMC fastSMC(params);
  fastSMC.run();

  std::vector<std::string> lines;
  FileUtils::AutoGzIfstream output;
  output.openOrExit(root.string() + ".1.1.FastSMC.ibd.gz");
  for (std::string line; getline(output, line);) {
    lines.push_back(line);
  }
  return lines;
}

} // namespace

TEST_CASE("FastSMC hashing reads alleles longer than one character", "[FastSMC]")
{
  const asmc::test::TempDirectory temp;

  writeExampleInput(temp.file("original"), [](std::vector<std::string>&) {});
  writeExampleInput(temp.file("long_alleles"), [](std::vector<std::string>& fields) {
    fields.at(3) = "AT";
    fields.at(4) = "GCC";
  });

  const std::vector<std::string> original = runHashingOnly(temp.file("original"), 0.f);
  const std::vector<std::string> longAlleles = runHashingOnly(temp.file("long_alleles"), 0.f);

  REQUIRE_FALSE(original.empty());
  CHECK(longAlleles == original);
}

TEST_CASE("FastSMC --min_maf keeps hashing words aligned with sites", "[FastSMC]")
{
  const asmc::test::TempDirectory temp;
  const float minMaf = 0.05f;

  // Reference input: sites outside the MAF range are blanked by hand. Filtering the original sites with --min_maf
  // should then give exactly the same matches, at the same positions.
  std::size_t numFilteredSites = 0;
  writeExampleInput(temp.file("filtered_by_hand"), [&numFilteredSites, minMaf](std::vector<std::string>& fields) {
    numFilteredSites += blankSiteOutsideMafRange(fields, minMaf) ? 1 : 0;
  });
  writeExampleInput(temp.file("filtered_by_option"), [](std::vector<std::string>&) {});

  const std::vector<std::string> byHand = runHashingOnly(temp.file("filtered_by_hand"), 0.f);
  const std::vector<std::string> byOption = runHashingOnly(temp.file("filtered_by_option"), minMaf);

  INFO("sites filtered: " << numFilteredSites);
  REQUIRE(numFilteredSites > 0);
  REQUIRE_FALSE(byHand.empty());
  CHECK(byOption == byHand);
}
