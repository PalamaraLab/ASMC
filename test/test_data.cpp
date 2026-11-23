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

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "test_fixtures.hpp"

using Catch::Matchers::ContainsSubstring;

namespace
{

constexpr const char* asmcExampleRoot = ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array";
constexpr const char* fastSmcExampleRoot = ASMC_DATA_DIR "/examples/fastsmc/example";
constexpr const char* asmcDecodingQuantities =
    ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";
constexpr const char* fastSmcDecodingQuantities =
    ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz";

/// FastSMC parameters for job `jobInd` of `jobs`, with or without the hashing step.
DecodingParams fastSmcParams(const bool hashing, const int jobs, const int jobInd)
{
  DecodingParams params(fastSmcExampleRoot, fastSmcDecodingQuantities, "out", /*_fastSMC=*/true);
  params.hashing = hashing;
  params.jobs = jobs;
  params.jobInd = jobInd;
  params.validateParamsFastSMC();
  return params;
}

} // namespace

TEST_CASE("ASMC loads every sample, whatever the job", "[Data]")
{
  const Data allJobs(DecodingParams(asmcExampleRoot, asmcDecodingQuantities));
  REQUIRE(allJobs.individuals.size() == allJobs.sampleSize);

  // HMM::decodeAll divides the pairs between jobs, so every job needs every sample, with the right genotypes
  const int jobs = 4;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(DecodingParams(asmcExampleRoot, asmcDecodingQuantities, "", jobs, jobInd));
    REQUIRE(oneJob.individuals.size() == allJobs.individuals.size());
    for (auto ind = 0ul; ind < allJobs.individuals.size(); ++ind) {
      CHECK(oneJob.individuals.at(ind).genotype1 == allJobs.individuals.at(ind).genotype1);
      CHECK(oneJob.individuals.at(ind).genotype2 == allJobs.individuals.at(ind).genotype2);
    }
  }
}

TEST_CASE("FastSMC without hashing loads every sample", "[Data]")
{
  const int jobs = 4;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(fastSmcParams(/*hashing=*/false, jobs, jobInd));
    CHECK(oneJob.individuals.size() == oneJob.sampleSize);
  }
}

TEST_CASE("FastSMC with hashing loads only its job's windows, and the jobs cover every sample", "[Data]")
{
  // With 4 jobs there are two windows. Jobs 1 and 4 compare a window with itself and load only that window; jobs 2
  // and 3 compare the two windows with each other and load both.
  const int jobs = 4;
  std::vector<bool> loadedByAnyJob;
  bool someJobLoadsASubset = false;
  for (int jobInd = 1; jobInd <= jobs; ++jobInd) {
    INFO("job " << jobInd << " of " << jobs);
    const Data oneJob(fastSmcParams(/*hashing=*/true, jobs, jobInd));
    loadedByAnyJob.resize(oneJob.sampleSize, false);

    auto numInJob = 0ul;
    for (auto sample = 0ul; sample < oneJob.sampleSize; ++sample) {
      if (oneJob.isSampleInJob(sample)) {
        ++numInJob;
        loadedByAnyJob.at(sample) = true;
      }
    }
    CHECK(oneJob.individuals.size() == numInJob);
    someJobLoadsASubset = someJobLoadsASubset || numInJob < oneJob.sampleSize;
  }

  CHECK(someJobLoadsASubset);
  CHECK(std::find(loadedByAnyJob.begin(), loadedByAnyJob.end(), false) == loadedByAnyJob.end());
}

TEST_CASE("blank lines in the samples and haps files are skipped", "[Data]")
{
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);
  asmc::test::writeTextFile(input.sampleFile, "ID_1 ID_2 missing\n0 0 0\n\nF1 I1 0\n   \nF2 I2 0\n\n");
  asmc::test::writeTextFile(input.hapFile, "1 rs1 100 A G 0 1 1 0\n\n"
                                           "1 rs2 200 A G 1 0 1 0\n"
                                           "1 rs3 300 A G 0 0 1 1\n\n");

  const Data data(asmc::test::makeTinyAsmcParams(input.root));
  CHECK(data.sampleSize == 2ul);
  CHECK(data.individuals.size() == 2ul);
  CHECK(data.sites == 3);
}

TEST_CASE("a haps record that cannot be read is reported", "[Data]")
{
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);
  asmc::test::writeTextFile(input.hapFile, "1 rs1 100 A G 0 1 1 0\n"
                                           "1 rs2 not_a_position A G 1 0 1 0\n"
                                           "1 rs3 300 A G 0 0 1 1\n");

  CHECK_THROWS_WITH(Data(asmc::test::makeTinyAsmcParams(input.root)), ContainsSubstring("could not read site 2 of 3"));
}

TEST_CASE("an explicit map file is used, and must exist", "[Data]")
{
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);
  const auto elsewhere = temp.file("elsewhere.map");
  std::filesystem::rename(input.mapFile, elsewhere);

  DecodingParams params = asmc::test::makeTinyAsmcParams(input.root);
  params.mapFile = elsewhere.string();
  CHECK(Data(params).geneticPositions.size() == 3ul);

  params.mapFile = temp.file("missing.map").string();
  CHECK_THROWS_WITH(Data(params), ContainsSubstring("missing.map does not exist"));
}

TEST_CASE("input files are found under every supported name", "[Data]")
{
  const asmc::test::TempDirectory temp;
  const std::string root = temp.file("input").string();

  CHECK_THROWS_AS(Data::findSamplesFile(root), std::runtime_error);
  CHECK_THROWS_AS(Data::findHapsFile(root), std::runtime_error);

  asmc::test::writeTextFile(root + ".sample", "");
  asmc::test::writeTextFile(root + ".haps", "");
  CHECK(Data::findSamplesFile(root) == root + ".sample");
  CHECK(Data::findHapsFile(root) == root + ".haps");
}

TEST_CASE("chromosome names are parsed to numbers", "[Data]")
{
  CHECK(Data::parseChromosomeNumber("1") == 1);
  CHECK(Data::parseChromosomeNumber("chr7") == 7);
  CHECK(Data::parseChromosomeNumber("CHR22") == 22);
  CHECK(Data::parseChromosomeNumber("1:12345") == 1);
  CHECK(Data::parseChromosomeNumber("chr1:5") == 1);

  // Anything else falls back to 0 rather than throwing
  CHECK(Data::parseChromosomeNumber("X") == 0);
  CHECK(Data::parseChromosomeNumber("chr") == 0);
  CHECK(Data::parseChromosomeNumber("0") == 0);
  CHECK(Data::parseChromosomeNumber("2000") == 0);
}

TEST_CASE("undistinguished counts without the CSFS stay within the samples", "[Data]")
{
  // With compress, the CSFS is not used, so there can be fewer haplotypes (here 4) than CSFS samples (here 300). The
  // counts are then undefined (-1) rather than drawn from past the end of the samples.
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);
  const DecodingParams params(input.root, ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz",
                              "", 1, 1, "array", false, false, /*_compress=*/true, false,
                              std::numeric_limits<float>::quiet_NaN());
  const Data data(params);

  for (const auto& countsAtSite : data.calculateUndistinguishedCounts(300)) {
    CHECK(countsAtSite == std::vector<int>{-1, -1, -1});
  }
}

TEST_CASE("too few samples for the CSFS is reported, with or without SNP IDs", "[Data]")
{
  // FastSMC's Data has no SNP IDs, so the message identifies the SNP by its index
  asmc::test::TempDirectory temp;
  const auto input = asmc::test::writeTinyAsmcInputRoot(temp);
  DecodingParams params(input.root, ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz", "out",
                        /*_fastSMC=*/true);
  params.hashing = false;
  const Data data(params);

  CHECK_THROWS_WITH(data.calculateUndistinguishedCounts(300), ContainsSubstring("SNP with numerical ID 0"));
}
