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

#ifndef ASMC_TEST_FIXTURES_HPP
#define ASMC_TEST_FIXTURES_HPP

#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"

#include <Eigen/Core>

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace asmc::test {

namespace fs = std::filesystem;

using RowMajorPosterior = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

class TempDirectory {
public:
  explicit TempDirectory(std::string_view prefix = "asmc_test") : m_path(makeUniquePath(prefix))
  {
    fs::create_directories(m_path);
  }

  TempDirectory(const TempDirectory&) = delete;
  TempDirectory& operator=(const TempDirectory&) = delete;

  ~TempDirectory()
  {
    std::error_code error;
    fs::remove_all(m_path, error);
  }

  [[nodiscard]] const fs::path& path() const
  {
    return m_path;
  }

  [[nodiscard]] fs::path file(std::string_view name) const
  {
    return m_path / std::string(name);
  }

private:
  static fs::path makeUniquePath(std::string_view prefix)
  {
    static std::atomic<unsigned long> counter{0ul};
    const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    return fs::temp_directory_path() / (std::string(prefix) + "_" + std::to_string(now) + "_" +
                                        std::to_string(counter.fetch_add(1ul)));
  }

  fs::path m_path;
};

inline void writeTextFile(const fs::path& path, std::string_view contents)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not create test fixture file: " + path.string());
  }
  out << contents;
}

struct TinyAsmcInput {
  std::string root;
  fs::path sampleFile;
  fs::path hapFile;
  fs::path mapFile;
};

inline TinyAsmcInput writeTinyAsmcInputRoot(TempDirectory& temp, std::string_view stem = "tiny")
{
  const fs::path rootPath = temp.file(stem);
  TinyAsmcInput input{rootPath.string(), rootPath.string() + ".samples", rootPath.string() + ".hap",
                      rootPath.string() + ".map"};

  writeTextFile(input.sampleFile, "ID_1 ID_2 missing\n0 0 0\nF1 I1 0\nF2 I2 0\n");
  writeTextFile(input.hapFile,
                "1 rs1 100 A G 0 1 1 0\n"
                "1 rs2 200 A G 1 0 1 0\n"
                "1 rs3 300 A G 0 0 1 1\n");
  writeTextFile(input.mapFile,
                "1\trs1\t0.0\t100\n"
                "1\trs2\t0.1\t200\n"
                "1\trs3\t0.3\t300\n");

  return input;
}

inline DecodingParams makeTinyAsmcParams(const std::string& inputRoot,
                                         std::string decodingQuantitiesFile =
                                             ASMC_DATA_DIR "/decoding_quantities/10-20-2000_CEU.decodingQuantities.gz")
{
  DecodingParams params(inputRoot, std::move(decodingQuantitiesFile));
  params.useKnownSeed = true;
  return params;
}

inline Data makeTinyAsmcData(TempDirectory& temp, std::string_view stem = "tiny")
{
  const TinyAsmcInput input = writeTinyAsmcInputRoot(temp, stem);
  return Data(makeTinyAsmcParams(input.root));
}

inline PairObservations makePairObservations(std::vector<bool> obsBits, std::vector<bool> homMinorBits,
                                             unsigned int iInd = 0u, unsigned int jInd = 1u,
                                             uint_least8_t iHap = 1u, uint_least8_t jHap = 1u)
{
  if (obsBits.size() != homMinorBits.size()) {
    throw std::runtime_error("Pair observation fixture vectors must have matching lengths");
  }

  PairObservations obs;
  obs.iHap = iHap;
  obs.jHap = jHap;
  obs.iInd = iInd;
  obs.jInd = jInd;
  obs.obsBits = std::move(obsBits);
  obs.homMinorBits = std::move(homMinorBits);
  return obs;
}

inline RowMajorPosterior makePosteriorMatrix(const std::vector<std::vector<float>>& valuesByState)
{
  if (valuesByState.empty()) {
    return RowMajorPosterior{};
  }

  const auto numSites = valuesByState.front().size();
  RowMajorPosterior posterior(static_cast<Eigen::Index>(valuesByState.size()), static_cast<Eigen::Index>(numSites));
  for (auto state = 0ul; state < valuesByState.size(); ++state) {
    if (valuesByState.at(state).size() != numSites) {
      throw std::runtime_error("Posterior fixture rows must all have the same number of sites");
    }
    for (auto site = 0ul; site < numSites; ++site) {
      posterior(static_cast<Eigen::Index>(state), static_cast<Eigen::Index>(site)) = valuesByState.at(state).at(site);
    }
  }
  return posterior;
}

inline RowMajorPosterior makeColumnNormalisedPosterior(unsigned int states, unsigned int sites)
{
  RowMajorPosterior posterior(states, sites);
  for (auto site = 0u; site < sites; ++site) {
    float sum = 0.f;
    for (auto state = 0u; state < states; ++state) {
      const float value = static_cast<float>((site + 1u) * (state + 1u));
      posterior(state, site) = value;
      sum += value;
    }
    posterior.col(site) /= sum;
  }
  return posterior;
}

inline Eigen::ArrayXf makeBatchedPosteriorBuffer(const std::vector<RowMajorPosterior>& posteriors,
                                                 int paddedBatchSize = -1)
{
  if (posteriors.empty()) {
    return Eigen::ArrayXf{};
  }

  const Eigen::Index states = posteriors.front().rows();
  const Eigen::Index sites = posteriors.front().cols();
  const int batchSize = paddedBatchSize < 0 ? static_cast<int>(posteriors.size()) : paddedBatchSize;
  if (batchSize < static_cast<int>(posteriors.size())) {
    throw std::runtime_error("Padded batch size cannot be smaller than the number of posterior matrices");
  }

  Eigen::ArrayXf buffer(sites * states * batchSize);
  buffer.setZero();

  for (auto batch = 0ul; batch < posteriors.size(); ++batch) {
    if (posteriors.at(batch).rows() != states || posteriors.at(batch).cols() != sites) {
      throw std::runtime_error("All batched posterior fixtures must have the same dimensions");
    }
    for (Eigen::Index site = 0; site < sites; ++site) {
      for (Eigen::Index state = 0; state < states; ++state) {
        buffer((site * states + state) * batchSize + static_cast<Eigen::Index>(batch)) =
            posteriors.at(batch)(state, site);
      }
    }
  }

  return buffer;
}

inline DecodePairsReturnStruct makeDecodePairsReturnStruct(const std::vector<unsigned long>& hapIndicesA,
                                                           const std::vector<unsigned long>& hapIndicesB,
                                                           long int numSites, long int numStates,
                                                           bool storeFullPosteriors = false,
                                                           bool storeSumOfPosteriors = false,
                                                           bool storePerPairPosteriors = false,
                                                           bool storePerPairMaps = false)
{
  DecodePairsReturnStruct result;
  result.initialise(hapIndicesA, hapIndicesB, numSites, numStates, storeFullPosteriors, storeSumOfPosteriors,
                    storePerPairPosteriors, storePerPairMaps);
  return result;
}

} // namespace asmc::test

#endif
