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

// Unit tests for ASMC::clusterChunk (the in-process ASMC-clust kernel). They use the first 20 haplotypes of the n300
// array example, held in memory, so that each call decodes 190 pairs rather than 44,850.

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "ASMC.hpp"
#include "Cluster.hpp"
#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"

namespace
{

constexpr int numHaps = 20;
const char* const decodingQuantFile = ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz";

/// An ASMC object for the first numHaps haplotypes of the n300 array example, without the CSFS (which needs at least
/// 300 haplotypes).
std::unique_ptr<ASMC::ASMC> makeSmallAsmc()
{
  const DecodingParams fileParams(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array", decodingQuantFile);
  const Data full(fileParams);

  // Laid out [site * numHaps + haplotype]
  std::vector<int8_t> haplotypes(static_cast<std::size_t>(full.sites) * numHaps);
  for (int site = 0; site < full.sites; ++site) {
    for (int ind = 0; ind < numHaps / 2; ++ind) {
      haplotypes[site * numHaps + 2 * ind] = full.individuals[ind].genotype1[site] ? 1 : 0;
      haplotypes[site * numHaps + 2 * ind + 1] = full.individuals[ind].genotype2[site] ? 1 : 0;
    }
  }

  DecodingParams params("", decodingQuantFile, "", 1, 1, "array", /*_decodingSequence=*/false, /*_usingCSFS=*/false,
                        /*_compress=*/true);
  Data data(full.sites, numHaps / 2, full.physicalPositions, full.geneticPositions, haplotypes.data(), params.usingCSFS,
            params.foldData);
  return std::make_unique<ASMC::ASMC>(params, std::move(data));
}

ASMC::ClusterChunkResult runChunk(ASMC::ASMC& asmc, const unsigned from, const unsigned to,
                                  const unsigned anchorStride = 1u, const unsigned siteBlock = 0u,
                                  const int numThreads = 0)
{
  return ASMC::clusterChunk(asmc, from, to, siteBlock, numThreads, /*cm_burn_in=*/0.5f, /*arg_offset=*/0,
                            /*arg_length=*/0.0, anchorStride);
}

/// One merge of a tree: its height and the sorted samples below it.
using Merge = std::pair<double, std::vector<int>>;

/// The merges of one site's tree in clusterChunk's output, in output order.
std::vector<Merge> mergesAtSite(const ASMC::ClusterChunkResult& result, const int site)
{
  const int n = result.num_haps;
  const int firstNode = n + site * (n - 1);
  std::vector<Merge> merges;
  for (int merge = 0; merge < n - 1; ++merge) {
    std::vector<int> samples;
    for (int edge = 0; edge < 2; ++edge) {
      const int child = result.edge_child[(static_cast<std::size_t>(site) * (n - 1) + merge) * 2 + edge];
      if (child < n) {
        samples.push_back(child);
      } else {
        const std::vector<int>& below = merges.at(child - firstNode).second;
        samples.insert(samples.end(), below.begin(), below.end());
      }
    }
    std::sort(samples.begin(), samples.end());
    merges.emplace_back(result.node_time[static_cast<std::size_t>(site) * (n - 1) + merge], std::move(samples));
  }
  return merges;
}

/// Reference UPGMA, written for clarity: repeatedly merge the two clusters with the smallest mean distance between
/// their samples. `distance(i, j)` is the distance between samples i and j.
template <typename Distance> std::vector<Merge> naiveUpgma(const int n, Distance distance)
{
  std::vector<std::vector<int>> clusters;
  for (int sample = 0; sample < n; ++sample) {
    clusters.push_back({sample});
  }
  auto meanDistance = [&](const std::vector<int>& a, const std::vector<int>& b) {
    double sum = 0.0;
    for (const int i : a) {
      for (const int j : b) {
        sum += distance(i, j);
      }
    }
    return sum / static_cast<double>(a.size() * b.size());
  };

  std::vector<Merge> merges;
  while (clusters.size() > 1) {
    std::tuple<double, std::size_t, std::size_t> closest{INFINITY, 0, 0};
    for (std::size_t a = 0; a < clusters.size(); ++a) {
      for (std::size_t b = a + 1; b < clusters.size(); ++b) {
        closest = std::min(closest, std::make_tuple(meanDistance(clusters[a], clusters[b]), a, b));
      }
    }
    const auto [height, a, b] = closest;
    std::vector<int> merged = clusters[a];
    merged.insert(merged.end(), clusters[b].begin(), clusters[b].end());
    std::sort(merged.begin(), merged.end());
    clusters.erase(clusters.begin() + static_cast<long>(b));
    clusters[a] = merged;
    merges.emplace_back(height, std::move(merged));
  }
  return merges;
}

} // namespace

TEST_CASE("clusterChunk: layout of the output", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto r = runChunk(*asmc, 100u, 108u);

  REQUIRE(r.num_haps == numHaps);
  REQUIRE(r.num_sites == 8);
  REQUIRE(r.num_internal_nodes() == 8 * (numHaps - 1));
  // Per site, (n - 1) internal nodes and 2 (n - 1) edges
  const std::size_t nodes = static_cast<std::size_t>(r.num_sites) * (numHaps - 1);
  REQUIRE(r.node_time.size() == nodes);
  REQUIRE(r.edge_left.size() == 2 * nodes);
  REQUIRE(r.edge_right.size() == 2 * nodes);
  REQUIRE(r.edge_parent.size() == 2 * nodes);
  REQUIRE(r.edge_child.size() == 2 * nodes);
  REQUIRE(r.site_boundaries.size() == static_cast<std::size_t>(r.num_sites + 1));
}

TEST_CASE("clusterChunk: heights strictly increase within each site", "[Cluster]")
{
  // clusterChunk nudges tied heights apart, because tskit needs a parent to be older than its child
  auto asmc = makeSmallAsmc();
  const auto r = runChunk(*asmc, 100u, 108u);
  for (int site = 0; site < r.num_sites; ++site) {
    for (int merge = 1; merge < numHaps - 1; ++merge) {
      const std::size_t node = static_cast<std::size_t>(site) * (numHaps - 1) + merge;
      CHECK(r.node_time[node] > r.node_time[node - 1]);
    }
  }
}

TEST_CASE("clusterChunk: edges join each site's own nodes over its interval", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto r = runChunk(*asmc, 100u, 108u);
  for (int site = 0; site < r.num_sites; ++site) {
    INFO("site " << site);
    const int firstNode = numHaps + site * (numHaps - 1);
    for (int edge = 0; edge < 2 * (numHaps - 1); ++edge) {
      const std::size_t index = static_cast<std::size_t>(site) * 2 * (numHaps - 1) + edge;
      const int parent = r.edge_parent[index];
      const int child = r.edge_child[index];
      CHECK(parent == firstNode + edge / 2);
      CHECK((child < numHaps || (child >= firstNode && child < parent)));
      CHECK(r.edge_left[index] == r.site_boundaries[site]);
      CHECK(r.edge_right[index] == r.site_boundaries[site + 1]);
    }
  }
}

TEST_CASE("clusterChunk: site boundaries increase", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto r = runChunk(*asmc, 100u, 108u);
  for (int site = 0; site < r.num_sites; ++site) {
    CHECK(r.site_boundaries[site] < r.site_boundaries[site + 1]);
  }
}

TEST_CASE("clusterChunk: each site's tree is the UPGMA tree of the pairs' posterior means", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const unsigned from = 100u;
  const unsigned to = 104u;
  const auto r = runChunk(*asmc, from, to);

  // The pairwise distances clusterChunk clusters: posterior means over the same window and burn-in. Pair (i, j)
  // with i < j is row pairRow[i][j].
  std::vector<unsigned long> hapsA;
  std::vector<unsigned long> hapsB;
  std::vector<std::vector<Eigen::Index>> pairRow(numHaps, std::vector<Eigen::Index>(numHaps));
  for (int i = 0; i < numHaps; ++i) {
    for (int j = i + 1; j < numHaps; ++j) {
      pairRow[i][j] = static_cast<Eigen::Index>(hapsA.size());
      hapsA.push_back(static_cast<unsigned long>(i));
      hapsB.push_back(static_cast<unsigned long>(j));
    }
  }
  asmc->decodePairs(hapsA, hapsB, from, to, 0.5f);
  const DecodePairsReturnStruct means = asmc->getCopyOfResults();

  for (int site = 0; site < r.num_sites; ++site) {
    INFO("site " << site);
    auto distance = [&](const int i, const int j) {
      return static_cast<double>(means.perPairPosteriorMeans(pairRow[std::min(i, j)][std::max(i, j)], site));
    };
    const std::vector<Merge> expected = naiveUpgma(numHaps, distance);
    const std::vector<Merge> actual = mergesAtSite(r, site);

    REQUIRE(actual.size() == expected.size());
    std::set<std::vector<int>> expectedClusters;
    std::set<std::vector<int>> actualClusters;
    for (std::size_t merge = 0; merge < expected.size(); ++merge) {
      CHECK(std::abs(actual[merge].first - expected[merge].first) <= 1e-9 * expected[merge].first);
      expectedClusters.insert(expected[merge].second);
      actualClusters.insert(actual[merge].second);
    }
    CHECK(actualClusters == expectedClusters);
  }
}

TEST_CASE("clusterChunk: sites between anchors use the nearest anchor's tree", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto everySite = runChunk(*asmc, 100u, 110u);
  const auto strided = runChunk(*asmc, 100u, 110u, /*anchorStride=*/4u);
  REQUIRE(strided.num_sites == everySite.num_sites);

  // Anchors are sites 0, 4, 8 and the last site, 9. Ties go to the left anchor: site 2 uses 0 and site 6 uses 4.
  const std::vector<int> expectedAnchor = {0, 0, 0, 4, 4, 4, 4, 8, 8, 9};
  for (int site = 0; site < strided.num_sites; ++site) {
    INFO("site " << site);
    CHECK(mergesAtSite(strided, site) == mergesAtSite(everySite, expectedAnchor[site]));
  }
}

TEST_CASE("clusterChunk: the result does not depend on the number of threads", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto oneThread = runChunk(*asmc, 100u, 110u, /*anchorStride=*/3u, /*siteBlock=*/4u, /*numThreads=*/0);
  const auto threeThreads = runChunk(*asmc, 100u, 110u, /*anchorStride=*/3u, /*siteBlock=*/4u, /*numThreads=*/3);
  CHECK(threeThreads.node_time == oneThread.node_time);
  CHECK(threeThreads.edge_parent == oneThread.edge_parent);
  CHECK(threeThreads.edge_child == oneThread.edge_child);
  CHECK(threeThreads.edge_left == oneThread.edge_left);
  CHECK(threeThreads.edge_right == oneThread.edge_right);
}

TEST_CASE("clusterChunk: site blocks give the same shape of output", "[Cluster]")
{
  // Each block is decoded separately, with its own burn-in, so the trees can differ slightly from one block
  auto asmc = makeSmallAsmc();
  const auto oneBlock = runChunk(*asmc, 100u, 108u, 1u, /*siteBlock=*/0u);
  const auto smallBlocks = runChunk(*asmc, 100u, 108u, 1u, /*siteBlock=*/3u);
  REQUIRE(smallBlocks.num_sites == oneBlock.num_sites);
  REQUIRE(smallBlocks.node_time.size() == oneBlock.node_time.size());
  CHECK(smallBlocks.site_boundaries == oneBlock.site_boundaries);
}

TEST_CASE("clusterChunk: an invalid range throws", "[Cluster]")
{
  auto asmc = makeSmallAsmc();
  const auto numSites = static_cast<unsigned>(asmc->getNumSites());
  CHECK_THROWS_AS(runChunk(*asmc, 10u, 10u), std::invalid_argument);
  CHECK_THROWS_AS(runChunk(*asmc, 0u, numSites + 1u), std::invalid_argument);
}
