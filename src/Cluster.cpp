// Per-site UPGMA over ASMC posterior-mean TMRCAs: see Cluster.hpp.
//
// A native version of the inner loop of ASMC-clust, ported from arg-needle's Python upgma_chunk
// (arg_needle/inference.py; Zhang et al., Nat Genet 2023).

#include "Cluster.hpp"

#include "ASMC.hpp"
#include "DecodePairsReturnStruct.hpp"

#include "../3rd_party/fastcluster/fastcluster_facade.h"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ASMC
{
namespace
{

/// The tree UPGMA builds at one anchor site, in scipy linkage order: merge k joins children[k] at heights[k]. A child
/// ID below numHaps is a sample; ID numHaps + j is the node made by merge j. Heights strictly increase.
struct SiteTree {
  std::vector<double> heights;
  std::vector<std::pair<int, int>> children;
};

/// Which sites of a chunk run UPGMA (the anchors), and whose tree every site uses. Indexed by site within the chunk.
struct AnchorLayout {
  std::vector<int> anchorAtSite;  // index of the anchor at this site, or -1 if the site is not an anchor
  std::vector<int> nearestAnchor; // index of the anchor whose tree this site uses
  int numAnchors = 0;
};

enum class TreeStatus { ok, linkageFailed, zeroDistance };

/// Number of the calling OpenMP thread, or 0 without OpenMP.
int currentThread()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

/**
 * Boundaries of the genomic interval each site's tree covers: chunkLength + 1 values, site s covering [boundary s,
 * boundary s + 1). Interior boundaries are midway between neighbouring sites' physical positions, less argOffset. The
 * outer two are 0 and argLength, as in upgma_chunk, or the first and last sites' positions if argLength is 0.
 */
std::vector<double> makeSiteBoundaries(const std::vector<int>& physicalPositions, const int chunkFirst,
                                       const int chunkLength, const int argOffset, const double argLength)
{
  std::vector<double> boundaries(static_cast<std::size_t>(chunkLength) + 1);
  if (argLength > 0.0) {
    boundaries.front() = 0.0;
    boundaries.back() = argLength;
  } else {
    boundaries.front() = static_cast<double>(physicalPositions[chunkFirst] - argOffset);
    boundaries.back() = static_cast<double>(physicalPositions[chunkFirst + chunkLength - 1] - argOffset);
  }
  for (int site = 1; site < chunkLength; ++site) {
    const int position = chunkFirst + site;
    boundaries[site] = 0.5 * (static_cast<double>(physicalPositions[position - 1] - argOffset) +
                              static_cast<double>(physicalPositions[position] - argOffset));
  }
  return boundaries;
}

/// Anchors at every stride-th site and at the last site. Every site uses the nearest anchor's tree, ties going left.
AnchorLayout makeAnchorLayout(const int chunkLength, const int stride)
{
  AnchorLayout layout;
  layout.anchorAtSite.assign(chunkLength, -1);
  std::vector<int> anchorSites;
  for (int site = 0; site < chunkLength; ++site) {
    if (site % stride == 0 || site == chunkLength - 1) {
      layout.anchorAtSite[site] = static_cast<int>(anchorSites.size());
      anchorSites.push_back(site);
    }
  }
  layout.numAnchors = static_cast<int>(anchorSites.size());

  layout.nearestAnchor.assign(chunkLength, 0);
  for (int site = 0; site < chunkLength; ++site) {
    // The first anchor at or after the site, and the one before it
    const auto after =
        static_cast<int>(std::lower_bound(anchorSites.begin(), anchorSites.end(), site) - anchorSites.begin());
    if (after == 0) {
      layout.nearestAnchor[site] = 0;
    } else if (after == layout.numAnchors) {
      layout.nearestAnchor[site] = layout.numAnchors - 1;
    } else {
      const bool beforeIsNearer = site - anchorSites[after - 1] <= anchorSites[after] - site;
      layout.nearestAnchor[site] = beforeIsNearer ? after - 1 : after;
    }
  }
  return layout;
}

/// All pairs of numHaps haplotypes in scipy's condensed order, (0, 1), (0, 2), ..., (0, n - 1), (1, 2), ..., so that
/// row k of decodePairs's results is entry k of a condensed distance matrix. Each pair is (larger, smaller) index.
std::pair<std::vector<unsigned long>, std::vector<unsigned long>> makeCondensedPairs(const int numHaps)
{
  std::vector<unsigned long> hapsA;
  std::vector<unsigned long> hapsB;
  const auto numPairs = static_cast<std::size_t>(numHaps) * (numHaps - 1) / 2;
  hapsA.reserve(numPairs);
  hapsB.reserve(numPairs);
  for (int smaller = 0; smaller < numHaps - 1; ++smaller) {
    for (int larger = smaller + 1; larger < numHaps; ++larger) {
      hapsA.push_back(static_cast<unsigned long>(larger));
      hapsB.push_back(static_cast<unsigned long>(smaller));
    }
  }
  return {std::move(hapsA), std::move(hapsB)};
}

/// Decode every pair over sites [from, to), keeping posterior means only. n_threads 0 means decodePairs.
void decodeBlock(ASMC& asmc, const std::vector<unsigned long>& hapsA, const std::vector<unsigned long>& hapsB,
                 const int from, const int to, const float cmBurnIn, const int numThreads)
{
  if (numThreads > 0) {
    asmc.decodePairsParallel(hapsA, hapsB, static_cast<unsigned>(from), static_cast<unsigned>(to), cmBurnIn,
                             numThreads);
  } else {
    asmc.decodePairs(hapsA, hapsB, static_cast<unsigned>(from), static_cast<unsigned>(to), cmBurnIn);
  }
}

/**
 * Build the UPGMA tree of one site from the pairs' posterior means in one column of decodePairs's results.
 *
 * @param distances scratch for the condensed distance matrix, numHaps (numHaps - 1) / 2 values
 * @param linkage scratch for fastcluster's linkage matrix, 4 (numHaps - 1) values
 * @param tree output, with numHaps - 1 merges
 */
TreeStatus buildSiteTree(const DecodePairsReturnStruct& results, const Eigen::Index column, const int numHaps,
                         std::vector<double>& distances, std::vector<double>& linkage, SiteTree& tree)
{
  for (std::size_t pair = 0; pair < distances.size(); ++pair) {
    distances[pair] = static_cast<double>(results.perPairPosteriorMeans(static_cast<Eigen::Index>(pair), column));
  }
  if (asmc_fastcluster_average_linkage(numHaps, distances.data(), linkage.data()) != 0) {
    return TreeStatus::linkageFailed;
  }

  // Tied merge heights would give a parent the same time as its child, which tskit rejects. As in upgma_chunk, each
  // height is nudged above the one before.
  double lastHeight = 0.0;
  for (int merge = 0; merge < numHaps - 1; ++merge) {
    const double* row = linkage.data() + static_cast<std::size_t>(merge) * 4; // [child a, child b, height, size]
    double height = row[2];
    if (height == 0.0) {
      return TreeStatus::zeroDistance;
    }
    if (height <= lastHeight) {
      height = lastHeight * 1.00001;
    }
    lastHeight = height;
    tree.heights[merge] = height;
    tree.children[merge] = std::make_pair(static_cast<int>(row[0]), static_cast<int>(row[1]));
  }
  return TreeStatus::ok;
}

/// Throw if any anchor's tree could not be built
void checkTrees(const std::vector<TreeStatus>& statuses)
{
  for (const TreeStatus status : statuses) {
    if (status == TreeStatus::linkageFailed) {
      throw std::runtime_error("clusterChunk: fastcluster failed at a site");
    }
    if (status == TreeStatus::zeroDistance) {
      throw std::runtime_error("clusterChunk: zero merge distance at a site");
    }
  }
}

/// Write the node times and edges of one site from its anchor's tree. The internal nodes of site s are numbered
/// numHaps + s (numHaps - 1) + merge; edges cover the site's interval.
void emitSite(const SiteTree& tree, const int numHaps, const int site, ClusterChunkResult& out)
{
  const double left = out.site_boundaries[site];
  const double right = out.site_boundaries[site + 1];
  const int firstNode = numHaps + site * (numHaps - 1);
  const std::size_t nodeBase = static_cast<std::size_t>(site) * (numHaps - 1);
  const std::size_t edgeBase = nodeBase * 2;

  for (int merge = 0; merge < numHaps - 1; ++merge) {
    out.node_time[nodeBase + merge] = tree.heights[merge];
    const std::array<int, 2> children = {tree.children[merge].first, tree.children[merge].second};
    for (std::size_t c = 0; c < children.size(); ++c) {
      const int childNode = children[c] < numHaps ? children[c] : firstNode + (children[c] - numHaps);
      const std::size_t edge = edgeBase + static_cast<std::size_t>(merge) * 2 + c;
      out.edge_left[edge] = left;
      out.edge_right[edge] = right;
      out.edge_parent[edge] = firstNode + merge;
      out.edge_child[edge] = childNode;
    }
  }
}

} // namespace

ClusterChunkResult clusterChunk(ASMC& asmc, unsigned from, unsigned to, unsigned site_block, int n_threads,
                                float cm_burn_in, int arg_offset, double arg_length, unsigned anchor_stride_k)
{
  if (asmc.getHaploidSampleSize() > static_cast<unsigned long>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("clusterChunk: numHaps exceeds int range");
  }
  const auto numHaps = static_cast<int>(asmc.getHaploidSampleSize());
  if (numHaps < 2) {
    throw std::invalid_argument("clusterChunk: need at least 2 haploids");
  }
  const int numSites = asmc.getNumSites();
  if (to == 0u) {
    to = static_cast<unsigned>(numSites);
  }
  if (from >= to || static_cast<int>(to) > numSites) {
    throw std::invalid_argument("clusterChunk: invalid [from, to)");
  }

  const auto chunkFirst = static_cast<int>(from);
  const auto chunkLength = static_cast<int>(to - from);
  const int blockLength =
      site_block == 0u ? chunkLength : static_cast<int>(std::min<unsigned>(site_block, chunkLength));
  const AnchorLayout layout = makeAnchorLayout(chunkLength, static_cast<int>(std::max(1u, anchor_stride_k)));
  const auto [hapsA, hapsB] = makeCondensedPairs(numHaps);

  ClusterChunkResult out;
  out.num_haps = numHaps;
  out.num_sites = chunkLength;
  out.site_boundaries =
      makeSiteBoundaries(asmc.getPhysicalPositions(), chunkFirst, chunkLength, arg_offset, arg_length);
  const std::size_t numNodes = static_cast<std::size_t>(chunkLength) * (numHaps - 1);
  out.node_time.assign(numNodes, 0.0);
  out.edge_left.assign(2 * numNodes, 0.0);
  out.edge_right.assign(2 * numNodes, 0.0);
  out.edge_parent.assign(2 * numNodes, 0);
  out.edge_child.assign(2 * numNodes, 0);

  // Only posterior means are needed
  asmc.setStorePerPairPosteriorMean(true);
  asmc.setStorePerPairMap(false);
  asmc.setStorePerPairPosterior(false);
  asmc.setStoreSumOfPosterior(false);

  // Pass 1, block by block: decode all pairs over the block, then build the tree at each anchor in it. Each anchor's
  // tree and status are written by one thread only, and each thread has its own scratch buffers.
  const int numThreads = std::max(1, n_threads);
  std::vector<SiteTree> anchorTrees(
      layout.numAnchors, SiteTree{std::vector<double>(numHaps - 1), std::vector<std::pair<int, int>>(numHaps - 1)});
  std::vector<TreeStatus> anchorStatus(layout.numAnchors, TreeStatus::ok);
  std::vector<std::vector<double>> distancesPerThread(numThreads, std::vector<double>(hapsA.size()));
  std::vector<std::vector<double>> linkagePerThread(numThreads, std::vector<double>(4 * (numHaps - 1)));

  for (int blockStart = 0; blockStart < chunkLength; blockStart += blockLength) {
    const int blockSites = std::min(blockLength, chunkLength - blockStart);
    decodeBlock(asmc, hapsA, hapsB, chunkFirst + blockStart, chunkFirst + blockStart + blockSites, cm_burn_in,
                n_threads);
    const DecodePairsReturnStruct& results = asmc.getRefOfResults();

#pragma omp parallel for num_threads(numThreads) schedule(static)
    for (int site = 0; site < blockSites; ++site) {
      const int anchor = layout.anchorAtSite[blockStart + site];
      if (anchor >= 0) {
        const int thread = currentThread();
        anchorStatus[anchor] = buildSiteTree(results, site, numHaps, distancesPerThread[thread],
                                             linkagePerThread[thread], anchorTrees[anchor]);
      }
    }
    checkTrees(anchorStatus);
  }

  // Pass 2: write every site from its nearest anchor's tree, which may have been in a later block
#pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int site = 0; site < chunkLength; ++site) {
    emitSite(anchorTrees[layout.nearestAnchor[site]], numHaps, site, out);
  }

  return out;
}

} // namespace ASMC
