// ASMC clusterChunk — per-site UPGMA over ASMC posterior-mean TMRCAs.
//
// Native re-implementation of arg-needle's "ASMC-clust" inner loop
// (Zhang et al., Nat Genet 2023). One ASMC decode_pairs[_parallel]
// call per site-block batches all n*(n-1)/2 pairs in scipy condensed
// order; an OMP-parallel UPGMA at every (anchor) site emits a flat
// edges/nodes table suitable for tskit / arg-needle-lib ingest.

#ifndef ASMC_CLUSTER_HPP
#define ASMC_CLUSTER_HPP

#include "ASMC.hpp"

#include <cstdint>
#include <vector>

namespace ASMC {

/// Result of clusterChunk — flat arrays suitable for direct conversion
/// to a tskit TableCollection or arg-needle-lib ARG.
///
/// Node IDs:
///   * [0, num_haps)        original sample nodes; nodeTime entry is 0.0
///   * [num_haps, num_haps + numInternalNodes) internal merge nodes
///
/// Edges (one row of edge_left/right/parent/child per entry) carry the
/// site-interval for which the edge is active. Each per-site UPGMA
/// emits 2(num_haps - 1) edges and (num_haps - 1) internal nodes.
struct ClusterChunkResult {
  /// num_haps × per-site count of internal nodes; equals num_internal_nodes
  /// per site for a strict bifurcating tree (numHaps - 1).
  int num_haps = 0;
  int num_sites = 0;
  /// Site bin boundaries in physical coordinates, length num_sites + 1.
  /// boundary[s] is the left edge for site s, boundary[s+1] is the right.
  std::vector<double> site_boundaries;
  /// All internal-node coalescence times in chunk order. Length =
  /// num_sites × (num_haps - 1). The first (num_haps - 1) are for site 0.
  std::vector<double> node_time;
  /// Edge tables, length num_sites × 2 × (num_haps - 1). Edges for site s
  /// have left == site_boundaries[s], right == site_boundaries[s+1].
  std::vector<double> edge_left;
  std::vector<double> edge_right;
  std::vector<int>   edge_parent;
  std::vector<int>   edge_child;

  /// Total number of internal nodes added across all sites
  /// (equals num_sites × (num_haps - 1)).
  int num_internal_nodes() const {
    return num_sites * (num_haps - 1);
  }
};

/// Run per-site UPGMA on ASMC posterior-mean TMRCAs across [from, to).
///
/// Behaviour matches arg_needle.inference.upgma_chunk for the same input,
/// up to the deterministic float ordering inside fastcluster — outputs
/// have identical edges and identical heights (within float rounding).
///
/// from           first site (0-based, inclusive)
/// to             last site (0-based, exclusive). 0 → use numSites.
/// site_block     UPGMA is run on blocks of this many sites at a time
///                so the resident pair-distance buffer is bounded by
///                numPairs × site_block × 4 B. 0 → one big block.
/// n_threads      threads for both the underlying decodePairsParallel and
///                the per-site UPGMA loop. 0 → decodePairs and one thread.
/// cm_burn_in     padding (cM) for the HMM at chunk boundaries (same
///                semantic as decode_pairs).
/// arg_offset     subtract from physical positions when computing
///                interior site_boundaries; matches upgma_chunk's
///                "arg_start" convention.
/// arg_length     total ARG length (used to anchor the trailing boundary
///                to the chunk-open / chunk-close edges; matches
///                upgma_chunk's `boundaries[0] = 0`, `boundaries[-1] =
///                arg_length`). If 0, fall back to physical-position
///                anchoring (legacy / debugging only).
/// anchor_stride_k Run UPGMA only at every K-th site (anchors). At
///                non-anchor sites, copy the nearest anchor's topology
///                and merge heights into that site's emit. K==1 (default)
///                = full per-site UPGMA, bit-identical to arg-needle's
///                Python `upgma_chunk` for the same input. K>1 saves
///                O((1 - 1/K) * num_sites * n²) of UPGMA work; topology
///                penalty is bounded by the truth-tree drift across K
///                sites in cM. K is in *sites*, not cM — for safety pick
///                K so K * (cM / site) ≲ 0.01 cM (empirical sweet spot at
///                array-density data; tighter for sequencing density).
///
/// Sets asmc's per-pair outputs to posterior means only, and leaves them so.
/// Throws std::invalid_argument for fewer than 2 haplotypes or an invalid
/// [from, to), and std::runtime_error if a site's tree cannot be built.
ClusterChunkResult clusterChunk(class ASMC& asmc,
                                unsigned from = 0u,
                                unsigned to = 0u,
                                unsigned site_block = 0u,
                                int n_threads = 0,
                                float cm_burn_in = 0.5f,
                                int arg_offset = 0,
                                double arg_length = 0.0,
                                unsigned anchor_stride_k = 1u);

} // namespace ASMC

#endif // ASMC_CLUSTER_HPP
