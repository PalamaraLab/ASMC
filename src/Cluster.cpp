// See Cluster.hpp.

#include "Cluster.hpp"

#include "ASMC.hpp"
#include "DecodePairsReturnStruct.hpp"

#include "../3rd_party/fastcluster/fastcluster_facade.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ASMC {

ClusterChunkResult clusterChunk(ASMC& asmc,
                                unsigned from, unsigned to,
                                unsigned site_block,
                                int n_threads,
                                float cm_burn_in,
                                int arg_offset,
                                double arg_length,
                                unsigned anchor_stride_k) {
  const auto num_haps_ul = asmc.getHaploidSampleSize();
  if (num_haps_ul > static_cast<unsigned long>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("clusterChunk: numHaps exceeds int range");
  }
  const int n = static_cast<int>(num_haps_ul);
  const int total_sites = asmc.getNumSites();
  if (n < 2) {
    throw std::runtime_error("clusterChunk: need at least 2 haploids");
  }
  if (to == 0u) to = static_cast<unsigned>(total_sites);
  if (from >= to || static_cast<int>(to) > total_sites) {
    throw std::runtime_error("clusterChunk: invalid [from, to)");
  }
  const int chunk_first = static_cast<int>(from);
  const int chunk_last  = static_cast<int>(to);
  const int chunk_len   = chunk_last - chunk_first;

  if (site_block == 0u) site_block = static_cast<unsigned>(chunk_len);
  const int sb = static_cast<int>(std::min<unsigned>(site_block,
                                                     static_cast<unsigned>(chunk_len)));
  const int K = std::max(1u, anchor_stride_k);

  const std::int64_t num_pairs =
      static_cast<std::int64_t>(n) * (n - 1) / 2;

  ClusterChunkResult out;
  out.num_haps = n;
  out.num_sites = chunk_len;
  out.site_boundaries.resize(static_cast<std::size_t>(chunk_len + 1));
  // Outputs are pre-sized so OMP threads can write into known offsets
  // without locks. Each site contributes (n-1) nodes + 2(n-1) edges.
  out.node_time.assign(static_cast<std::size_t>(chunk_len) * (n - 1), 0.0);
  out.edge_left.assign(static_cast<std::size_t>(chunk_len) * 2 * (n - 1), 0.0);
  out.edge_right.assign(static_cast<std::size_t>(chunk_len) * 2 * (n - 1), 0.0);
  out.edge_parent.assign(static_cast<std::size_t>(chunk_len) * 2 * (n - 1), 0);
  out.edge_child.assign(static_cast<std::size_t>(chunk_len) * 2 * (n - 1), 0);

  // Site boundaries: interior bins are midpoints of adjacent physical
  // positions shifted by arg_offset. The first boundary is 0 and the
  // last boundary is arg_length, mirroring arg-needle's upgma_chunk
  // convention so chunk-open / chunk-close edges anchor to the ARG
  // extents instead of the first/last marker. If arg_length == 0 the
  // caller wants legacy physical-anchoring behaviour.
  const auto phys = asmc.getPhysicalPositions();
  if (arg_length > 0.0) {
    out.site_boundaries[0] = 0.0;
    out.site_boundaries[chunk_len] = arg_length;
  } else {
    out.site_boundaries[0] = static_cast<double>(phys[chunk_first] - arg_offset);
    out.site_boundaries[chunk_len] =
        static_cast<double>(phys[chunk_last - 1] - arg_offset);
  }
  for (int s = 1; s < chunk_len; ++s) {
    const int gs = chunk_first + s;
    out.site_boundaries[s] = 0.5 * (
        static_cast<double>(phys[gs - 1] - arg_offset) +
        static_cast<double>(phys[gs]     - arg_offset));
  }

  // Anchor map: K==1 → every site is its own anchor (bit-identical to
  // full per-site UPGMA). K>1 → anchors at chunk-local indices
  // {0, K, 2K, …, chunk_len-1}; every other site copies the *nearest*
  // anchor by chunk-local site index, ties broken left.
  std::vector<int> anchor_local_idx;
  std::vector<int> anchor_for_site(static_cast<std::size_t>(chunk_len), -1);
  for (int s = 0; s < chunk_len; ++s) {
    if (s % K == 0 || s == chunk_len - 1) {
      anchor_for_site[s] = static_cast<int>(anchor_local_idx.size());
      anchor_local_idx.push_back(s);
    }
  }
  const int n_anchors = static_cast<int>(anchor_local_idx.size());
  if (n_anchors == 0) {
    throw std::runtime_error("clusterChunk: empty chunk");
  }
  // Round-to-nearest anchor for every site (ties → left).
  std::vector<int> nearest_anchor(static_cast<std::size_t>(chunk_len), 0);
  for (int s = 0; s < chunk_len; ++s) {
    // Binary-search anchor_local_idx for s; ties broken left.
    auto it = std::lower_bound(anchor_local_idx.begin(),
                               anchor_local_idx.end(), s);
    int hi_idx = static_cast<int>(it - anchor_local_idx.begin());
    int chosen;
    if (hi_idx == 0) {
      chosen = 0;
    } else if (hi_idx == n_anchors) {
      chosen = n_anchors - 1;
    } else {
      const int lo_site = anchor_local_idx[hi_idx - 1];
      const int hi_site = anchor_local_idx[hi_idx];
      // Tie (s - lo == hi - s) → pick LEFT (lo)
      chosen = ((s - lo_site) <= (hi_site - s)) ? (hi_idx - 1) : hi_idx;
    }
    nearest_anchor[s] = chosen;
  }

  // Per-anchor storage. For K==1 we still allocate (n_anchors == chunk_len);
  // memory matches the previous single-pass behaviour. For K>1 this is
  // (chunk_len/K) × (n - 1) × 8 B for heights and a similar (n - 1) × 8 B
  // for the topology pairs — much smaller than the (n - 1) × chunk_len of
  // node_time we still emit.
  struct AnchorBuf {
    std::vector<double> heights;            // size n - 1
    std::vector<std::pair<int,int>> topo;   // size n - 1, scipy-linkage IDs
  };
  std::vector<AnchorBuf> anchor_data(static_cast<std::size_t>(n_anchors));
  for (auto& a : anchor_data) {
    a.heights.assign(static_cast<std::size_t>(n - 1), 0.0);
    a.topo.assign(static_cast<std::size_t>(n - 1), std::make_pair(0, 0));
  }

  // Per-OMP-thread scratch for the per-site UPGMA. n_threads == 0 means
  // the caller hasn't asked for parallelism so we use a single thread
  // (the per-block decode_pairs above also stays sequential in that case).
  const int upgma_threads = std::max(1, n_threads);
  std::vector<std::vector<double>> D_per_thread(
      upgma_threads, std::vector<double>(num_pairs));
  std::vector<std::vector<double>> Z_per_thread(
      upgma_threads,
      std::vector<double>(static_cast<std::size_t>(n - 1) * 4));

  // All-pairs decode args, in scipy condensed order so the decoder's
  // perPairPosteriorMeans rows align 1:1 with our pair indices and we
  // skip a per-iteration scatter step.
  std::vector<unsigned long> all_a(static_cast<std::size_t>(num_pairs));
  std::vector<unsigned long> all_b(static_cast<std::size_t>(num_pairs));
  {
    std::int64_t k = 0;
    for (int j = 0; j < n - 1; ++j) {
      for (int i = j + 1; i < n; ++i) {
        all_a[static_cast<std::size_t>(k)] = static_cast<unsigned long>(i);
        all_b[static_cast<std::size_t>(k)] = static_cast<unsigned long>(j);
        ++k;
      }
    }
  }

  // Force ASMC to surface posterior means; clusterChunk doesn't need MAP.
  asmc.setStorePerPairPosteriorMean(true);
  asmc.setStorePerPairMap(false);
  asmc.setStorePerPairPosterior(false);
  asmc.setStoreSumOfPosterior(false);

  for (int s_local = 0; s_local < chunk_len; s_local += sb) {
    const int b_size = std::min(sb, chunk_len - s_local);
    const int s_global_lo = chunk_first + s_local;
    const int s_global_hi = s_global_lo + b_size;

    // One decode_pairs call per site-block, batching all n*(n-1)/2 pairs
    // together. ASMC's batched decoder is far cheaper than (n-1)
    // separate calls — same forward-backward total work, but per-call
    // setup amortised over the whole pair batch.
    if (n_threads > 0) {
      asmc.decodePairsParallel(all_a, all_b,
                               static_cast<unsigned>(s_global_lo),
                               static_cast<unsigned>(s_global_hi),
                               cm_burn_in, n_threads);
    } else {
      asmc.decodePairs(all_a, all_b,
                       static_cast<unsigned>(s_global_lo),
                       static_cast<unsigned>(s_global_hi),
                       cm_burn_in);
    }
    const auto& res = asmc.getRefOfResults();
    const auto& M = res.perPairPosteriorMeans;  // (num_pairs, b_size)
    // Pairs are emitted in scipy condensed order so M's row k lines up
    // with the k-th condensed pair index. PASS 1 reads M directly per
    // anchor site (one column per UPGMA call) — no D_block transpose,
    // saves numPairs × b_size × 8 B of intermediate storage (e.g.
    // 47 GB at n=2400, sb=2048).

    // PASS 1 — UPGMA only at anchor sites in this block, OMP-parallel.
    // Each thread writes into its own AnchorBuf slot (disjoint by anchor
    // index), so no synchronisation is needed.
#ifdef _OPENMP
    bool fastcluster_failed = false;
    bool zero_distance_seen = false;
#pragma omp parallel for num_threads(upgma_threads) schedule(static) \
    shared(fastcluster_failed, zero_distance_seen)
#endif
    for (int s = 0; s < b_size; ++s) {
      const int site_local = s_local + s;
      const int ai = anchor_for_site[site_local];
      if (ai < 0) continue;  // not an anchor → skip UPGMA
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
#else
      const int tid = 0;
#endif
      auto& D_one = D_per_thread[tid];
      auto& Z_one = Z_per_thread[tid];

      // Copy column s of M (Eigen Matrix<float, Dynamic, Dynamic, ColMajor>;
      // contiguous in column) into D_one as double, in scipy condensed
      // order. Each thread does this for a different site so cache lines
      // don't bounce.
      const auto col = M.col(s);
      for (std::int64_t k = 0; k < num_pairs; ++k) {
        D_one[static_cast<std::size_t>(k)] =
            static_cast<double>(col(static_cast<Eigen::Index>(k)));
      }
      const int rc = asmc_fastcluster_average_linkage(n, D_one.data(),
                                                      Z_one.data());
      if (rc != 0) {
#ifdef _OPENMP
#pragma omp atomic write
        fastcluster_failed = true;
        continue;
#else
        throw std::runtime_error("clusterChunk: fastcluster failed at site");
#endif
      }
      // Store linkage row-by-row (heights + (left_id, right_id) pairs in
      // scipy convention: ids in [0, n) are samples, ids in [n, 2n-1) are
      // internal node indices created in earlier rows).
      auto& abuf = anchor_data[static_cast<std::size_t>(ai)];
      double last_distance = 0.0;
      bool local_zero = false;
      for (int k = 0; k < n - 1; ++k) {
        const double* row = Z_one.data() + static_cast<std::size_t>(k) * 4;
        double height = row[2];
        if (height == 0.0) { local_zero = true; break; }
        if (height <= last_distance) height = last_distance * 1.00001;
        last_distance = height;
        abuf.heights[static_cast<std::size_t>(k)] = height;
        abuf.topo[static_cast<std::size_t>(k)] = std::make_pair(
            static_cast<int>(row[0]), static_cast<int>(row[1]));
      }
      if (local_zero) {
#ifdef _OPENMP
#pragma omp atomic write
        zero_distance_seen = true;
#else
        throw std::runtime_error("clusterChunk: zero merge distance");
#endif
      }
    }

#ifdef _OPENMP
    if (fastcluster_failed) {
      throw std::runtime_error("clusterChunk: fastcluster failed at a site");
    }
    if (zero_distance_seen) {
      throw std::runtime_error("clusterChunk: zero merge distance at a site");
    }
#endif

  } // end per-block PASS 1 loop

  // PASS 2 — emit all sites in the chunk now that every anchor is filled.
  // Round-to-nearest anchor (computed above as `nearest_anchor`) means a
  // site can need an anchor that lives in a *later* block, so we defer
  // emission to here rather than running it inline per block. OMP-parallel
  // since each site writes to disjoint output slabs.
#ifdef _OPENMP
#pragma omp parallel for num_threads(upgma_threads) schedule(static)
#endif
  for (int site_local = 0; site_local < chunk_len; ++site_local) {
    const int ai = nearest_anchor[site_local];
    const auto& abuf = anchor_data[static_cast<std::size_t>(ai)];

    const double left  = out.site_boundaries[site_local];
    const double right = out.site_boundaries[site_local + 1];
    const int node_offset = site_local * (n - 1);
    const std::size_t node_base = static_cast<std::size_t>(site_local) * (n - 1);
    const std::size_t edge_base = node_base * 2;

    for (int k = 0; k < n - 1; ++k) {
      out.node_time[node_base + static_cast<std::size_t>(k)] =
          abuf.heights[static_cast<std::size_t>(k)];
      const int parent = n + node_offset + k;
      const auto& tp = abuf.topo[static_cast<std::size_t>(k)];
      for (int c = 0; c < 2; ++c) {
        const int child_id = (c == 0) ? tp.first : tp.second;
        const int child_node = (child_id < n)
            ? child_id
            : n + node_offset + (child_id - n);
        const std::size_t out_idx = edge_base
            + static_cast<std::size_t>(k) * 2
            + static_cast<std::size_t>(c);
        out.edge_left[out_idx]   = left;
        out.edge_right[out_idx]  = right;
        out.edge_parent[out_idx] = parent;
        out.edge_child[out_idx]  = child_node;
      }
    }
  }

  return out;
}

} // namespace ASMC
