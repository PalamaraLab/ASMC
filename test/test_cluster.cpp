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

// Unit tests for ASMC::clusterChunk (the in-process ASMC-clust kernel).
// These run on the small bundled n=300 array fixture to keep CI fast;
// large-n properties (memory, parallel scaling) are exercised in the
// out-of-tree regression suite at tasks/asmc_clust_native/.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <unordered_set>

#include "ASMC.hpp"
#include "Cluster.hpp"

namespace {

// Build a small ASMC instance from the bundled array fixture and run
// clusterChunk over a [from, to) range. The returned struct has
// (chunk_len = to - from) sites, (n - 1) internal nodes per site, and
// 2*(n - 1) edges per site.
ASMC::ClusterChunkResult runChunk(unsigned from, unsigned to, unsigned k,
                                  unsigned site_block = 0u) {
  ASMC::ASMC asmc(ASMC_DATA_DIR "/examples/asmc/exampleFile.n300.array",
                  ASMC_DATA_DIR "/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz");
  return ASMC::clusterChunk(asmc, /*from=*/from, /*to=*/to,
                            /*site_block=*/site_block,
                            /*n_threads=*/0, /*cm_burn_in=*/0.5f,
                            /*arg_offset=*/0, /*arg_length=*/0.0,
                            /*anchor_stride_k=*/k);
}

} // namespace

TEST_CASE("clusterChunk: layout invariants at K=1", "[Cluster]")
{
  const unsigned from = 0u;
  const unsigned to = 8u;  // 8 sites is enough to exercise the per-site loop without long wallclock
  auto r = runChunk(from, to, /*k=*/1u);

  const int n = r.num_haps;
  // The fixture has 150 diploid samples → 300 haploids; the file is named
  // "n300" because the haploid count is 300 (consistent with the rest of
  // the ASMC test suite).
  REQUIRE(n == 300);
  REQUIRE(r.num_sites == static_cast<int>(to - from));

  // Expected sizes: per site, (n - 1) internal nodes + 2*(n - 1) edges.
  const std::size_t per_site_nodes = static_cast<std::size_t>(n - 1);
  const std::size_t per_site_edges = 2 * per_site_nodes;
  REQUIRE(r.node_time.size() == per_site_nodes * r.num_sites);
  REQUIRE(r.edge_left.size() == per_site_edges * r.num_sites);
  REQUIRE(r.edge_right.size() == per_site_edges * r.num_sites);
  REQUIRE(r.edge_parent.size() == per_site_edges * r.num_sites);
  REQUIRE(r.edge_child.size() == per_site_edges * r.num_sites);
  REQUIRE(r.site_boundaries.size() == static_cast<std::size_t>(r.num_sites + 1));
}

TEST_CASE("clusterChunk: heights are strictly increasing per site", "[Cluster]")
{
  // Strict monotonicity is enforced by the *= 1.00001 nudge inside
  // clusterChunk (matches arg-needle's upgma_chunk behaviour). Equal-height
  // ties from fastcluster would otherwise break tskit's parent-after-child
  // time invariant.
  auto r = runChunk(/*from=*/0u, /*to=*/4u, /*k=*/1u);
  const int n = r.num_haps;
  for (int s = 0; s < r.num_sites; ++s) {
    const std::size_t base = static_cast<std::size_t>(s) * (n - 1);
    double prev = -1.0;
    for (int k = 0; k < n - 1; ++k) {
      const double h = r.node_time[base + static_cast<std::size_t>(k)];
      CHECK(h > prev);
      prev = h;
    }
  }
}

TEST_CASE("clusterChunk: edges reference valid nodes", "[Cluster]")
{
  auto r = runChunk(/*from=*/0u, /*to=*/4u, /*k=*/1u);
  const int n = r.num_haps;
  const int total_internal = r.num_internal_nodes();
  // Children: any sample [0, n) or any internal node already created
  // before this row in the per-site sweep. Parents: only internal IDs.
  for (std::size_t i = 0; i < r.edge_parent.size(); ++i) {
    const int parent = r.edge_parent[i];
    const int child = r.edge_child[i];
    CHECK(parent >= n);
    CHECK(parent < n + total_internal);
    CHECK(child >= 0);
    CHECK(child < n + total_internal);
    CHECK(child != parent);
  }
}

TEST_CASE("clusterChunk: site boundaries are monotonic", "[Cluster]")
{
  auto r = runChunk(/*from=*/0u, /*to=*/8u, /*k=*/1u);
  for (int s = 0; s < r.num_sites; ++s) {
    CHECK(r.site_boundaries[s] < r.site_boundaries[s + 1]);
  }
}

TEST_CASE("clusterChunk: anchor_stride_k>1 emits the same site count",
          "[Cluster]")
{
  // K=1 = full UPGMA at every site. K=2 = UPGMA at every other site;
  // remaining sites copy the nearest anchor's topology and heights.
  // The output table must still cover every site (same num_sites); only
  // the underlying merge work is reduced.
  auto r1 = runChunk(/*from=*/0u, /*to=*/8u, /*k=*/1u);
  auto r2 = runChunk(/*from=*/0u, /*to=*/8u, /*k=*/2u);
  REQUIRE(r1.num_sites == r2.num_sites);
  REQUIRE(r1.node_time.size() == r2.node_time.size());
  REQUIRE(r1.edge_parent.size() == r2.edge_parent.size());
}

TEST_CASE("clusterChunk: site_block bounds resident memory but yields same shape",
          "[Cluster]")
{
  // sb=0 is "one big block over the whole chunk" (cheapest wallclock,
  // biggest RAM). sb>0 splits into per-block decode_pairs calls. The
  // per-site outputs (counts, monotonicity) must be identical in shape.
  auto r_big   = runChunk(/*from=*/0u, /*to=*/8u, /*k=*/1u, /*sb=*/0u);
  auto r_small = runChunk(/*from=*/0u, /*to=*/8u, /*k=*/1u, /*sb=*/2u);
  REQUIRE(r_big.num_sites == r_small.num_sites);
  REQUIRE(r_big.node_time.size() == r_small.node_time.size());
}
