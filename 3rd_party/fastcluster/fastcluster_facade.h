// Thin C-style facade over Daniel Müllner's fastcluster library, vendored
// under ASMC_dev/3rd_party/fastcluster/ for ASMC::clusterChunk. Exposes
// only the entry points ASMC actually needs (currently: average linkage on
// a precomputed dissimilarity matrix); skips the Python / R glue.
//
// fastcluster is dual-licensed BSD-2 OR GPL-2-or-later (see COPYING.txt
// in this directory); ASMC vendors under BSD-2.

#ifndef ASMC_FASTCLUSTER_FACADE_H
#define ASMC_FASTCLUSTER_FACADE_H

#ifdef __cplusplus
extern "C" {
#endif

// average_linkage — UPGMA on a condensed distance matrix.
//
// Inputs:
//   N            number of items being clustered (≥ 2).
//   D            condensed distance matrix in scipy convention: a flat
//                array of N*(N-1)/2 doubles, where pair (i, j) with i<j
//                lives at index N*(N-1)/2 - (N-i)*(N-i-1)/2 + (j-i-1).
//                The buffer is mutated in-place by the algorithm; copy
//                first if you need to preserve it.
//
// Output:
//   Z_out        scipy-style linkage matrix, (N-1) × 4 doubles in row-
//                major order. Row k = [a, b, dist, count]:
//                  - a, b: cluster IDs being merged (smallest first)
//                  - dist: merge distance
//                  - count: number of original items in the new cluster
//                Cluster IDs in [0, N) are original items; IDs in
//                [N, 2N-1) are merged clusters created in earlier rows.
//                Output is sorted ascending by `dist`, matching
//                scipy.cluster.hierarchy.linkage(method='average').
//
// Returns:
//   0 on success, non-zero on error (e.g. NaN dissimilarity).
int asmc_fastcluster_average_linkage(int N, double* D, double* Z_out);

#ifdef __cplusplus
}
#endif

#endif // ASMC_FASTCLUSTER_FACADE_H
