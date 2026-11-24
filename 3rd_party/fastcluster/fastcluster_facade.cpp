// fastcluster facade — single TU that owns the algorithm symbols.
//
// We #include fastcluster.cpp directly (no header guard) so the
// file-static templates (NN_chain_core, generic_linkage) and helper
// classes (cluster_result, union_find, auto_array_ptr) are visible to
// the wrapper code below. This is unusual but matches the upstream
// pattern, which inlines the same source into fastcluster_python.cpp.
//
// The wrapper exposes a clean extern "C" entry point so callers can
// link against this object without seeing any of the Müllner internals.

// fastcluster.cpp uses fc_isnan but does not define it; upstream's
// fastcluster_python.cpp defines it before #include "fastcluster.cpp".
// We mirror that here so the algorithm sources see a definition.
#define fc_isnan(X) ((X) != (X))

#include "fastcluster.cpp"  // defines NN_chain_core, cluster_result, etc.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace asmc_fc {

// Lifted from fastcluster_python.cpp: convenience wrapper around the
// (N-1) × 4 output buffer with a column-major-ish append helper.
class linkage_output {
 private:
  t_float* Z;

 public:
  explicit linkage_output(t_float* const Z_) : Z(Z_) {}

  void append(const t_index node1, const t_index node2,
              const t_float dist, const t_float size) {
    if (node1 < node2) {
      *(Z++) = static_cast<t_float>(node1);
      *(Z++) = static_cast<t_float>(node2);
    } else {
      *(Z++) = static_cast<t_float>(node2);
      *(Z++) = static_cast<t_float>(node1);
    }
    *(Z++) = dist;
    *(Z++) = size;
  }
};

// size_(r) returns the cluster-size of cluster id `r`. r in [0, N) → 1
// (original sample), r ≥ N → looked up from the linkage matrix's count
// column for the cluster created in row r-N.
//
// We use an inline helper instead of the upstream macro so it doesn't
// shadow scope outside this translation unit.
template <typename T>
inline T cluster_size(const T r, const T N, const t_float* const Z) {
  return (r < N) ? T{1} : static_cast<T>(Z[(r - N) * 4 + 3]);
}

// Generate the SciPy-format dendrogram from a cluster_result. Sorts
// in-place by merge distance (stable_sort, matches upstream sorted=false
// path used for NN_chain output). Identical algorithm to
// generate_SciPy_dendrogram<false> in fastcluster_python.cpp.
static void generate_scipy_dendrogram(t_float* const Z,
                                      cluster_result& Z2,
                                      const t_index N) {
  union_find nodes(N);
  std::stable_sort(Z2[0], Z2[N - 1]);

  linkage_output output(Z);
  for (node const* NN = Z2[0]; NN != Z2[N - 1]; ++NN) {
    const t_index node1 = nodes.Find(NN->node1);
    const t_index node2 = nodes.Find(NN->node2);
    nodes.Union(node1, node2);
    const t_float size = cluster_size(node1, N, Z) + cluster_size(node2, N, Z);
    output.append(node1, node2, NN->dist, size);
  }
}

}  // namespace asmc_fc

extern "C" int asmc_fastcluster_average_linkage(int N, double* D, double* Z_out) {
  if (N < 2 || D == nullptr || Z_out == nullptr) {
    return 1;
  }
  try {
    cluster_result Z2(static_cast<t_index>(N - 1));
    auto_array_ptr<t_index> members;
    members.init(static_cast<t_index>(N), 1);
    NN_chain_core<METHOD_METR_AVERAGE, t_index>(
        static_cast<t_index>(N), D, members, Z2);
    asmc_fc::generate_scipy_dendrogram(Z_out, Z2, static_cast<t_index>(N));
    return 0;
  } catch (const nan_error&) {
    return 2;  // NaN distance encountered
  } catch (const std::bad_alloc&) {
    return 3;
  } catch (...) {
    return 4;
  }
}
