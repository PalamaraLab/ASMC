// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef DATA_MODULE_VECTOR_UTILS_HPP
#define DATA_MODULE_VECTOR_UTILS_HPP

#include <algorithm>
#include <functional>
#include <vector>

namespace asmc {

/**
 * Determine whether a vector is strictly increasing.
 *
 * @tparam T the type over which the vector is templated
 * @param vec vector to test for monotonicity
 * @return whether the vector is strictly increasing
 */
template <typename T> bool isStrictlyIncreasing(std::vector<T> vec) {
  return std::adjacent_find(vec.begin(), vec.end(), std::greater_equal<T>()) == vec.end();
}

/**
 * Determine whether a vector is increasing.
 *
 * @tparam T the type over which the vector is templated
 * @param vec vector to test for monotonicity
 * @return whether the vector is increasing
 */
template <typename T> bool isIncreasing(std::vector<T> vec) {
  return std::adjacent_find(vec.begin(), vec.end(), std::greater<T>()) == vec.end();
}

} // namespace asmc

#endif // DATA_MODULE_VECTOR_UTILS_HPP
