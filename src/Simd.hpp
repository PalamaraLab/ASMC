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

#ifndef ASMC_SIMD
#define ASMC_SIMD

#include <Eigen/Core>

namespace asmc
{

/**
 * Gets the number of SIMD float lanes available, detected at runtime.
 *
 * @return an int representing the number of SIMD float lanes available
 */
int getNumSimdLanes();

/**
 * Validates that the given batch size is a multiple of the SIMD lane width.
 *
 * @param batchSize the number of elements processed in a batch
 * @throws std::invalid_argument if batchSize is not a multiple of the SIMD lane width
 */
void validateBatchSize(int batchSize);

/**
 * Prints runtime SIMD information, including the number of float lanes
 * and the detected target SIMD instruction set.
 */
void printRuntimeSimdInfo();


/**
 * Calculate scaling factors for a contiguous array of curBatchSize * numStates floats. The computed scaling factors
 * are normalized across all states for each item in the batch.
 *
 * This method takes three buffers: one with data, one for storing calculated scaling factors, and one for temporarily
 * storing intermediate sums. It is assumed that all three buffers are appropriately allocated and deallocated outside
 * this function.
 *
 * The sums buffer is set to zeros in this method.
 *
 * @param vec an array of length curBatchSize * numStates containing the data
 * @param scalings an array of length curBatchSize to write calculated scaling factors into
 * @param sums an array of length curBatchSize for storing intermediate sums
 * @param batchSize the number of items in the batch
 * @param numStates the number of states over which to normalize
 */
void calculateScalingBatch(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings,
                           Eigen::Ref<Eigen::ArrayXf> sums, int batchSize, int numStates);

/**
 * Apply scaling factors to a contiguous array of curBatchSize * numStates floats. The scale factors are applied per
 * item to each state, and each state is separated in the data array by a stride length of batchSize.
 *
 * This method takes two buffers; one with data that is modified in-place, and one with scale factors that is just read
 * from.
 *
 * @param vec the buffer of length curBatchSize * numStates containing the data
 * @param scalings an array of length curBatchSize containing the prescribed scale factors
 * @param batchSize the number of items in the batch
 * @param numStates the number of states over which to normalize
 */
void applyScalingBatch(Eigen::Ref<Eigen::ArrayXf> vec, Eigen::Ref<Eigen::ArrayXf> scalings, int batchSize,
                       int numStates);

/**
 * Multiply two buffers elementwise (alpha *= beta), accumulate the result across states into a scaling buffer,
 * normalize the scaling buffer by computing its reciprocal, and apply it to the alpha buffer in-place.
 *
 * The alpha and beta buffers are of length batchSize * numStates * sequenceLength (flattened 3D layout). For each
 * position, the alpha and beta values for all states are multiplied, accumulated into a per-batch scaling factor, and
 * then alpha is normalized using the inverse of these factors.
 *
 * All arrays are accessed with strides of batchSize; alpha and beta are modified in-place, scale is both written and
 * read.
 *
 * @param alpha the buffer to modify in-place (length: batchSize * numStates * sequenceLength)
 * @param beta the buffer to multiply with alpha (same shape as alpha)
 * @param scale the scaling buffer (length: batchSize * sequenceLength), written to and used for normalization
 * @param batchSize the number of items in the batch
 * @param numStates the number of HMM states
 * @param from the start position (inclusive) along the sequence dimension
 * @param to the end position (exclusive) along the sequence dimension
 */
void normalizeAlphaWithBeta(Eigen::Ref<Eigen::ArrayXf> alpha, Eigen::Ref<Eigen::ArrayXf> beta,
                            Eigen::Ref<Eigen::ArrayXf> scale, int batchSize, int numStates, int from, int to);

/**
 * Update a column of the alpha buffer by adding the next column and the corresponding values in previousAlpha.
 *
 * The data is laid out row-major in a flattened 2D array of shape [numStates, batchSize]. This function updates
 * the k-th row (in-place) as: alphaC[k] = alphaC[k + 1] + previousAlpha[k].
 *
 * @param alphaC buffer of shape [numStates, batchSize] (flattened); modified in-place
 * @param previousAlpha buffer of same shape; read-only
 * @param batchSize number of items in a single state (stride size)
 * @param k index of the state to update
 */
void updateAlphaColumn(Eigen::Ref<Eigen::ArrayXf> alphaC, Eigen::Ref<Eigen::ArrayXf> previousAlpha, int batchSize,
                       int k);

} // namespace asmc

#endif // ASMC_SIMD
