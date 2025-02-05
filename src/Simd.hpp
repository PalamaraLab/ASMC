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

namespace asmc
{

/**
 * Gets the number of SIMD float lanes available, detected at runtime.
 *
 * @return an int representing the number of SIMD float lanes availalble
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

} // namespace asmc

#endif // ASMC_SIMD
