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

#ifndef ASMC_HASHING_UTILS_HPP
#define ASMC_HASHING_UTILS_HPP

#include <vector>

namespace asmc
{

/**
 * Convenience function to compute genetic distance between two words (start of w1 and end of w2)
 *
 * @param w1 the first word
 * @param w2 the second word
 * @param geneticPositions vector of genetic positions
 * @param wordSize number of locations per word
 * @return the number of centimorgans between start of w1 and end of w2
 */
double cmBetween(int w1, int w2, const std::vector<float>& geneticPositions, int wordSize);

} // namespace asmc

#endif // ASMC_HASHING_UTILS_HPP
