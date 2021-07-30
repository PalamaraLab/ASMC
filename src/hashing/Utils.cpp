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

#include "hashing/Utils.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

double asmc::cmBetween(const int w1, const int w2, const std::vector<float>& geneticPositions, const int wordSize)
{
  assert(!geneticPositions.empty());
  assert(wordSize * w1 < geneticPositions.size());
  assert(w1 >= 0);
  assert(w2 >= 0);
  assert(w2 >= w1);

  const std::size_t start = wordSize * w1;
  const std::size_t end = std::min<std::size_t>(wordSize * w2 + wordSize - 1, geneticPositions.size() - 1ul);

  return 100.0 * (geneticPositions[end] - geneticPositions[start]);
}
