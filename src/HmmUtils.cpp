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

#include "HmmUtils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>

#include <fmt/format.h>

namespace asmc
{

std::vector<bool> subsetXorVec(const std::vector<bool>& v1, const std::vector<bool>& v2, const unsigned long from,
                               const unsigned long to) noexcept
{
  const auto min_to = std::min<unsigned long>(v1.size(), to);

  assert(v1.size() == v2.size());
  assert(from < min_to);

  std::vector<bool> ret(min_to - from);

  for (unsigned long i = from; i < min_to; ++i) {
    ret[i - from] = v1[i] ^ v2[i];
  }

  return ret;
}

std::vector<bool> subsetAndVec(const std::vector<bool>& v1, const std::vector<bool>& v2, const unsigned long from,
                               const unsigned long to) noexcept
{
  const auto min_to = std::min<unsigned long>(v1.size(), to);

  assert(v1.size() == v2.size());
  assert(from < min_to);

  std::vector<bool> ret(min_to - from);

  for (unsigned i = from; i < min_to; i++) {
    ret[i - from] = v1[i] & v2[i];
  }

  return ret;
}

float roundMorgans(const float value, const int precision, const float min) noexcept
{
  assert(precision >= 0);
  assert(min > 0.f);

  if (value <= min) {
    return min;
  }

  const float correction = 10.f - static_cast<float>(precision);
  const float L10 = std::max<float>(0.f, floorf(log10f(value)) + correction);
  const float factor = powf(10.f, 10.f - L10);

  return roundf(value * factor) / factor;
}

int roundPhysical(const int value, const int precision) noexcept
{
  assert(precision >= 0);
  assert(value > -2); // Since HMM for sequence uses distance-1, it can be -1

  if (value <= 1) {
    return 1;
  }

  int L10 = std::max<int>(0, static_cast<int>(floor(log10(value))) - precision);
  int factor = static_cast<int>(pow(10, L10));

  return static_cast<int>(round(value / static_cast<double>(factor))) * factor;
}

void printPctTime(const std::string& str, double fracTime, std::ostream& os)
{
  os << "Time in " << std::setw(14) << std::left << str << " : " << std::setw(5) << std::right
     << std::setprecision(1) << std::fixed << 100.0 * fracTime << "%\n" << std::flush;
}

unsigned getFromPosition(const std::vector<float>& geneticPositions, unsigned from, const float cmDist)
{
  assert(cmDist > 0.f);
  assert(geneticPositions.size() > from);

  float cumGenDist = 0.f;
  while (cumGenDist < cmDist && from > 0u) {
    from--;
    cumGenDist += (geneticPositions[from + 1u] - geneticPositions[from]) * 100.f;
  }
  return from;
}

unsigned getToPosition(const std::vector<float>& geneticPositions, unsigned to, const float cmDist)
{
  assert(cmDist > 0.f);
  assert(geneticPositions.size() >= to);

  float cumGenDist = 0.f;
  while (cumGenDist < cmDist && to + 1u < geneticPositions.size()) {
    to++;
    cumGenDist += (geneticPositions[to] - geneticPositions[to - 1u]) * 100.f;
  }
  return std::min<unsigned>(to + 1u, static_cast<unsigned>(geneticPositions.size()));
}

std::pair<unsigned long, unsigned long> hapToDipId(unsigned long hapId)
{
  return std::make_pair(hapId / 2ul, 1ul + (hapId % 2ul));
}

unsigned long dipToHapId(unsigned long ind, unsigned long hap)
{
  assert(hap > 0ul);
  return 2ul * ind + hap - 1ul;
}

std::string indPlusHapToCombinedId(std::string_view indId, unsigned long hap)
{
  if (indId.empty() || !(hap == 1ul || hap == 2ul)) {
    throw std::runtime_error(
        fmt::format("Expected an individual ID and either 1 or 2, but got {} and {}\n", indId, hap));
  }

  return fmt::format("{}_{}", indId, hap);
}

std::pair<std::string, unsigned long> combinedIdToIndPlusHap(std::string_view combinedId)
{
  if (combinedId.length() < 3 || !(combinedId.substr(combinedId.length() - 2, 2) == "_1" ||
                                   combinedId.substr(combinedId.length() - 2, 2) == "_2")) {
    throw std::runtime_error(fmt::format("Expected combined ID in form <id>_1 OR <id>_2, but got {}\n", combinedId));
  }
  return std::make_pair(std::string{combinedId.substr(0, combinedId.length() - 2)},
                        combinedId.back() == '1' ? 1ul : 2ul);
}

unsigned long getIndIdxFromIdString(const std::vector<std::string>& idStrings, std::string_view idString)
{
  auto it = std::find(idStrings.begin(), idStrings.end(), idString);
  if (it == idStrings.end()) {
    throw std::runtime_error(fmt::format("The ID string {} is not in the list of IDs\n", idString));
  }
  return std::distance(idStrings.begin(), it);
}

} // namespace asmc
