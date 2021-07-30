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

#ifndef ASMC_HASHING_INDIVIDUALS_HPP
#define ASMC_HASHING_INDIVIDUALS_HPP

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "boost/dynamic_bitset.hpp"

class Individuals
{
  unsigned mIdNum;
  unsigned long mWordSize = 64ul;
  unsigned long mNumReadAhead = 10ul;

  std::vector<boost::dynamic_bitset<>> mHap{mNumReadAhead, boost::dynamic_bitset<>(mWordSize, 0ul)};

public:
  explicit Individuals(const unsigned long wordSize, const unsigned long numReadAhead, const unsigned idNum)
      : mIdNum{idNum}, mWordSize{wordSize}, mNumReadAhead{numReadAhead}
  {
    assert(wordSize > 0ul);
    assert(numReadAhead > 0ul);

    mHap.resize(numReadAhead);
    std::fill(mHap.begin(), mHap.end(), boost::dynamic_bitset<>(wordSize, 0ul));
  }

  void clear(const int w)
  {
    assert(w >= 0);
    mHap.at(w % mNumReadAhead).reset();
  }

  void setMarker(const int w, const std::size_t bit)
  {
    assert(w >= 0);
    assert(bit < mWordSize);
    mHap.at(w % mNumReadAhead).set(bit);
  }

  unsigned long getWordHash(const int w)
  {
    assert(w >= 0);
    return mHap.at(w % mNumReadAhead).to_ulong();
  }

  std::string getWordString(const int w)
  {
    assert(w >= 0);
    std::string buffer;
    boost::to_string(mHap.at(w % mNumReadAhead), buffer);
    return buffer;
  }

  unsigned int getIdNum() const
  {
    return mIdNum;
  }

  unsigned long getWordSize() const
  {
    return mHap.front().size();
  }

  unsigned long getNumReadAhead() const
  {
    return mHap.size();
  }
};

#endif // ASMC_HASHING_INDIVIDUALS_HPP
