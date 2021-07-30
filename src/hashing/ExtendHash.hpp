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

#ifndef ASMC_HASHING_EXTEND_HASH_HPP
#define ASMC_HASHING_EXTEND_HASH_HPP

#include <utility>

#include <boost/unordered_map.hpp>

#include "hashing/Match.hpp"

/* Object for storing extension between pairs of Individuals */
class ExtendHash
{

  boost::unordered_map<unsigned long int, Match> extend_hash;

  unsigned long mWordSize;
  unsigned long num;

  bool mParHaploid;

  // Empty Match to insert into hash
  Match m;

  // Iterator for testing insertion
  std::pair<boost::unordered::iterator_detail::iterator<
                boost::unordered::detail::ptr_node<std::pair<const unsigned long, Match>>>,
            bool>
      extend_ret;

public:
  explicit ExtendHash(const unsigned long wordSize, const unsigned long num, const bool PAR_HAPLOID)
      : mWordSize(wordSize), num(num), mParHaploid(PAR_HAPLOID), m(wordSize)
  {
  }

  // Compute pair of Individuals from location indicator
  std::pair<unsigned, unsigned> locationToPair(unsigned long loc)
  {
    const unsigned second = mParHaploid ? loc % num : 2 * (loc % num);
    const unsigned first = mParHaploid ? (loc - second) / num : 2 * ((loc - second / 2) / num);

    return std::make_pair(first, second);
  }

  // Compute location from pair of Individuals
  unsigned long pairToLocation(unsigned int i, unsigned int j)
  {
    if (!mParHaploid) {
      // round everyone down to the nearest haplotype
      i = (i - (i % 2)) / 2;
      j = (j - (j % 2)) / 2;
    }
    unsigned long loc = (i > j) ? j * num + i : i * num + j;
    return loc;
  }

  // Extend or add a given pair in the current hash
  // unsigned int i,j : identifiers for the two Individuals
  // int w : current word # to extend or add
  void extendPair(unsigned int i, unsigned int j, int w, const int GLOBAL_CURRENT_WORD)
  {
    m.getModifiableInterval()[0] = GLOBAL_CURRENT_WORD;
    // Find/extend this location in the hash
    extend_ret = extend_hash.insert(std::pair<unsigned long int, Match>(pairToLocation(i, j), m));
    (extend_ret.first->second).extend(w);
  }

  // Remove all pairs that were not extended beyond w
  // int w : word # to remove prior to
  void clearPairsPriorTo(int w, const int GLOBAL_CURRENT_WORD, const double PAR_MIN_MATCH,
                         const std::vector<float>& geneticPositions, HMM& hmm)
  {
    for (auto it = extend_hash.begin(); it != extend_hash.end();) {
      if (it->second.getInterval()[1] < w) {
        it->second.print(locationToPair(it->first), PAR_MIN_MATCH, geneticPositions, hmm);
        it = extend_hash.erase(it);
      } else {
        if (it->second.getInterval()[1] < GLOBAL_CURRENT_WORD)
          it->second.addGap();
        it++;
      }
    }
  }

  // Remove all pairs that were not extended beyond w
  // int w : word # to remove prior to
  void extendAllPairsTo(int w)
  {
    for (auto it = extend_hash.begin(); it != extend_hash.end(); it++)
      it->second.getModifiableInterval()[1] = w;
  }

  // Remove all pairs
  // int w : word # to remove prior to
  void clearAllPairs(const double PAR_MIN_MATCH, const std::vector<float>& geneticPositions, HMM& hmm)
  {
    for (auto it = extend_hash.begin(); it != extend_hash.end();) {
      it->second.print(locationToPair(it->first), PAR_MIN_MATCH, geneticPositions, hmm);
      it = extend_hash.erase(it);
    }
  }

  std::size_t size() const
  {
    return extend_hash.size();
  }

  unsigned long getWordSize() const
  {
    return mWordSize;
  }

};

#endif // ASMC_HASHING_EXTEND_HASH_HPP
