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

#ifndef ASMC_HASHING_SEEDHASH_HPP
#define ASMC_HASHING_SEEDHASH_HPP

#include <utility>
#include <vector>

#include <boost/unordered_map.hpp>

#include "hashing/ExtendHash.hpp"
#include "hashing/Individuals.hpp"
#include "Types.hpp"

/* Object for storing initial word seeds */
class SeedHash
{

  using ind_vec = std::vector<Individuals>;

  boost::unordered_map<hash_size, std::vector<unsigned int>> seed_hash;
  // Empty vector to insert into the seed hash
  std::vector<unsigned int> vec;
  // Iterator for testing insertion of elements
  // std::pair<boost::unordered::iterator_detail::iterator<boost::unordered::detail::ptr_node<std::pair<const unsigned
  // int, vector<unsigned int> > > >, bool> seed_ret;
public:
  void insertIndividuals(unsigned int i, hash_size word)
  {
    auto seed_ret = seed_hash.insert(std::pair<hash_size, std::vector<unsigned int>>(word, vec));
    (seed_ret.first->second).push_back(i);
  }
  void clear()
  {
    seed_hash.clear();
  }
  int size()
  {
    return seed_hash.size();
  }

  // Generate a new hash for this vector of Individualss
  static unsigned long subHash(ExtendHash* e, std::vector<unsigned> v, int w, ind_vec all_ind, const int MAX_seeds,
                               const int jobID, const int jobs, const unsigned w_i, const unsigned w_j,
                               const unsigned windowSize, const int GLOBAL_READ_WORDS, int& GLOBAL_SKIPPED_WORDS,
                               const int GLOBAL_CURRENT_WORD, const bool is_j_above_diag)
  {
    SeedHash cur_sh;
    // seed the next word from this subset of Individualss
    for (unsigned int& i : v) {
      cur_sh.insertIndividuals(i, all_ind[i].getWordHash(w));
    }
    // recursion:
    return cur_sh.extendAllPairs(e, w, all_ind, MAX_seeds, jobID, jobs, w_i, w_j, windowSize, GLOBAL_READ_WORDS,
                                 GLOBAL_SKIPPED_WORDS, GLOBAL_CURRENT_WORD, is_j_above_diag);
  }

  // Extend/save all pairs in the current hash
  // ExtendHash * e : Pointer to ExtendHash which will be called for each pair
  // returns : number of pairs evaluated
  unsigned long extendAllPairs(ExtendHash* e, int w, ind_vec all_ind, const int MAX_seeds, const int jobID,
                               const int jobs, const unsigned w_i, const unsigned w_j, const unsigned windowSize,
                               const int GLOBAL_READ_WORDS, int& GLOBAL_SKIPPED_WORDS, const int GLOBAL_CURRENT_WORD,
                               const bool is_j_above_diag)
  {
    unsigned long tot_pairs = 0;
    for (auto it = seed_hash.begin(); it != seed_hash.end(); ++it) {

      // *** As long as the # of pairs is high, generate a sub-hash for the next word
      // *** Only store pairs of Individuals that have collision in a small hash
      // *** Extend only to the haplotypes that seeded here
      if (MAX_seeds != 0 && it->second.size() > static_cast<unsigned long>(MAX_seeds) && w + 1 < GLOBAL_READ_WORDS) {
        // recursively generate a sub-hash
        // IMPORTANT: if we run out of buffered words then this seed does not get analyzed
        if (w + 1 < GLOBAL_READ_WORDS) {
          tot_pairs += subHash(e, it->second, w + 1, all_ind, MAX_seeds, jobID, jobs, w_i, w_j, windowSize,
                               GLOBAL_READ_WORDS, GLOBAL_SKIPPED_WORDS, GLOBAL_CURRENT_WORD, is_j_above_diag);
        } else {
          GLOBAL_SKIPPED_WORDS++;
        }
      } else {
        // tot_pairs += it->second.size() * (it->second.size() - 1) / 2;
        for (auto i = 0ul; i < it->second.size(); i++) {
          for (auto ii = i + 1ul; ii < it->second.size(); ii++) {

            unsigned int ind_i = std::max(it->second[i], it->second[ii]);
            unsigned int ind_j = std::min(it->second[i], it->second[ii]);

            // for the last job only
            if (jobID == jobs) {
              if (all_ind[ind_i].getIdNum() >= (w_i - 1) * windowSize &&
                  all_ind[ind_j].getIdNum() >= (w_j - 1) * windowSize) {
                if (all_ind[ind_j].getIdNum() <
                    (w_j - 1) * windowSize + (all_ind[ind_i].getIdNum() - (w_i - 1) * windowSize)) {
                  e->extendPair(ind_j, ind_i, w, GLOBAL_CURRENT_WORD);
                  tot_pairs++;
                }
              }
            }

            // for all other jobs
            else if ((all_ind[ind_i].getIdNum() >= (w_i - 1) * windowSize &&
                      all_ind[ind_i].getIdNum() < w_i * windowSize) &&
                     (all_ind[ind_j].getIdNum() >= (w_j - 1) * windowSize &&
                      all_ind[ind_j].getIdNum() < w_j * windowSize)) {
              if (is_j_above_diag && all_ind[ind_j].getIdNum() < (w_j - 1) * windowSize + (all_ind[ind_i].getIdNum() -
                                                                                           (w_i - 1) * windowSize)) {
                e->extendPair(ind_j, ind_i, w, GLOBAL_CURRENT_WORD);
                tot_pairs++;
              } else if (!is_j_above_diag &&
                         all_ind[ind_j].getIdNum() >=
                             (w_j - 1) * windowSize + (all_ind[ind_i].getIdNum() - (w_i - 1) * windowSize)) {
                e->extendPair(ind_j, ind_i, w, GLOBAL_CURRENT_WORD);
                tot_pairs++;
              }
            }
          }
        }
      }
    }
    return tot_pairs;
  }
};

#endif // ASMC_HASHING_SEEDHASH_HPP
