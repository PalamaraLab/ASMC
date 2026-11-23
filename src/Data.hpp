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

#ifndef ASMC_DATA_HPP
#define ASMC_DATA_HPP

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Individual.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"

class Data
{

public:

  std::vector<std::string> FamIDList = {};
  std::vector<std::string> IIDList = {};
  std::vector<std::string> famAndIndNameList = {};
  std::vector<Individual> individuals = {};

  unsigned long sampleSize = 0ul;
  unsigned long haploidSampleSize = 0ul;
  int sites = 0;
  bool decodingUsesCSFS = false;

  /// Whether each job loads only the samples in its windows (w_i, w_j). True only for FastSMC with hashing, whose
  /// hashing step compares samples within a job's windows. Everything else loads all samples, and HMM::decodeAll
  /// divides the pairs between jobs.
  bool mJobbing = false;
  bool foldToMinorAlleles = false;
  std::vector<float> geneticPositions = {};
  std::vector<int> physicalPositions = {};
  std::vector<bool> siteWasFlippedDuringFolding = {};
  std::vector<float> recRateAtMarker = {};

  // Variables relating to FastSMC
  int chrNumber = 0;
  unsigned int windowSize = 0u; // window size in triangles for each job
  unsigned int w_i = 0u;        // window id for ind_i for jobs
  unsigned int w_j = 0u;        // window id for ind_j for jobs
  bool is_j_above_diag = false;
  std::unordered_map<int, unsigned int> physicalPositionsMap = {}; // map where key=physicalPosition, value=indexPosition

  /**
   * Construct the data object, which also constructs the decoding quantities that will be owned by this object
   *
   * @param params the decoding params
   */
  explicit Data(const DecodingParams& params);

  static int countHapLines(std::string inFileRoot);
  static int countSamplesLines(std::string inFileRoot);

  /**
   * Calculate the undistinguished counts
   *
   * @param numCsfsSamples the number of CSFS samples
   * @return the undistinguished counts
   */
  std::vector<std::vector<int>> calculateUndistinguishedCounts(int numCsfsSamples) const;

  const std::vector<std::string>& getSnpIDs() const;

  /**
   * Whether a sample is loaded by this job: always true unless mJobbing, in which case the sample must lie in one of
   * the job's windows.
   *
   * @param sampleIndex index of the (diploid) sample in the samples file, counting from 0 after any header
   */
  bool isSampleInJob(unsigned long sampleIndex) const;

private:
  int mJobInd = 1; // this job's index, 1..mJobs
  int mJobs = 1;   // the number of jobs

  /**
   * Read the samples file and populate members `FamIDList`, `IIDList` and `famAndIndNameList` with the samples in this
   * job (see isSampleInJob).
   *
   * @param inFileRoot location of input files
   */
  void readSamplesList(const std::string& inFileRoot);

  /// Read haplotypes for ASMC, which always loads every sample.
  void readHaps(std::string inFileRoot, bool foldToMinorAlleles);

  /// Read haplotypes and the genetic map for FastSMC, loading only the samples in this job (see isSampleInJob).
  void readHaps(const std::string& inFileRoot, bool foldToMinorAlleles,
                std::vector<std::pair<unsigned long int, double>>& genetic_map);

  /**
   * Read Plink-format map file
   * @param inFileRoot
   * @param mapFile: optional direct path to map file, used if map file is not in inFileRoot
   */
  void readMap(const std::string& inFileRoot, const std::string& mapFile = "");

  /**
   * Subsumed functionality from FastSMC to read genetic map as a vector of pairs.
   * TODO: can this be harmonised with the other readMap method?
   * @param inFileRoot
   * @param mapFile: optional direct path to map file, used if map file is not in inFileRoot
   * @return
   */
  static std::vector<std::pair<unsigned long, double>> readMapFastSMC(const std::string& inFileRoot,
                                                                      const std::string& mapFile = "");

  std::vector<int> totalSamplesCount;
  std::vector<int> derivedAlleleCounts;
  std::vector<std::string> SNP_IDs;

  static int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);


  void readGeneticMap(unsigned long int bp, std::vector<std::pair<unsigned long int, double>>& genetic_map,
                      unsigned int& cur_g, unsigned int pos);

  void addMarker(unsigned long int physicalPosition, double geneticPosition, unsigned int pos);




};

#endif // ASMC_DATA_HPP
