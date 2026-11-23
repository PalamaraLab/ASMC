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

#include <cstdint>
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

  /**
   * Construct from in-memory arrays, bypassing all file I/O.
   *
   * @param numSites number of variant sites
   * @param numDiploidSamples number of diploid individuals
   * @param physPositions physical positions (bp), length numSites
   * @param genPositions genetic positions (Morgans), length numSites
   * @param haplotypes flat row-major genotype matrix (0/1), length numSites * 2*numDiploidSamples
   * @param usingCSFS whether to use CSFS
   * @param foldData whether to fold to minor allele
   * @throws std::invalid_argument if a genotype is not 0 or 1
   */
  Data(int numSites, unsigned long numDiploidSamples,
       std::vector<int> physPositions,
       std::vector<float> genPositions,
       const int8_t* haplotypes,
       bool usingCSFS, bool foldData);

  /**
   * Construct from precomputed allele counts, without loading genotypes.
   * Genotypes can be loaded later via loadSubsetGenotypes(). With foldData, sites where the derived allele is the
   * majority are folded, as in the constructor from genotypes.
   * @throws std::invalid_argument if the count vectors differ in length, or a derived count is outside [0, total]
   */
  Data(std::vector<int> derivedCounts, std::vector<int> totalCounts,
       std::vector<int> physPositions, std::vector<float> genPositions,
       bool usingCSFS, bool foldData);

  /**
   * Load genotypes for a subset of haplotypes (for streaming pair decoding), as described for setGenotypes.
   * Does NOT recompute emissions — allele counts must already be set.
   * @throws std::invalid_argument if a genotype is not 0 or 1
   */
  void loadSubsetGenotypes(const int8_t* haplotypes, int numSubsetHaps);

  /**
   * Clear loaded genotypes to free memory between batches.
   */
  void clearGenotypes();

  /// Number of sites in the haplotypes file for `inFileRoot`, i.e. its non-blank lines.
  static int countHapLines(std::string inFileRoot);

  /// Number of samples in the samples file for `inFileRoot` (see isSampleLine).
  static int countSamplesLines(std::string inFileRoot);

  /**
   * Path of the samples file for `inFileRoot`: `.samples`, or else `.sample`.
   * @throws std::runtime_error if neither exists
   */
  static std::string findSamplesFile(const std::string& inFileRoot);

  /**
   * Path of the haplotypes file for `inFileRoot`: the first of `.hap.gz`, `.hap`, `.haps.gz` and `.haps` that exists.
   * @throws std::runtime_error if none exists
   */
  static std::string findHapsFile(const std::string& inFileRoot);

  /**
   * Path of the genetic map: `mapFile` if it is not empty, otherwise `.map.gz` or `.map` next to `inFileRoot`.
   * @throws std::runtime_error if the chosen file does not exist
   */
  static std::string findMapFile(const std::string& inFileRoot, const std::string& mapFile);

  /**
   * Whether a line of a samples file, split into whitespace-separated fields, names a sample. Blank lines and the two
   * header lines ("ID_1 ID_2 missing" and "0 0 0") do not.
   * @throws std::runtime_error if any other line has fewer than two fields (family ID and individual ID)
   */
  static bool isSampleLine(const std::vector<std::string>& fields);

  /**
   * Chromosome number from the first column of a haplotypes file, e.g. 1 from "1", "chr1" or "1:12345". Returns 0 if
   * the name is not a number between 1 and 1260 (the most chromosomes of any known organism), e.g. "X".
   */
  static int parseChromosomeNumber(const std::string& chromosome);

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

  /**
   * Check that readHaps parsed a record for every site counted by countHapLines.
   * @throws std::runtime_error naming the first site that could not be read
   */
  void requireAllSitesRead(unsigned long sitesRead) const;

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

  /// Set recRateAtMarker from geneticPositions and physicalPositions: the genetic distance per base pair from the
  /// previous site, or 0 where the genetic distance is 0. Site 0 takes site 1's rate.
  void computeRecombinationRates();

  /**
   * Replace the individuals with `numHaps / 2` with placeholder names and genotypes from `haplotypes`, laid out
   * [site * numHaps + hap] with 1 for the derived allele. Sites in siteWasFlippedDuringFolding are flipped, so that 1
   * then means the minor allele, as it does for the allele counts.
   * @throws std::invalid_argument if a genotype is not 0 or 1
   */
  void setGenotypes(const int8_t* haplotypes, int numHaps);

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

  /**
   * Draw `sampleSize` items without replacement from `populationSize` items, `numberOfSuccesses` of which are
   * successes, and return how many successes were drawn. Returns -1 if the draw is impossible: more successes or more
   * draws than items. When the CSFS is not used, the draw can be larger than the population.
   */
  static int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);


  void readGeneticMap(unsigned long int bp, std::vector<std::pair<unsigned long int, double>>& genetic_map,
                      unsigned int& cur_g, unsigned int pos);

  void addMarker(unsigned long int physicalPosition, double geneticPosition, unsigned int pos);




};

#endif // ASMC_DATA_HPP
