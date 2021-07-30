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

#ifndef ASMC_HMM
#define ASMC_HMM

#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "Individual.hpp"
#include "Types.hpp"

#include <Eigen/Core>

#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

#include <zlib.h>

// individual ids and XOR/AND of genotypes
struct PairObservations {
  uint_least8_t iHap;
  uint_least8_t jHap;

  unsigned int iInd;
  unsigned int jInd;

  std::vector<bool> obsBits;
  std::vector<bool> homMinorBits;
};

// individual ids and XOR of genotypes
struct DecodingReturnValues {

  /// output for sum over all pairs
  Eigen::ArrayXXf sumOverPairs;

  /// output for sum over all pairs with genotype 00
  Eigen::ArrayXXf sumOverPairs00;

  /// output for sum over all pairs with genotype 01 or 10
  Eigen::ArrayXXf sumOverPairs01;

  /// output for sum over all pairs with genotype 11
  Eigen::ArrayXXf sumOverPairs11;

  int sites = 0;
  unsigned int states = 0;
  std::vector<bool> siteWasFlippedDuringFolding = {};
};

// does the linear-time decoding
class HMM
{

  int m_batchSize = 64;

  Eigen::ArrayXf m_alphaBuffer;
  Eigen::ArrayXf m_betaBuffer;
  Eigen::ArrayXf m_scalingBuffer;
  Eigen::ArrayXf m_allZeros;

  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> meanPost;
  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> MAP;
  Eigen::ArrayXf currentMAPValue;

  // for decoding
  Data data;

  /**
   * The decoding quantities object, that is owned by Data. This member can be accessed publicly by const ref using
   * Data::getDecodingQuantities
   */
  DecodingQuantities m_decodingQuant;
  DecodingParams decodingParams;

  std::string outFileRoot;
  std::string expectedCoalTimesFile;

  long int sequenceLength;
  uint stateThreshold;
  uint ageThreshold;
  int states;

  int scalingSkip;

  std::vector<float> expectedCoalTimes;
  std::vector<bool> useCSFSatThisPosition;

  std::vector<std::vector<float>> emission1AtSite;
  std::vector<std::vector<float>> emission0minus1AtSite;
  std::vector<std::vector<float>> emission2minus0AtSite;

  /** Calculated by the data module in the HMM constructor */
  std::vector<std::vector<int>> undistinguishedCounts;

  bool noBatches;
  uint64 currPair = 0;

  double timeASMC = 0.0;

  // New members copied from fastSMC
  std::vector<PairObservations> batchObservations;
  unsigned int startBatch;
  unsigned int endBatch;
  std::vector<unsigned int> fromBatch;
  std::vector<unsigned int> toBatch;
  unsigned long int cpt = 0;
  unsigned long int nbSegmentsDetected = 0;
  unsigned long int nbBatch;
  float probabilityThreshold;

  const int precision = 2;
  const float minGenetic = 1e-10f;

  std::vector<PairObservations> m_observationsBatch;

  // output
  DecodingReturnValues m_decodingReturnValues;
  FileUtils::AutoGzOfstream foutPosteriorMeanPerPair;
  FileUtils::AutoGzOfstream foutMAPPerPair;

  // output when decoding specific pairs
  DecodePairsReturnStruct m_decodePairsReturnStruct;

  // output options
  bool m_storePerPairPosteriorMean = false;
  bool m_writePerPairPosteriorMean = false;
  bool m_calculatePerPairPosteriorMean = false;

  bool m_storePerPairMAP = false;
  bool m_writePerPairMAP = false;
  bool m_calculatePerPairMAP = false;

  bool m_storePerPairPosterior = false;
  bool m_storeSumOfPosterior = false;

  Eigen::IOFormat m_eigenOutputFormat{Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n"};

  gzFile gzoutIBD;

  // timing
  std::chrono::duration<double> t1sum = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> t2sum = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> ticksForward = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> ticksBackward = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> ticksCombine = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> ticksSumOverPairs = std::chrono::high_resolution_clock::duration::zero();
  std::chrono::duration<double> ticksOutputPerPair = std::chrono::high_resolution_clock::duration::zero();

  /// Update the size of output structures for per pair posterior mean / MAP
  void updateOutputStructures();

public:
  // constructor
  HMM(Data _data, const DecodingParams& _decodingParams, int _scalingSkip = 1);

  /// Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x
  /// states).
  void decodeAll(int jobs, int jobInd);

  /// The FastSMC version that takes a from and a to
  std::vector<std::vector<float>> decode(const PairObservations& observations, unsigned from, unsigned to);

  /// Original decode version, which calls the three-parameter version with sensible defaults
  std::vector<std::vector<float>> decode(const PairObservations& observations);

  std::pair<std::vector<float>, std::vector<float>> decodeSummarize(const PairObservations& observations);

  PairObservations makePairObs(int_least8_t iHap, unsigned int ind1, int_least8_t jHap, unsigned int ind2);

  /// decode a single pair
  ///
  /// i and j must be a valid index in `individuals`
  /// if noBatches is not set then the pair is saved and processing is delayed until the
  /// observationBatch array is full
  ///
  /// @param i index of first individual
  /// @param j index of second individual
  ///
  void decodePair(const uint i, const uint j);

  /// decode a single hap pair
  ///
  /// i and j must be a valid index in 2 * num_individuals
  /// if noBatches is not set then the pair is saved and processing is delayed until the
  /// observationBatch array is full
  ///
  /// @param i hap index of first hap
  /// @param j hap index of second hap
  ///
  void decodeHapPair(unsigned long i, unsigned long j);

  /// decode a list of pairs
  ///
  /// the input pairs are described using two vectors of indices of `individuals`.
  /// if noBatches is not set then the pairs are processed in batches for efficiency.
  /// This means that after the call to `decodePairs` there might be unprocessed pairs
  /// waiting for the buffer to be full. Use `finishDecoding` to process these.
  ///
  /// @param individualsA vector of indices of first individual
  /// @param individualsB vector of indices of second individual
  ///
  void decodePairs(const std::vector<uint>& individualsA, const std::vector<uint>& individualsB);

  /// decode a list of pairs index by their haplotype indices
  ///
  /// the input pairs are described using two vectors of indices of `haplotypes`.
  /// if noBatches is not set then the pairs are processed in batches for efficiency.
  /// This means that after the call to `decodePairs` there might be unprocessed pairs
  /// waiting for the buffer to be full. Use `finishDecoding` to process these.
  ///
  /// @param individualsA vector of indices of first hap
  /// @param individualsB vector of indices of second hap
  ///
  void decodeHapPairs(const std::vector<unsigned long>& individualsA, const std::vector<unsigned long>& individualsB);

  /// decode a single pair over a segment of the genome
  ///
  /// i and j must be a valid index in `individuals`, `fromPosition` and `toPosition`
  /// must be less than the size of ?
  /// if noBatches is not set then the pair is saved and processing is delayed until the
  /// observationBatch array is full. This is only efficient if subseqent pairs overlap
  /// in the range `fromPosition` -> `toPosition`
  ///
  /// @param i index of first individual
  /// @param j index of second individual
  ///
  void decodeFromHashing(uint i, uint j, uint fromPosition, uint toPosition);

  /// convert generation threshold into state threshold
  uint getStateThreshold();

  /// returns the current buffer of pair observations
  const std::vector<PairObservations>& getBatchBuffer()
  {
    return m_observationsBatch;
  }

  /// returns the decoding quantities calculated thus far
  const DecodingReturnValues& getDecodingReturnValues()
  {
    return m_decodingReturnValues;
  }

  /// get the struct holding return values for each pair decoded using decodePairs()
  DecodePairsReturnStruct& getDecodePairsReturnStruct();

  /// finish decoding pairs
  ///
  /// tells HMM object to finish processing whatever pairs are stored in the
  /// observationsBatch buffer and close output files
  void finishDecoding();

  /**
   * Close the gzipped file stream for writing IBD data
   */
  void closeIBDFile();

  void finishFromHashing();

  /**
   * @return const ref to the decoding quantities owned by this object
   */
  const DecodingQuantities& getDecodingQuantities() const;

  /// Set to true to store per pair posterior mean
  void setStorePerPairPosteriorMean(bool storePerPairPosteriorMean = true);

  /// Set to true to write per pair posterior mean to file
  void setWritePerPairPosteriorMean(bool writePerPairPosteriorMean = true);

  /// Set to true to store per pair MAP
  void setStorePerPairMap(bool storePerPairMAP = true);

  /// Set to true to write per pair MAP to file
  void setWritePerPairMap(bool writePerPairMAP = true);

  /// Set to true to store per pair posterior
  void setStorePerPairPosterior(bool storePerPairPosterior = true);

  /// Set to true to store the sum of posteriors
  void setStoreSumOfPosterior(bool storeSumOfPosterior = true);

private:
  void writeBinaryInfoIntoFile();

  void makeBits(PairObservations& obs, unsigned from, unsigned to);

  /// resets the internal state of HMM to a clean state
  void resetDecoding();

  void prepareEmissions();

  // add pair to batch and run if we have enough
  void addToBatch(std::vector<PairObservations>& obsBatch, const PairObservations& observations);

  // complete with leftover pairs
  void runLastBatch(std::vector<PairObservations>& obsBatch);

  // decode a batch
  void decodeBatch(const std::vector<PairObservations>& obsBatch, unsigned from, unsigned to);

  // forward step
  void forwardBatch(Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                    int curBatchSize, unsigned from, unsigned to);

  // compute next alpha vector in linear time
  void getNextAlphaBatched(float recDistFromPrevious, Eigen::Ref<Eigen::ArrayXf> alphaC, int curBatchSize,
                           Eigen::Ref<Eigen::ArrayXf> previousAlpha, uint pos,
                           Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                           Eigen::Ref<Eigen::ArrayXf> AU, Eigen::Ref<Eigen::ArrayXf> nextAlpha,
                           const std::vector<float>& emission1AtSite, const std::vector<float>& emission0minus1AtSite,
                           const std::vector<float>& emission2minus0AtSite);
  // backward step
  void backwardBatch(Eigen::ArrayXf obsIsZeroBatch, Eigen::ArrayXf obsIsTwoBatch, int curBatchSize, unsigned from,
                     unsigned to);

  // compute previous beta vector in linear time
  void getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, Eigen::Ref<Eigen::ArrayXf> lastComputedBeta,
                              int pos, Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch,
                              Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch, Eigen::Ref<Eigen::ArrayXf> vec,
                              Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> BL,
                              Eigen::Ref<Eigen::ArrayXf> currentBeta, const std::vector<float>& emission1AtSite,
                              const std::vector<float>& emission0minus1AtSite,
                              const std::vector<float>& emission2minus0AtSite);

  // --posteriorSums
  void augmentSumOverPairs(std::vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize);

  float getMAP(std::vector<float> posterior);

  float getPosteriorMean(const std::vector<float>& posterior);

  void writePairIBD(const PairObservations& obs, unsigned int posStart, unsigned int posEnd, float prob,
                    std::vector<float>& posterior, int v, int paddedBatchSize);

  // will eventually write binary output instead of gzipped
  void writePerPairOutput(int actualBatchSize, int paddedBatchSize, const std::vector<PairObservations>& obsBatch);
  void writePerPairOutputFastSMC(int actualBatchSize, int paddedBatchSize,
                                 const std::vector<PairObservations>& obsBatch);

  // *********************************************************************
  // non-batched computations (for debugging and pedagogical reasons only)
  // *********************************************************************

  std::vector<float> getEmission(int pos, int distinguished, int undistinguished, int emissionIndex);

  // forward step
  std::vector<std::vector<float>> forward(const PairObservations& observations, unsigned from, unsigned to);

  void getNextAlpha(float recDistFromPrevious, std::vector<float>& alphaC, std::vector<float>& previousAlpha,
                    std::vector<float>& nextAlpha, std::vector<float>& emission1AtSite,
                    std::vector<float>& emission0minus1AtSite, std::vector<float>& emission2minus0AtSite,
                    float obsIsZero, float obsIsHomMinor);

  // backward step
  std::vector<std::vector<float>> backward(const PairObservations& observations, unsigned from, unsigned to);

  void getPreviousBeta(float recDistFromPrevious, std::vector<float>& lastComputedBeta, std::vector<float>& BL,
                       std::vector<float>& BU, std::vector<float>& currentBeta, std::vector<float>& emission1AtSite,
                       std::vector<float>& emission0minus1AtSite, std::vector<float>& emission2minus0AtSite,
                       float obsIsZero, float obsIsHomMinor);

};

#endif // ASMC_HMM
