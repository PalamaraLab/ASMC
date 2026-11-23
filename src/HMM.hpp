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
#include <cstddef>
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
  float mCmBurnIn = 0.5f;
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

  /// Non-batched forward pass over sites [from, to), used as a reference for the batched code. Returns alpha indexed
  /// [state][site] over all sites, zero outside [from, to). Each site's alpha is rescaled to sum to 1 when the site
  /// index is a multiple of the constructor's `_scalingSkip`, and at `from`.
  std::vector<std::vector<float>> forward(const PairObservations& observations, unsigned from, unsigned to);

  /// Non-batched backward pass over sites [from, to), used as a reference for the batched code. Returns beta indexed
  /// [state][site] over all sites, zero outside [from, to). Rescaled as in `forward`, and at `to - 1`.
  std::vector<std::vector<float>> backward(const PairObservations& observations, unsigned from, unsigned to);

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
  /// @param from the start index (defaults to zero)
  /// @param to the end index (defaults to the sequence length)
  /// @param cmBurnIn the target distance in centimorgans before `from` and after `to`
  ///
  void decodeHapPairs(const std::vector<unsigned long>& individualsA, const std::vector<unsigned long>& individualsB,
                      unsigned from = 0u, unsigned to = 0u, float cmBurnIn = 0.5f);

  /// Memory-efficient parallel pair decoding with forward checkpointing.
  ///
  /// Uses O(sqrt(T) * S * batchSize) memory per thread instead of O(T * S * batchSize) by:
  ///   1. Forward pass: save alpha at sqrt(T) checkpoint intervals (rolling buffers)
  ///   2. Backward pass: process segments right-to-left, re-compute alpha from checkpoint,
  ///      fuse posterior gamma = alpha * beta inline, extract mean/MAP immediately
  ///
  /// This approach also improves cache locality since the working set (~sqrt(T) * S * bs floats)
  /// fits in L2/L3 cache, unlike the full-table approach which strides through multi-GB buffers.
  ///
  /// Supports posterior mean and MAP output only. For features requiring full posterior tables
  /// (sum-over-pairs, per-pair posterior storage, file output), use the sequential path.
  /// Throws std::invalid_argument if per-pair posteriors, the sum of posteriors or per-pair file output are
  /// requested, or in FastSMC mode.
  ///
  /// Batches are decoded on nThreads OpenMP threads (one if nThreads < 1, or without OpenMP). The posterior means and
  /// MAPs are bit-identical to decodeHapPairs's for any number of threads.
  void decodeHapPairsParallel(const std::vector<unsigned long>& individualsA,
                              const std::vector<unsigned long>& individualsB,
                              unsigned from, unsigned to, float cmBurnIn, int nThreads);

  /// Whether decodeHapPairsParallel can produce every per-pair output currently requested: true unless per-pair
  /// posteriors, the sum of posteriors or per-pair file output are requested, or in FastSMC mode.
  bool canDecodeCheckpointed() const;

  /// Number of sites per segment in decodeHapPairsParallel's checkpointing, for a window of `numSites` sites
  /// including burn-in. About sqrt(numSites), which minimises the checkpoints plus one segment of alpha. At least 256,
  /// so that short windows are not split into many segments.
  static unsigned checkpointSegmentLength(unsigned numSites);

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

  /**
   * Open the gzipped file stream for writing IBD data
   */
  void openIBDFile();

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

  bool getStorePerPairPosteriorMean() const;
  bool getWritePerPairPosteriorMean() const;
  bool getStorePerPairMap() const;
  bool getWritePerPairMap() const;
  bool getStorePerPairPosterior() const;
  bool getStoreSumOfPosterior() const;

  void writePairIBD(const PairObservations& obs, unsigned int posStart, unsigned int posEnd, float prob = 0.f,
                    const std::vector<float>& posterior = {});

  const DecodingParams& getDecodingParams() const;

  /// The data being decoded, including any genotypes loaded with loadSubsetGenotypes
  const Data& getData() const
  {
    return data;
  }

  /// Load subset genotypes into the Data object (for streaming pair decoding)
  void loadSubsetGenotypes(const int8_t* haplotypes, int numSubsetHaps);

  /// Clear loaded genotypes to free memory between batches
  void clearGenotypes();

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

  /// Observation indicators for one batch of pairs, as filled by fillBatchObservations in HMM.cpp. Each array is laid
  /// out [site * batchSize + pair] over the whole chromosome. The step functions take non-const references, but none
  /// of these arrays is modified.
  struct BatchObservations {
    Eigen::Ref<Eigen::ArrayXf> isZero;   // 1 where neither haplotype carries the minor allele, else 0
    Eigen::Ref<Eigen::ArrayXf> isTwo;    // 1 where both haplotypes carry it, else 0
    Eigen::Ref<Eigen::ArrayXf> allZeros; // zeros: the observations in the homozygous step of sequence mode
    int batchSize;                       // pairs in the batch, a multiple of the SIMD lane count
  };

  /// Scratch vectors for one batch, used by the batched step functions. All are laid out [state * batchSize + pair]
  /// or [pair].
  struct BatchStepScratch {
    BatchStepScratch(int numStates, int batchSize);
    Eigen::ArrayXf alphaC;       // [state][pair]: suffix sums of the previous alpha over states
    Eigen::ArrayXf AU;           // [pair]: running "above" term of the forward step
    Eigen::ArrayXf vec;          // [state][pair]: next beta times the emission
    Eigen::ArrayXf BU;           // [state][pair]: "above" terms of the backward step
    Eigen::ArrayXf BL;           // [pair]: "below" term of the backward step
    Eigen::ArrayXf sums;         // [pair]: sums over states, for rescaling
    Eigen::ArrayXf scaling;      // [pair]: rescaling factors
    Eigen::ArrayXf betweenSteps; // [state][pair]: alpha or beta between the two steps of a sequence-mode site
  };

  /// Set alpha at site `from` to the initial state probabilities times the emission there, rescaled to sum to 1.
  void initialAlphaBatched(const BatchObservations& obs, unsigned from, Eigen::Ref<Eigen::ArrayXf> alpha,
                           BatchStepScratch& scratch);

  /// Compute alpha at site `pos` from alpha at `pos - 1`, rescaling it if `pos` is a multiple of scalingSkip. In
  /// sequence mode this takes two steps: across the homozygous gap before the marker, then the marker itself.
  void forwardStepBatched(long int pos, Eigen::Ref<Eigen::ArrayXf> previousAlpha, Eigen::Ref<Eigen::ArrayXf> nextAlpha,
                          const BatchObservations& obs, BatchStepScratch& scratch);

  /// Set beta at the last site of a window to 1, rescaled to sum to 1.
  void initialBetaBatched(Eigen::Ref<Eigen::ArrayXf> beta, BatchStepScratch& scratch);

  /// Compute beta at site `pos` from beta at `pos + 1`, rescaling it if `pos` is a multiple of scalingSkip. In
  /// sequence mode this undoes forwardStepBatched's two steps from `pos` to `pos + 1` in reverse order.
  void backwardStepBatched(long int pos, Eigen::Ref<Eigen::ArrayXf> nextBeta, Eigen::Ref<Eigen::ArrayXf> beta,
                           const BatchObservations& obs, BatchStepScratch& scratch);

  /// Full forward pass over sites [from, to), storing every site's alpha in alphaBuffer, laid out
  /// [(site * states + state) * batchSize + pair].
  void forwardPassInto(Eigen::Ref<Eigen::ArrayXf> alphaBuffer, const BatchObservations& obs, unsigned from,
                       unsigned to);

  /// Full backward pass over sites [from, to), storing every site's beta in betaBuffer, laid out as in forwardPassInto.
  void backwardPassInto(Eigen::Ref<Eigen::ArrayXf> betaBuffer, const BatchObservations& obs, unsigned from,
                        unsigned to);

  /// The sites decoded for a batch, including burn-in, and the sites whose results are kept.
  struct DecodingWindow {
    unsigned decodeFrom; // first site decoded
    unsigned decodeTo;   // one past the last site decoded
    unsigned keepFrom;   // first site whose results are kept
    unsigned keepTo;     // one past the last site whose results are kept
  };

  /// A batch of pairs for decodeHapPairsParallel: pairs [firstPair, firstPair + numPairs) of the request, padded to a
  /// multiple of the SIMD lane count by repeating the last pair.
  struct PairBatch {
    std::vector<PairObservations> pairs;
    std::size_t firstPair = 0;
    std::size_t numPairs = 0;
  };

  /// Posterior means and MAP states for one batch, indexed (pair in batch, site - keepFrom). Only the arrays for the
  /// requested outputs are sized.
  struct BatchSummaries {
    Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> posteriorMeans;
    Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> maps;
  };

  /// Buffers for decoding batches with forward checkpointing, owned by one thread and reused across its batches.
  struct CheckpointedDecodeBuffers {
    CheckpointedDecodeBuffers(long int numSites, int numStates, int batchSize, unsigned numSegments,
                              unsigned segmentLength);
    Eigen::ArrayXf obsIsZero;       // [site][pair], see BatchObservations
    Eigen::ArrayXf obsIsTwo;        // [site][pair]
    Eigen::ArrayXf allZeros;        // [site][pair]
    Eigen::ArrayXf checkpointAlpha; // [segment][state][pair]: alpha at the first site of each segment
    Eigen::ArrayXf segmentAlpha;    // [site in segment][state][pair]: alpha recomputed for the current segment
    Eigen::ArrayXf beta;            // [state][pair]: beta at the current site
    Eigen::ArrayXf nextBeta;        // [state][pair]: beta at the site after it
    Eigen::ArrayXf posteriorSums;   // [pair]: scratch for normalising the posterior at one site
    Eigen::ArrayXf mapValue;        // [pair]: largest posterior so far at one site
  };

  /// Make the observations for the pairs (hapsA[i], hapsB[i]) and group them into padded batches of m_batchSize.
  std::vector<PairBatch> makePairBatches(const std::vector<unsigned long>& hapsA,
                                         const std::vector<unsigned long>& hapsB, const DecodingWindow& window);

  /// Decode one batch with forward checkpointing and store its posterior means and MAPs. Alpha is kept only at the
  /// first site of each segment and recomputed one segment at a time, while beta is computed from right to left.
  /// Uses O((numSegments + segmentLength) * states * batchSize) memory instead of O(sites * states * batchSize). The
  /// results are bit-identical to decodeBatch followed by writePerPairOutput.
  void decodeBatchCheckpointed(const PairBatch& batch, const DecodingWindow& window, unsigned segmentLength,
                               CheckpointedDecodeBuffers& buffers, BatchSummaries& summaries);

  /// Forward pass of decodeBatchCheckpointed: store alpha at the first site of every segment in checkpointAlpha.
  void forwardPassCheckpointed(const BatchObservations& obs, const DecodingWindow& window, unsigned segmentLength,
                               CheckpointedDecodeBuffers& buffers, BatchStepScratch& scratch);

  /// Recompute alpha at sites [segmentStart, segmentEnd) into segmentAlpha, from the checkpoint at segmentStart.
  void recomputeSegmentAlpha(const BatchObservations& obs, unsigned segment, unsigned segmentStart, unsigned segmentEnd,
                             CheckpointedDecodeBuffers& buffers, BatchStepScratch& scratch);

  /// Turn alpha at one site into the posterior (alpha * beta, normalised over states), in place, and add it to the
  /// posterior means and MAPs of the batch's real (unpadded) pairs at column `offset`.
  void summarisePosteriorAtSite(Eigen::Ref<Eigen::ArrayXf> alpha, int batchSize, std::size_t numPairs, long int offset,
                                CheckpointedDecodeBuffers& buffers, BatchSummaries& summaries);

  /// Copy each batch's summaries and pair indices into m_decodePairsReturnStruct, in the order the pairs were given.
  void storeBatchSummaries(const std::vector<PairBatch>& batches, const std::vector<BatchSummaries>& summaries);

  // Sequential wrappers using class-level buffers
  void forwardBatch(Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                    int curBatchSize, unsigned from, unsigned to);
  void backwardBatch(Eigen::ArrayXf obsIsZeroBatch, Eigen::ArrayXf obsIsTwoBatch, int curBatchSize, unsigned from,
                     unsigned to);

  // compute next alpha vector in linear time
  void getNextAlphaBatched(float recDistFromPrevious, Eigen::Ref<Eigen::ArrayXf> alphaC, int curBatchSize,
                           Eigen::Ref<Eigen::ArrayXf> previousAlpha, uint pos,
                           Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                           Eigen::Ref<Eigen::ArrayXf> AU, Eigen::Ref<Eigen::ArrayXf> nextAlpha,
                           const std::vector<float>& emission1AtSite, const std::vector<float>& emission0minus1AtSite,
                           const std::vector<float>& emission2minus0AtSite);

  // compute previous beta vector in linear time
  void getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, Eigen::Ref<Eigen::ArrayXf> lastComputedBeta,
                              int pos, Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch,
                              Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch, Eigen::Ref<Eigen::ArrayXf> vec,
                              Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> BL,
                              Eigen::Ref<Eigen::ArrayXf> currentBeta, const std::vector<float>& emission1AtSite,
                              const std::vector<float>& emission0minus1AtSite,
                              const std::vector<float>& emission2minus0AtSite);

  // --posteriorSums
  void augmentSumOverPairs(std::vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize, unsigned from, unsigned to);

  float getMAP(std::vector<float> posterior);

  float getPosteriorMean(const std::vector<float>& posterior);

  // will eventually write binary output instead of gzipped
  void writePerPairOutput(int actualBatchSize, int paddedBatchSize, const std::vector<PairObservations>& obsBatch);
  void writePerPairOutputFastSMC(int actualBatchSize, int paddedBatchSize,
                                 const std::vector<PairObservations>& obsBatch);

  // *********************************************************************
  // non-batched computations (for debugging and pedagogical reasons only)
  // *********************************************************************

  std::vector<float> getEmission(int pos, int distinguished, int undistinguished, int emissionIndex);

  void getNextAlpha(float recDistFromPrevious, std::vector<float>& alphaC, std::vector<float>& previousAlpha,
                    std::vector<float>& nextAlpha, std::vector<float>& emission1AtSite,
                    std::vector<float>& emission0minus1AtSite, std::vector<float>& emission2minus0AtSite,
                    float obsIsZero, float obsIsHomMinor);

  void getPreviousBeta(float recDistFromPrevious, std::vector<float>& lastComputedBeta, std::vector<float>& BL,
                       std::vector<float>& BU, std::vector<float>& currentBeta, std::vector<float>& emission1AtSite,
                       std::vector<float>& emission0minus1AtSite, std::vector<float>& emission2minus0AtSite,
                       float obsIsZero, float obsIsHomMinor);

};

#endif // ASMC_HMM
