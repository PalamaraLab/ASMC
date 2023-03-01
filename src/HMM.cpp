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

#include "HMM.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

#include <Eigen/Core>

#include <fmt/format.h>

#include "AvxDefinitions.hpp"
#include "HmmUtils.hpp"
#include "MemoryUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"

using namespace std;

// read expected times from a file
vector<float> readExpectedTimesFromIntervalsFile(const char* fileName)
{
  FileUtils::AutoGzIfstream br;
  br.openOrExit(fileName);
  string line;
  vector<float> expCoalTimes;
  while (getline(br, line)) {
    vector<string> splitString;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitString.push_back(buf);
    if (splitString.size() != 3) {
      cerr << fileName
           << " should have \"intervalStart\texpectedCoalescentTime\tintervalEnd\" at "
              "each line."
           << endl;
      exit(1);
    }
    expCoalTimes.push_back(StringUtils::stof(splitString[1]));
  }
  return expCoalTimes;
}

// constructor
HMM::HMM(Data _data, const DecodingParams& _decodingParams, int _scalingSkip)
    : data(std::move(_data)), m_decodingQuant(_decodingParams.decodingQuantFile), decodingParams(_decodingParams),
      scalingSkip(_scalingSkip), noBatches(_decodingParams.noBatches)
{
  if (decodingParams.hashing && !decodingParams.FastSMC) {
    cerr << "Identification only is not yet supported Cannot have hashing==true and FastSMC==false.\n";
    exit(1);
  }

  cout << "Will decode using " << MODE << " instruction set.\n\n";
  outFileRoot = decodingParams.outFileRoot;
  expectedCoalTimesFile = decodingParams.expectedCoalTimesFile;
  sequenceLength = data.sites;
  states = m_decodingQuant.states;
  useCSFSatThisPosition = vector<bool>(sequenceLength, false);
  emission1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  emission0minus1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  emission2minus0AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  prepareEmissions();

  m_batchSize = decodingParams.batchSize;

  for (int i = 0; i < m_batchSize; i++) {
    fromBatch.push_back(0);
    toBatch.push_back(sequenceLength);
  }

  // get state threshold
  stateThreshold = getStateThreshold();
  // probabilityThreshold = (1./decodingQuant->states)*stateThreshold;

  probabilityThreshold = 0.f;
  for (int i = 0; i < stateThreshold; i++) {
    probabilityThreshold += m_decodingQuant.initialStateProb.at(i);
  }

  if (decodingParams.noConditionalAgeEstimates) {
    ageThreshold = states;
  } else {
    ageThreshold = stateThreshold;
  }

  startBatch = sequenceLength;
  endBatch = 0;

  // allocate buffers
  m_scalingBuffer.resize(m_batchSize);
  m_alphaBuffer.resize(sequenceLength * states * m_batchSize);
  m_betaBuffer.resize(sequenceLength * states * m_batchSize);

  m_allZeros.resize(sequenceLength * m_batchSize);
  m_allZeros.setZero();

  // Output variables
  m_calculatePerPairPosteriorMean = decodingParams.doPerPairPosteriorMean;
  m_calculatePerPairMAP = decodingParams.doPerPairMAP;
  updateOutputStructures();

  // output for python interface (TODO: not sure if this is the right place)
  m_decodingReturnValues.sites = data.sites;
  m_decodingReturnValues.states = m_decodingQuant.states;
  m_decodingReturnValues.siteWasFlippedDuringFolding = data.siteWasFlippedDuringFolding;
}

PairObservations HMM::makePairObs(int_least8_t iHap, unsigned int ind1, int_least8_t jHap, unsigned int ind2)
{
  PairObservations ret;
  ret.iHap = iHap;
  ret.jHap = jHap;
  ret.iInd = ind1;
  ret.jInd = ind2;

  //\todo: ideally all calls to makeBits would be in one place, but hashing calls are in addToBatch and runLastBatch
  // because this is where from and to can be calculated
  const bool makeBitsOnWholeSequence = !(decodingParams.FastSMC && decodingParams.hashing);
  if (makeBitsOnWholeSequence || noBatches) {
    makeBits(ret, 0, sequenceLength);
  }

  return ret;
}

void HMM::makeBits(PairObservations& obs, unsigned from, unsigned to)
{
  unsigned iInd = obs.iInd;
  unsigned jInd = obs.jInd;
  obs.obsBits =
      asmc::subsetXorVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
                         obs.jHap == 1 ? data.individuals[jInd].genotype1 : data.individuals[jInd].genotype2, from, to);
  obs.homMinorBits =
      asmc::subsetAndVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
                         obs.jHap == 1 ? data.individuals[jInd].genotype1 : data.individuals[jInd].genotype2, from, to);
}

void HMM::prepareEmissions()
{
  undistinguishedCounts = data.calculateUndistinguishedCounts(m_decodingQuant.CSFSSamples);

  if (decodingParams.skipCSFSdistance < std::numeric_limits<float>::infinity()) {
    useCSFSatThisPosition[0] = true;
    float lastGenCSFSwasUsed = 0.f;
    for (int pos = 1; pos < sequenceLength; pos++) {
      if (data.geneticPositions[pos] - lastGenCSFSwasUsed >= decodingParams.skipCSFSdistance) {
        // this position is a CSFS position
        useCSFSatThisPosition[pos] = true;
        lastGenCSFSwasUsed = data.geneticPositions[pos];
      }
    }
  }
  for (int pos = 0; pos < sequenceLength; pos++) {
    if (useCSFSatThisPosition[pos]) {
      int undistAtThisSiteFor0dist = undistinguishedCounts[pos][0];
      int undistAtThisSiteFor1dist = undistinguishedCounts[pos][1];
      int undistAtThisSiteFor2dist = undistinguishedCounts[pos][2];
      if (decodingParams.foldData) {
        // working with folded data
        for (int k = 0; k < states; k++) {
          if (undistAtThisSiteFor1dist >= 0) {
            emission1AtSite[pos][k] = decodingParams.decodingSequence
                                          ? m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor1dist][1][k]
                                          : m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
          } else {
            emission1AtSite[pos][k] = 0.f;
          }
          emission0minus1AtSite[pos][k] =
              decodingParams.decodingSequence
                  ? (m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k] - emission1AtSite[pos][k])
                  : (m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k] -
                     emission1AtSite[pos][k]);
          if (undistAtThisSiteFor2dist >= 0) {
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                       m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                    : (m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                       m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
          } else {
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (0 - m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                    : (0 - m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
          }
        }
      } else {
        // working with unfolded data
        for (int k = 0; k < states; k++) {
          if (undistAtThisSiteFor1dist >= 0) {
            emission1AtSite[pos][k] = decodingParams.decodingSequence
                                          ? m_decodingQuant.CSFSmap[undistAtThisSiteFor1dist][1][k]
                                          : m_decodingQuant.ascertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
          } else {
            emission1AtSite[pos][k] = 0.f;
          }
          float emission0AtThisSiteAndState = 0.f;
          if (undistAtThisSiteFor0dist >= 0) {
            emission0AtThisSiteAndState = decodingParams.decodingSequence
                                              ? m_decodingQuant.CSFSmap[undistAtThisSiteFor0dist][0][k]
                                              : m_decodingQuant.ascertainedCSFSmap[undistAtThisSiteFor0dist][0][k];
          }
          emission0minus1AtSite[pos][k] = emission0AtThisSiteAndState - emission1AtSite[pos][k];
          if (undistAtThisSiteFor2dist >= 0) {
            int dist = 2;
            int undist = undistAtThisSiteFor2dist;
            if (undistAtThisSiteFor2dist == m_decodingQuant.CSFSSamples - 2) {
              // for monomorphic derived, fold to CSFS[0][0]
              dist = 0;
              undist = 0;
            }
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (m_decodingQuant.CSFSmap[undist][dist][k] - emission0AtThisSiteAndState)
                    : (m_decodingQuant.ascertainedCSFSmap[undist][dist][k] - emission0AtThisSiteAndState);
          } else {
            emission2minus0AtSite[pos][k] = 0 - emission0AtThisSiteAndState;
          }
        }
      }
    } else {
      // this position is not a CSFS position
      for (int k = 0; k < states; k++) {
        emission1AtSite[pos][k] = decodingParams.decodingSequence ? m_decodingQuant.classicEmissionTable[1][k]
                                                                  : m_decodingQuant.compressedEmissionTable[1][k];
        emission0minus1AtSite[pos][k] =
            decodingParams.decodingSequence
                ? (m_decodingQuant.classicEmissionTable[0][k] - m_decodingQuant.classicEmissionTable[1][k])
                : (m_decodingQuant.compressedEmissionTable[0][k] - m_decodingQuant.compressedEmissionTable[1][k]);
        // emission2 = emission0
        emission2minus0AtSite[pos][k] = 0.f;
      }
    }
  }
}

void HMM::resetDecoding()
{
  const long int seqFrom = static_cast<long int>(*std::min_element(fromBatch.begin(), fromBatch.end()));
  const long int seqTo = static_cast<long int>(*std::max_element(toBatch.begin(), toBatch.end()));
  const unsigned seqLen = seqTo - seqFrom;

  if (m_writePerPairPosteriorMean && !decodingParams.FastSMC) {
    if (foutPosteriorMeanPerPair) {
      foutPosteriorMeanPerPair.close();
    }
    foutPosteriorMeanPerPair.openOrExit(outFileRoot + ".perPairPosteriorMeans.gz");
  }
  if (m_writePerPairMAP && !decodingParams.FastSMC) {
    if (foutMAPPerPair) {
      foutMAPPerPair.close();
    }
    foutMAPPerPair.openOrExit(outFileRoot + ".perPairMAP.gz");
  }

  m_decodingReturnValues.sumOverPairs = Eigen::ArrayXXf::Zero(seqLen, states);

  if (decodingParams.doMajorMinorPosteriorSums) {
    m_decodingReturnValues.sumOverPairs00 = Eigen::ArrayXXf::Zero(seqLen, states);
    m_decodingReturnValues.sumOverPairs01 = Eigen::ArrayXXf::Zero(seqLen, states);
    m_decodingReturnValues.sumOverPairs11 = Eigen::ArrayXXf::Zero(seqLen, states);
  }
}

// Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x states).
void HMM::decodeAll(int jobs, int jobInd)
{

  auto t0 = std::chrono::high_resolution_clock::now();
  Timer timer;

  resetDecoding();

  const vector<Individual>& individuals = data.individuals;
  uint64 lastPercentage = -1;
  uint64 N = individuals.size();
  uint64 pairs = 0, pairsJob = 0;

  if (decodingParams.FastSMC) {
    // create IBD output file
    if (!decodingParams.BIN_OUT) {
      gzoutIBD = gzopen( (decodingParams.outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".FastSMC.ibd.gz").c_str(), "w" );
    } else {
      gzoutIBD = gzopen( (decodingParams.outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".FastSMC.bibd.gz").c_str(), "wb" );
      writeBinaryInfoIntoFile();
    }

    if (decodingParams.hashing) {
      return;
    }
  }

  // calculate total number of pairs to decode
  uint64 totPairs;
  if (!decodingParams.withinOnly) {
    totPairs = 2 * N * N - N;
  } else {
    totPairs = N;
  }

  // figure out the range of pairs for this job number
  uint64 pairsStart = totPairs * (jobInd - 1) / jobs;
  uint64 pairsEnd = totPairs * jobInd / jobs;
  uint64 totPairsJob = pairsEnd - pairsStart;

  fmt::print("jobInd: {}   jobs: {}\n", jobInd, jobs);
  fmt::print("pairsStart: {}   pairsEnd: {}   totPairsJob: {}\n", pairsStart, pairsEnd, totPairsJob);

  // alloc
  m_observationsBatch.clear();
  for (uint i = 0; i < individuals.size(); i++) {
    if (!decodingParams.withinOnly) {
      for (uint j = 0; j < i; j++) {
        // different individuals; decode 2 haps x 2 haps
        for (int iHap = 1; iHap <= 2; iHap++) {
          for (int jHap = 1; jHap <= 2; jHap++) {
            if (pairsStart <= pairs && pairs < pairsEnd) {
              PairObservations observations = makePairObs(jHap, j, iHap, i);
              if (noBatches) {
                decode(observations);
              } else {
                addToBatch(m_observationsBatch, observations);
              }
              pairsJob++;
            }
            pairs++;
            currPair = pairs;
          }
        }
      }
    }
    // this is the same individual; only decode across chromosomes
    if (pairsStart <= pairs && pairs < pairsEnd) {
      PairObservations observations = makePairObs(1, i, 2, i);
      if (noBatches) {
        decode(observations);
      } else {
        addToBatch(m_observationsBatch, observations);
      }
      pairsJob++;
    }
    pairs++;
    currPair = pairs;
    uint64 percentage = 100 * pairsJob / totPairsJob;
    if (percentage != lastPercentage) {
      cout << "\rDecoding progress: " << percentage << "%"
           << "  (" << pairsJob << "/" << totPairsJob << ")" << flush;
    }
    lastPercentage = percentage;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  auto ticksDecodeAll = t1 - t0;
  printf("\nDecoded %" PRIu64 " pairs in %.9f seconds.\n", pairsJob, timer.update_time());

  // print some stats
  asmc::printPctTime("forward", ticksForward / ticksDecodeAll);
  asmc::printPctTime("backward", ticksBackward / ticksDecodeAll);
  asmc::printPctTime("combine", ticksCombine / ticksDecodeAll);
  asmc::printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
  asmc::printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
  asmc::printPctTime("other",
                     1 - (ticksForward + ticksBackward + ticksCombine + ticksSumOverPairs + ticksOutputPerPair) /
                             ticksDecodeAll);

  finishDecoding();
}

void HMM::writeBinaryInfoIntoFile()
{
  const bool outputIbdScore = !decodingParams.hashingOnly;
  gzwrite(gzoutIBD, (char*)&decodingParams.outputIbdSegmentLength, sizeof(bool));
  gzwrite(gzoutIBD, (char*)&outputIbdScore, sizeof(bool));
  gzwrite(gzoutIBD, (char*)&decodingParams.doPerPairPosteriorMean, sizeof(bool));
  gzwrite(gzoutIBD, (char*)&decodingParams.doPerPairMAP, sizeof(bool));
  gzwrite(gzoutIBD, (char*)&data.chrNumber, sizeof(int));
  unsigned int nbInd = data.individuals.size();
  unsigned int lengthFamid;
  unsigned int lengthIid;
  gzwrite(gzoutIBD, (char*)&nbInd, sizeof(unsigned int));
  for (unsigned int i = 0; i < nbInd; i++) {
    lengthFamid = data.FamIDList[i].size();
    gzwrite(gzoutIBD, (char*)&lengthFamid, sizeof(unsigned int));
    gzwrite(gzoutIBD, data.FamIDList[i].c_str(), lengthFamid);
    lengthIid = data.IIDList[i].size();
    gzwrite(gzoutIBD, (char*)&lengthIid, sizeof(unsigned int));
    gzwrite(gzoutIBD, data.IIDList[i].c_str(), lengthIid);
  }
}

void HMM::decodePairs(const vector<uint>& individualsA, const vector<uint>& individualsB)
{
  if (individualsA.size() != individualsB.size()) {
    throw runtime_error("vector of A indicies must be the same size as vector of B indicies");
  }
  for (size_t i = 0; i < individualsA.size(); ++i) {
    decodePair(individualsA[i], individualsB[i]);
  }
}

void HMM::decodePair(const uint i, const uint j)
{
  const vector<Individual>& individuals = data.individuals;
  assert(i < individuals.size());
  assert(j < individuals.size());

  if (i != j) {
    // different individuals; decode 2 haps x 2 haps
    for (int iHap = 1; iHap <= 2; iHap++) {
      for (int jHap = 1; jHap <= 2; jHap++) {
        PairObservations observations = makePairObs(iHap, i, jHap, j);
        if (noBatches) {
          decode(observations);
        } else {
          addToBatch(m_observationsBatch, observations);
        }
      }
    }
  } else {
    // this is the same individual; only decode across chromosomes
    PairObservations observations = makePairObs(1, i, 2, i);
    if (noBatches) {
      decode(observations);
    } else {
      addToBatch(m_observationsBatch, observations);
    }
  }
}

void HMM::decodeHapPair(const unsigned long i, const unsigned long j)
{
  auto numHaps = 2ul * data.individuals.size();
  assert(i < numHaps);
  assert(j < numHaps);

  auto [iInd, iHap] = asmc::hapToDipId(i);
  auto [jInd, jHap] = asmc::hapToDipId(j);

  PairObservations obs = makePairObs(static_cast<int_least8_t>(iHap), iInd, static_cast<int_least8_t>(jHap), jInd);

  if (noBatches) {
    decode(obs);
  } else {
    addToBatch(m_observationsBatch, obs);
  }
}

void HMM::decodeHapPairs(const std::vector<unsigned long>& hapsA, const std::vector<unsigned long>& hapsB,
                         const unsigned from, unsigned to, const float cmBurnIn)
{
  const unsigned sequenceLength = static_cast<unsigned>(data.sites);

  if (to == 0u) {
    to = sequenceLength;
  }

  if (from >= to || to > sequenceLength) {
    throw std::runtime_error(
      fmt::format("Require 0 <= from < to <= sequenceLength but got 0 <= {} < {} <= {}\n", from, to, sequenceLength)
    );
  }

  if (cmBurnIn < 0.f) {
        throw std::runtime_error(fmt::format("Burn-in dist in cM should be >= 0.0 but got cmBirnIn = {}\n", cmBurnIn));
  }

  std::fill(fromBatch.begin(), fromBatch.end(), from);
  std::fill(toBatch.begin(), toBatch.end(), to);
  mCmBurnIn = cmBurnIn;

  if (hapsA.size() != hapsB.size()) {
    throw std::runtime_error("vector of A indices must be the same size as vector of B indices");
  }
  for (size_t i = 0; i < hapsA.size(); ++i) {
    decodeHapPair(hapsA[i], hapsB[i]);
  }
}

void HMM::decodeFromHashing(const uint indivID1, const uint indivID2, const uint fromPosition, const uint toPosition)
{
  const vector<Individual>& individuals = data.individuals;

  // indivID1 & indivID2 are indices of chromosomes corresponding to individuals indivID1/2 and indivID2/2
  assert(indivID1 / 2 < individuals.size());
  assert(indivID2 / 2 < individuals.size());
  assert(fromPosition < sequenceLength);
  assert(toPosition < sequenceLength);

  Timer timerASMC;

  // ID of individual j must be smaller than ID of individual i
  unsigned int jInd = indivID1 / 2;
  unsigned int iInd = indivID2 / 2;

  PairObservations observation = makePairObs(indivID1 % 2 == 0 ? 1 : 2, jInd, indivID2 % 2 == 0 ? 1 : 2, iInd);

  if (noBatches) {
    decode(observation, fromPosition, toPosition);
  } else {
    nbBatch = cpt % static_cast<unsigned long>(m_batchSize);
    fromBatch[nbBatch] = fromPosition;
    toBatch[nbBatch] = toPosition;
    addToBatch(batchObservations, observation);
    cpt++;
  }
  if (cpt % 10000 == 0) {
    cout << "\rnumber of decoded segments: " << cpt << "\t"
         << "\tdetected segments: " << nbSegmentsDetected << flush;
  }
  timeASMC += timerASMC.update_time();
}

uint HMM::getStateThreshold()
{
  uint result = 0u;
  const vector<float>& disc = m_decodingQuant.discretization;

  while (disc[result] < static_cast<float>(decodingParams.time) && result < m_decodingQuant.states) {
    result++;
  }
  return result;
}

void HMM::finishDecoding()
{
  runLastBatch(m_observationsBatch);
  if (m_writePerPairPosteriorMean) {
    foutPosteriorMeanPerPair.close();
  }
  if (m_writePerPairMAP) {
    foutMAPPerPair.close();
  }
}

void HMM::closeIBDFile()
{
  gzclose(gzoutIBD);
}

void HMM::finishFromHashing()
{
  if (!noBatches && !decodingParams.hashingOnly) {
    runLastBatch(batchObservations);
  }

  closeIBDFile();
}

// add pair to batch and run if we have enough
void HMM::addToBatch(vector<PairObservations>& obsBatch, const PairObservations& observations)
{
  obsBatch.push_back(observations);
  if (static_cast<int>(obsBatch.size()) == m_batchSize) {

    // taking the maximum 'to' position and the minimum 'from' position in the batch
    startBatch = *std::min_element(fromBatch.begin(), fromBatch.end());
    endBatch = *std::max_element(toBatch.begin(), toBatch.end());

    unsigned int from = asmc::getFromPosition(data.geneticPositions, startBatch, mCmBurnIn);
    unsigned int to = asmc::getToPosition(data.geneticPositions, endBatch, mCmBurnIn);

    const bool fastsmcWithHashing = decodingParams.FastSMC && decodingParams.hashing;
    const bool fromOrTo = from > 0u || to < sequenceLength;
    const bool makeBitsOnlyOnSubsequence = fastsmcWithHashing || fromOrTo;
    if (makeBitsOnlyOnSubsequence) {
      for (auto& obs : obsBatch) {
        makeBits(obs, from, to);
      }
    }

    // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x m_batchSize]
    decodeBatch(obsBatch, from, to);

    augmentSumOverPairs(obsBatch, m_batchSize, m_batchSize, from, to);
    if ((m_calculatePerPairMAP || m_calculatePerPairPosteriorMean) && !decodingParams.FastSMC) {
      writePerPairOutput(m_batchSize, m_batchSize, obsBatch);
    }

    if (decodingParams.FastSMC) {
      writePerPairOutputFastSMC(m_batchSize, m_batchSize, obsBatch);
    }

    // reinitializing batch variables
    startBatch = sequenceLength;
    endBatch = 0;
    obsBatch.clear();
  }
}

// complete with leftover pairs
void HMM::runLastBatch(vector<PairObservations>& obsBatch)
{
  if (obsBatch.empty()) {
    return;
  }

  auto actualBatchSize = obsBatch.size();

  // taking the maximum To position and the minimum From position in the batch
  startBatch = *std::min_element(fromBatch.begin(), fromBatch.begin() + actualBatchSize);
  endBatch = *std::max_element(toBatch.begin(), toBatch.begin() + actualBatchSize);

  unsigned int from = asmc::getFromPosition(data.geneticPositions, startBatch, mCmBurnIn);
  unsigned int to = asmc::getToPosition(data.geneticPositions, endBatch, mCmBurnIn);

  const bool fastsmcWithHashing = decodingParams.FastSMC && decodingParams.hashing;
  const bool fromOrTo = from > 0u || to < sequenceLength;
  const bool makeBitsOnlyOnSubsequence = fastsmcWithHashing || fromOrTo;
  if (makeBitsOnlyOnSubsequence) {
    for (auto& obs : obsBatch) {
      makeBits(obs, from, to);
    }
  }

  // fill to size divisible by VECX
  while (obsBatch.size() % VECX != 0) {
    obsBatch.push_back(obsBatch.back());
  }

  auto paddedBatchSize = obsBatch.size();

  // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x paddedBatchSize]
  decodeBatch(obsBatch, from, to);
  augmentSumOverPairs(obsBatch, actualBatchSize, paddedBatchSize, from, to);

  if ((m_calculatePerPairMAP || m_calculatePerPairPosteriorMean) && !decodingParams.FastSMC) {
    writePerPairOutput(actualBatchSize, paddedBatchSize, obsBatch);
  }

  if (decodingParams.FastSMC) {
    writePerPairOutputFastSMC(actualBatchSize, paddedBatchSize, obsBatch);
  }

  obsBatch.clear();
}

// decode a batch
void HMM::decodeBatch(const vector<PairObservations>& obsBatch, const unsigned from, const unsigned to)
{

  int curBatchSize = static_cast<int>(obsBatch.size());

  Eigen::ArrayXf obsIsZeroBatch(sequenceLength * curBatchSize);
  Eigen::ArrayXf obsIsTwoBatch(sequenceLength * curBatchSize);

  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      obsIsZeroBatch[pos * curBatchSize + v] = (!obsBatch[v].obsBits[pos - from] ? 1.0f : 0.0f);
      obsIsTwoBatch[pos * curBatchSize + v] = (obsBatch[v].homMinorBits[pos - from] ? 1.0f : 0.0f);
    }
  }

  auto t0 = std::chrono::high_resolution_clock::now();

  // run forward
  forwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksForward += t1 - t0;

  // run backward
  backwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);

  auto t2 = std::chrono::high_resolution_clock::now();
  ticksBackward += t2 - t1;

  // combine (alpha * beta), normalize and store
  Eigen::Map<Eigen::ArrayXf> scale(obsIsZeroBatch.data(), obsIsZeroBatch.size()); // reuse buffer but rename to be less confusing
  scale.setZero();
#ifdef NO_SSE
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        long int ind = (pos * states + k) * curBatchSize + v;
        m_alphaBuffer[ind] *= m_betaBuffer[ind];
        scale[pos * curBatchSize + v] += m_alphaBuffer[ind];
      }
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      scale[pos * curBatchSize + v] = 1.0f / scale[pos * curBatchSize + v];
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        m_alphaBuffer[(pos * states + k) * curBatchSize + v] *= scale[pos * curBatchSize + v];
      }
    }
  }
#else
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        FLOAT prod = MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&m_betaBuffer[ind]));
        STORE(&m_alphaBuffer[ind], prod);
        STORE(&scale[pos * curBatchSize + v], ADD(LOAD(&scale[pos * curBatchSize + v]), prod));
      }
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v += VECX) {
      long int ind = pos * curBatchSize + v;
      STORE(&scale[ind], RECIPROCAL(LOAD(&scale[ind])));
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        STORE(&m_alphaBuffer[ind], MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&scale[pos * curBatchSize + v])));
      }
    }
  }
#endif

  auto t3 = std::chrono::high_resolution_clock::now();
  ticksCombine += t3 - t2;
}

// forward step
void HMM::forwardBatch(Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                       int curBatchSize, const unsigned from, const unsigned to)
{

  assert(curBatchSize % VECX == 0);

  Eigen::ArrayXf alphaC(states * curBatchSize);
  Eigen::ArrayXf AU(curBatchSize);
  Eigen::ArrayXf sums(curBatchSize);

  // fill pos=0 in alpha
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      float firstEmission = emission1AtSite[from][k] +
                            emission0minus1AtSite[from][k] * obsIsZeroBatch[from * curBatchSize + v] +
                            emission2minus0AtSite[from][k] * obsIsTwoBatch[from * curBatchSize + v];
      m_alphaBuffer[(states * from + k) * curBatchSize + v] = m_decodingQuant.initialStateProb[k] * firstEmission;
    }
  }

  Eigen::Map<Eigen::ArrayXf> currentAlpha(&m_alphaBuffer[states * from * curBatchSize], states * curBatchSize);
  asmc::calculateScalingBatch(currentAlpha, m_scalingBuffer, sums, curBatchSize, states);
  asmc::applyScalingBatch(currentAlpha, m_scalingBuffer, curBatchSize, states);

  // Induction Step:
  float lastGeneticPos = data.geneticPositions[from];
  int lastPhysicalPos = data.physicalPositions[from];

  for (long int pos = from + 1; pos < to; pos++) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(data.geneticPositions[pos] - lastGeneticPos, precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    Eigen::Map<Eigen::ArrayXf> previousAlpha(&m_alphaBuffer[(pos - 1) * states * curBatchSize], states * curBatchSize);
    Eigen::Map<Eigen::ArrayXf> nextAlpha(&m_alphaBuffer[pos * states * curBatchSize], states * curBatchSize);

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getNextAlphaBatched(recDistFromPreviousMinusOne, alphaC, curBatchSize, previousAlpha, pos, m_allZeros, m_allZeros,
                          AU, nextAlpha, homozEmission, homozEmission, homozEmission);
      previousAlpha = nextAlpha;
      getNextAlphaBatched(currentRecRate, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch, AU,
                          nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
    } else {
      getNextAlphaBatched(recDistFromPrevious, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch,
                          AU, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
    }

    if (pos % scalingSkip == 0) {
      asmc::calculateScalingBatch(nextAlpha, m_scalingBuffer, sums, curBatchSize, states);
      asmc::applyScalingBatch(nextAlpha, m_scalingBuffer, curBatchSize, states);
    }
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
}

// compute next alpha vector in linear time
void HMM::getNextAlphaBatched(float recDistFromPrevious, Eigen::Ref<Eigen::ArrayXf> alphaC, int curBatchSize,
                              Eigen::Ref<Eigen::ArrayXf> previousAlpha, uint pos,
                              Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch, Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch,
                              Eigen::Ref<Eigen::ArrayXf> AU, Eigen::Ref<Eigen::ArrayXf> nextAlpha,
                              const std::vector<float>& emission1AtSite, const std::vector<float>& emission0minus1AtSite,
                              const std::vector<float>& emission2minus0AtSite)
{

  const float* B = &m_decodingQuant.Bvectors.at(recDistFromPrevious)[0];
  const float* U = &m_decodingQuant.Uvectors.at(recDistFromPrevious)[0];
  const float* D = &m_decodingQuant.Dvectors.at(recDistFromPrevious)[0];

  memcpy(&alphaC[(states - 1) * curBatchSize], &previousAlpha[(states - 1) * curBatchSize],
         curBatchSize * sizeof(alphaC[0]));

  for (int k = states - 2; k >= 0; k--) {
#ifdef NO_SSE
    for (int v = 0; v < curBatchSize; v++) {
      alphaC[k * curBatchSize + v] = alphaC[(k + 1) * curBatchSize + v] + previousAlpha[k * curBatchSize + v];
    }
#else
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT alphaC_kp1_v = LOAD(&alphaC[(k + 1) * curBatchSize + v]);
      FLOAT prevAlpha_k_v = LOAD(&previousAlpha[k * curBatchSize + v]);
      STORE(&alphaC[k * curBatchSize + v], ADD(alphaC_kp1_v, prevAlpha_k_v));
    }
#endif
  }

  AU.setZero();
  for (int k = 0; k < states; k++) {

#ifdef NO_SSE
    for (int v = 0; v < curBatchSize; v++) {
      if (k)
        AU[v] = U[k - 1] * previousAlpha[(k - 1) * curBatchSize + v] + m_decodingQuant.columnRatios[k - 1] * AU[v];
      float term = AU[v] + D[k] * previousAlpha[k * curBatchSize + v];
      if (k < states - 1) {
        term += B[k] * alphaC[(k + 1) * curBatchSize + v];
      }
      float currentEmission_k = emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZeroBatch[pos * curBatchSize + v] +
                                emission2minus0AtSite[k] * obsIsTwoBatch[pos * curBatchSize + v];
      nextAlpha[k * curBatchSize + v] = currentEmission_k * term;
    }
#else
#ifdef AVX512
    FLOAT D_k = LOAD1(D[k]);
    FLOAT B_k = LOAD1(B[k]);
    FLOAT em1Prob_k = LOAD1(emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(emission2minus0AtSite[k]);
#else
    FLOAT D_k = LOAD1(&D[k]);
    FLOAT B_k = LOAD1(&B[k]);
    FLOAT em1Prob_k = LOAD1(&emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(&emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(&emission2minus0AtSite[k]);
#endif

    FLOAT Ukm1, colRatios_km1;
    if (k) {
#ifdef AVX512
      Ukm1 = LOAD1(U[k - 1]);
      colRatios_km1 = LOAD1(m_decodingQuant.columnRatios[k - 1]);
#else
      Ukm1 = LOAD1(&U[k - 1]);
      colRatios_km1 = LOAD1(&m_decodingQuant.columnRatios[k - 1]);
#endif
    }
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT AU_v;
      if (k) {
        FLOAT term1 = MULT(Ukm1, LOAD(&previousAlpha[(k - 1) * curBatchSize + v]));
        FLOAT term2 = MULT(colRatios_km1, LOAD(&AU[v]));
        AU_v = ADD(term1, term2);
        STORE(&AU[v], AU_v);
      } else
        AU_v = LOAD(&AU[v]);

      FLOAT term = ADD(AU_v, MULT(D_k, LOAD(&previousAlpha[k * curBatchSize + v])));
      if (k < states - 1) { // TODO: just extend B and alphaC?
        term = ADD(term, MULT(B_k, LOAD(&alphaC[(k + 1) * curBatchSize + v])));
      }
      FLOAT currentEmission_k =
          ADD(ADD(em1Prob_k, MULT(em0minus1Prob_k, LOAD(&obsIsZeroBatch[pos * curBatchSize + v]))),
              MULT(em2minus0Prob_k, LOAD(&obsIsTwoBatch[pos * curBatchSize + v])));
      if (v == 0) {
      }
      STORE(&nextAlpha[k * curBatchSize + v], MULT(currentEmission_k, term));
    }
#endif
  }
}

// backward step
void HMM::backwardBatch(Eigen::ArrayXf obsIsZeroBatch, Eigen::ArrayXf obsIsTwoBatch, int curBatchSize,
                        const unsigned from, const unsigned to)
{

  // fill pos=sequenceLenght-1 in beta
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      m_betaBuffer[((to - 1) * states + k) * curBatchSize + v] = 1.0f;
    }
  }

  Eigen::ArrayXf sums(curBatchSize);
  Eigen::Map<Eigen::ArrayXf> currentBeta(&m_betaBuffer[states * (to - 1) * curBatchSize], states * curBatchSize);

  asmc::calculateScalingBatch(currentBeta, m_scalingBuffer, sums, curBatchSize, states);
  asmc::applyScalingBatch(currentBeta, m_scalingBuffer, curBatchSize, states);

  // Induction Step:
  Eigen::ArrayXf BL(curBatchSize);
  Eigen::ArrayXf BU(states * curBatchSize);
  Eigen::ArrayXf vec(states * curBatchSize);

  float lastGeneticPos = data.geneticPositions[to - 1];
  int lastPhysicalPos = data.physicalPositions[to - 1];

  for (long int pos = to - 2; pos >= from; pos--) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(lastGeneticPos - data.geneticPositions[pos], precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);

    Eigen::Map<Eigen::ArrayXf> previousBeta(&m_betaBuffer[pos * states * curBatchSize], states * curBatchSize);
    Eigen::Map<Eigen::ArrayXf> lastComputedBeta(&m_betaBuffer[(pos + 1) * states * curBatchSize], states * curBatchSize);

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1, precision);
      float recDistFromPreviousMinusOne = asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getPreviousBetaBatched(recDistFromPreviousMinusOne, curBatchSize, lastComputedBeta, pos, m_allZeros, m_allZeros,
                             vec, BU, BL, previousBeta, homozEmission, homozEmission, homozEmission);
      lastComputedBeta = previousBeta;
      getPreviousBetaBatched(currentRecRate, curBatchSize, lastComputedBeta, pos, obsIsZeroBatch, obsIsTwoBatch, vec,
                             BU, BL, previousBeta, emission1AtSite[pos + 1], emission0minus1AtSite[pos + 1],
                             emission2minus0AtSite[pos + 1]);
    } else {
      getPreviousBetaBatched(recDistFromPrevious, curBatchSize, lastComputedBeta, pos, obsIsZeroBatch, obsIsTwoBatch,
                             vec, BU, BL, previousBeta, emission1AtSite[pos + 1], emission0minus1AtSite[pos + 1],
                             emission2minus0AtSite[pos + 1]);
    }
    if (pos % scalingSkip == 0) {
      // normalize betas using alpha scaling
      asmc::calculateScalingBatch(previousBeta, m_scalingBuffer, sums, curBatchSize, states);
      asmc::applyScalingBatch(previousBeta, m_scalingBuffer, curBatchSize, states);
    }
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
}

// compute previous beta vector in linear time
void HMM::getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, Eigen::Ref<Eigen::ArrayXf> lastComputedBeta,
                                 int pos, Eigen::Ref<Eigen::ArrayXf> obsIsZeroBatch,
                                 Eigen::Ref<Eigen::ArrayXf> obsIsTwoBatch, Eigen::Ref<Eigen::ArrayXf> vec,
                                 Eigen::Ref<Eigen::ArrayXf> BU, Eigen::Ref<Eigen::ArrayXf> BL,
                                 Eigen::Ref<Eigen::ArrayXf> currentBeta, const std::vector<float>& emission1AtSite,
                                 const std::vector<float>& emission0minus1AtSite,
                                 const std::vector<float>& emission2minus0AtSite)
{
  const vector<float>& B = m_decodingQuant.Bvectors.at(recDistFromPrevious);
  const vector<float>& U = m_decodingQuant.Uvectors.at(recDistFromPrevious);
  const vector<float>& RR = m_decodingQuant.rowRatioVectors.at(recDistFromPrevious);
  const vector<float>& D = m_decodingQuant.Dvectors.at(recDistFromPrevious);
#ifdef NO_SSE

  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      float currentEmission_k = emission1AtSite[k] +
                                emission0minus1AtSite[k] * obsIsZeroBatch[(pos + 1) * curBatchSize + v] +
                                emission2minus0AtSite[k] * obsIsTwoBatch[(pos + 1) * curBatchSize + v];
      vec[k * curBatchSize + v] = lastComputedBeta[k * curBatchSize + v] * currentEmission_k;
    }
  }

#else
  for (int k = 0; k < states; k++) {
#ifdef AVX512
    FLOAT em1Prob_k = LOAD1(emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(emission2minus0AtSite[k]);
#else
    FLOAT em1Prob_k = LOAD1(&emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(&emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(&emission2minus0AtSite[k]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT currentEmission_k =
          ADD(ADD(em1Prob_k, MULT(em0minus1Prob_k, LOAD(&obsIsZeroBatch[(pos + 1) * curBatchSize + v]))),
              MULT(em2minus0Prob_k, LOAD(&obsIsTwoBatch[(pos + 1) * curBatchSize + v])));
      STORE(&vec[k * curBatchSize + v], MULT(LOAD(&lastComputedBeta[k * curBatchSize + v]), currentEmission_k));
    }
  }
#endif

  BU.setZero();
#ifdef NO_SSE
  for (int k = states - 2; k >= 0; k--)
    for (int v = 0; v < curBatchSize; v++)
      BU[k * curBatchSize + v] = U[k] * vec[(k + 1) * curBatchSize + v] + RR[k] * BU[(k + 1) * curBatchSize + v];
#else
  for (int k = states - 2; k >= 0; k--) {
#ifdef AVX512
    FLOAT U_k = LOAD1(U[k]);
    FLOAT RR_k = LOAD1(RR[k]);
#else
    FLOAT U_k = LOAD1(&U[k]);
    FLOAT RR_k = LOAD1(&RR[k]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT term1 = MULT(U_k, LOAD(&vec[(k + 1) * curBatchSize + v]));
      FLOAT term2 = MULT(RR_k, LOAD(&BU[(k + 1) * curBatchSize + v]));
      STORE(&BU[k * curBatchSize + v], ADD(term1, term2));
    }
  }
#endif

  BL.setZero();
#ifdef NO_SSE
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      if (k)
        BL[v] += B[k - 1] * vec[(k - 1) * curBatchSize + v];
      currentBeta[k * curBatchSize + v] = BL[v] + D[k] * vec[k * curBatchSize + v] + BU[k * curBatchSize + v];
    }
  }
#else
  for (int k = 0; k < states; k++) {
#ifdef AVX512
    FLOAT D_k = LOAD1(D[k]);
    FLOAT B_km1;
    if (k)
      B_km1 = LOAD1(B[k - 1]);
#else
    FLOAT D_k = LOAD1(&D[k]);
    FLOAT B_km1;
    if (k)
      B_km1 = LOAD1(&B[k - 1]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT BL_v = LOAD(&BL[v]);
      if (k) {
        BL_v = ADD(BL_v, MULT(B_km1, LOAD(&vec[(k - 1) * curBatchSize + v])));
        STORE(&BL[v], BL_v);
      }
      STORE(&currentBeta[k * curBatchSize + v],
            ADD(BL_v, ADD(MULT(D_k, LOAD(&vec[k * curBatchSize + v])), LOAD(&BU[k * curBatchSize + v]))));
    }
  }
#endif
}

// --posteriorSums
void HMM::augmentSumOverPairs(vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize, unsigned from, unsigned to)
{

  auto t0 = std::chrono::high_resolution_clock::now();

  if (!decodingParams.doPosteriorSums && !decodingParams.doMajorMinorPosteriorSums)
    return;

  for (long int pos = from; pos < to; pos++) {

    long int offset = pos - static_cast<long int>(from);

    for (int k = 0; k < states; k++) {
      float sum = 0;
      float sum00 = 0;
      float sum01 = 0;
      float sum11 = 0;
      for (int v = 0; v < actualBatchSize; v++) { // only loop over actual (not padding) pairs
        float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
        if (decodingParams.doPosteriorSums) {
          sum += posterior_pos_state_pair;
        }
        if (decodingParams.doMajorMinorPosteriorSums) {
          if (obsBatch.at(v).homMinorBits.at(offset) == 1)
            sum11 += posterior_pos_state_pair;
          else if (obsBatch.at(v).homMinorBits.at(offset) == 0)
            sum00 += posterior_pos_state_pair;
          else
            sum01 += posterior_pos_state_pair;
        }
      }
      if (decodingParams.doPosteriorSums) {
        m_decodingReturnValues.sumOverPairs(offset, k) += sum;
      }
      if (decodingParams.doMajorMinorPosteriorSums) {
        m_decodingReturnValues.sumOverPairs00(offset, k) += sum00;
        m_decodingReturnValues.sumOverPairs01(offset, k) += sum01;
        m_decodingReturnValues.sumOverPairs11(offset, k) += sum11;
      }
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksSumOverPairs += t1 - t0;
}

float HMM::getPosteriorMean(const vector<float>& posterior)
{

  float normalization = 1.f / std::accumulate(posterior.begin(), posterior.end(), 0.f);

  float posteriorMean = 0.f;
  for (auto k = 0ul; k < posterior.size(); k++) {
    posteriorMean += normalization * posterior[k] * m_decodingQuant.expectedTimes[k];
  }
  return posteriorMean;
}

float HMM::getMAP(vector<float> posterior)
{
  vector<float> ratioPriorPosterior(posterior.size());
  std::transform(posterior.begin(), posterior.end(), m_decodingQuant.initialStateProb.begin(),
                 ratioPriorPosterior.begin(), std::divides<>());
  const auto MAP_location = std::distance(ratioPriorPosterior.begin(),
                                          std::max_element(ratioPriorPosterior.begin(), ratioPriorPosterior.end()));
  return m_decodingQuant.expectedTimes[MAP_location];
}

// write an IBD segment into output file
void HMM::writePairIBD(const PairObservations& obs, unsigned int posStart, unsigned int posEnd, float prob,
                       const vector<float>& posterior)
{
  nbSegmentsDetected++;
  if (!decodingParams.BIN_OUT) {

    std::stringstream record;
    record << std::setprecision(std::numeric_limits<float>::digits10 + 1);

    record << data.FamIDList[obs.iInd] << '\t' << data.IIDList[obs.iInd] << '\t' << static_cast<int>(obs.iHap) << '\t'
           << data.FamIDList[obs.jInd] << '\t' << data.IIDList[obs.jInd] << '\t' << static_cast<int>(obs.jHap) << '\t'
           << data.chrNumber;

    record << '\t' << data.physicalPositions[posStart] << '\t' << data.physicalPositions[posEnd];

    if (decodingParams.outputIbdSegmentLength) {
      const float length_cM = 100.f * (data.geneticPositions[posEnd] - data.geneticPositions[posStart]);
      record << '\t' << length_cM;
    }

    if (!decodingParams.hashingOnly) {
      const double ibd_score = prob / static_cast<double>(posEnd - posStart + 1u);
      record << '\t' << ibd_score;
    }

    if (decodingParams.doPerPairPosteriorMean) {
      record << '\t' << getPosteriorMean(posterior);
    }

    if (decodingParams.doPerPairMAP) {
      record << '\t' << getMAP(posterior);
    }

    record << '\n';

    const std::string record_str = record.str();
    gzwrite(gzoutIBD, record_str.c_str(), record_str.size());

  } else {
    unsigned int ind[2];
    ind[0] = obs.iInd;
    ind[1] = obs.jInd;
    std::uint_least8_t hap[2];
    hap[0] = obs.iHap;
    hap[1] = obs.jHap;
    int pos[2];
    pos[0] = data.physicalPositions[posStart];
    pos[1] = data.physicalPositions[posEnd];
    gzwrite(gzoutIBD, (char*)&ind[0], sizeof(unsigned int));
    gzwrite(gzoutIBD, &hap[0], sizeof(std::uint_least8_t));
    gzwrite(gzoutIBD, (char*)&ind[1], sizeof(unsigned int));
    gzwrite(gzoutIBD, &hap[1], sizeof(std::uint_least8_t));
    gzwrite(gzoutIBD, (char*)&pos[0], sizeof(int));
    gzwrite(gzoutIBD, (char*)&pos[1], sizeof(int));
    if (decodingParams.outputIbdSegmentLength) {
      const float length_cM = 100.f * (data.geneticPositions[posEnd] - data.geneticPositions[posStart]);
      gzwrite(gzoutIBD, (char*)&length_cM, sizeof(float));
    }
    if (!decodingParams.hashingOnly) {
      const auto ibd_score = static_cast<float>(prob / static_cast<double>(posEnd - posStart + 1u));
      gzwrite(gzoutIBD, (char*)&ibd_score, sizeof(float));
    }
    if (decodingParams.doPerPairPosteriorMean) {
      float postMean = getPosteriorMean(posterior);
      gzwrite(gzoutIBD, (char*)&postMean, sizeof(float));
    }
    if (decodingParams.doPerPairMAP) {
      float map = getMAP(posterior);
      gzwrite(gzoutIBD, (char*)&map, sizeof(float));
    }
  }
}

void HMM::writePerPairOutputFastSMC(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
{
  for (int v = 0; v < actualBatchSize; v++) {
    bool isIBD = false, isIBD1 = false, isIBD2 = false,
         isIBD3 = false; // true if previous position is IBD, false otherwise
    unsigned int startIBD = 0, startIBD1 = 0, startIBD2 = 0, startIBD3 = 0;
    unsigned int endIBD = 0, endIBD1 = 0, endIBD2 = 0, endIBD3 = 0;
    float posteriorIBD = 0; // cumulative posterior on an IBD segment
    vector<float> posterior;
    vector<float> sum_posterior_per_state;
    vector<float> prev_sum_posterior_per_state;

    if (decodingParams.doPerPairPosteriorMean || decodingParams.doPerPairMAP) {
      for (uint k = 0; k < ageThreshold; k++) {
        posterior.push_back(0);
        sum_posterior_per_state.push_back(0);
        prev_sum_posterior_per_state.push_back(0);
      }
    }

    if (decodingParams.hashing) {
      // remove these 2 lines if you want the preprocessing step to be less permissive
      // TODO : add a flag for this option
      fromBatch[v] = startBatch;
      toBatch[v] = endBatch;
    }

    for (unsigned int pos = fromBatch[v]; pos < toBatch[v]; pos++) {
      float sum = 0;

      if (decodingParams.doPerPairPosteriorMean || decodingParams.doPerPairMAP) {
        for (uint k = 0; k < ageThreshold; k++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
          posterior[k] = posterior_pos_state_pair;
          prev_sum_posterior_per_state[k] = sum_posterior_per_state[k];
          sum_posterior_per_state[k] += posterior_pos_state_pair;
          if (k < stateThreshold) {
            sum += posterior_pos_state_pair;
          }
        }
      } else {
        for (uint k = 0; k < stateThreshold; k++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
          sum += posterior_pos_state_pair;
        }
      }

      if (sum >= 1000 * probabilityThreshold) {
        if (!isIBD) {
          startIBD = pos;
          sum_posterior_per_state = posterior;
          if (pos > fromBatch[v] && isIBD1) {
            endIBD1 = pos - 1;
            writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD2) {
            endIBD2 = pos - 1;
            writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD3) {
            endIBD3 = pos - 1;
            writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state);
          }
          posteriorIBD = sum;
        } else {
          posteriorIBD += sum;
        }
        if (pos == toBatch[v] - 1) {
          endIBD = toBatch[v] - 1;
          writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, sum_posterior_per_state);
          posteriorIBD = 0;
        }
        isIBD = true;
        isIBD1 = false, isIBD2 = false, isIBD3 = false;
      } else if (sum >= 100 * probabilityThreshold) {
        if (!isIBD1) {
          startIBD1 = pos;
          sum_posterior_per_state = posterior;
          if (pos > fromBatch[v] && isIBD) {
            endIBD = pos - 1;
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD2) {
            endIBD2 = pos - 1;
            writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD3) {
            endIBD3 = pos - 1;
            writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state);
          }
          posteriorIBD = sum;
        } else {
          posteriorIBD += sum;
        }
        if (pos == toBatch[v] - 1) {
          endIBD1 = toBatch[v] - 1;
          writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, sum_posterior_per_state);
          posteriorIBD = 0;
        }
        isIBD = false, isIBD2 = false, isIBD3 = false;
        isIBD1 = true;
      } else if (sum >= 10 * probabilityThreshold) {
        if (!isIBD2) {
          startIBD2 = pos;
          sum_posterior_per_state = posterior;
          if (pos > fromBatch[v] && isIBD1) {
            endIBD1 = pos - 1;
            writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD) {
            endIBD = pos - 1;
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD3) {
            endIBD3 = pos - 1;
            writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state);
          }
          posteriorIBD = sum;
        } else {
          posteriorIBD += sum;
        }
        if (pos == toBatch[v] - 1) {
          endIBD2 = toBatch[v] - 1;
          writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, sum_posterior_per_state);
          posteriorIBD = 0;
        }
        isIBD = false, isIBD1 = false, isIBD3 = false;
        isIBD2 = true;
      } else if (sum >= probabilityThreshold) {
        if (!isIBD3) {
          startIBD3 = pos;
          sum_posterior_per_state = posterior;
          if (pos > fromBatch[v] && isIBD1) {
            endIBD1 = pos - 1;
            writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD) {
            endIBD = pos - 1;
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state);
          } else if (pos > fromBatch[v] && isIBD2) {
            endIBD2 = pos - 1;
            writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state);
          }
          posteriorIBD = sum;
        } else {
          posteriorIBD += sum;
        }
        if (pos == toBatch[v] - 1) {
          endIBD3 = toBatch[v] - 1;
          writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, sum_posterior_per_state);
          posteriorIBD = 0;
        }
        isIBD = false, isIBD1 = false, isIBD2 = false;
        isIBD3 = true;
      } else {
        if (isIBD) {
          endIBD = pos - 1;
          writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state);
          posteriorIBD = 0;
        } else if (isIBD1) {
          endIBD1 = pos - 1;
          writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state);
          posteriorIBD = 0;
        } else if (isIBD2) {
          endIBD2 = pos - 1;
          writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state);
          posteriorIBD = 0;
        } else if (isIBD3) {
          endIBD3 = pos - 1;
          writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state);
          posteriorIBD = 0;
        }
        isIBD = false, isIBD1 = false, isIBD2 = false, isIBD3 = false;
      }
    }
  }
}

// will eventually write binary output instead of gzipped
void HMM::writePerPairOutput(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
{
  auto t0 = std::chrono::high_resolution_clock::now();

  const long int seqFrom = static_cast<long int>(*std::min_element(fromBatch.begin(), fromBatch.end()));
  const long int seqTo = static_cast<long int>(*std::max_element(toBatch.begin(), toBatch.end()));
  assert(seqTo > seqFrom);
  const unsigned seqLen = seqTo - seqFrom;

  Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1> posteriors;

  // Calculate per pair posterior mean
  if (m_calculatePerPairPosteriorMean) {
    meanPost.setZero();

    // Store posteriors for batch, but only if requested (large amount of data)
    if (m_storePerPairPosterior) {
      posteriors.resize(actualBatchSize);
      for (Eigen::Index batchIdx = 0ll; batchIdx < actualBatchSize; ++batchIdx) {
        posteriors(batchIdx).resize(states, seqLen);
      }
    }

    for (long int pos = seqFrom; pos < seqTo; pos++) {
      const long int offset = pos - seqFrom;
      for (int k = 0; k < states; k++) {
        for (int batchIdx = 0; batchIdx < actualBatchSize; batchIdx++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + batchIdx];
          float postValue = posterior_pos_state_pair;
          meanPost(batchIdx, offset) += postValue * expectedCoalTimes[k];
          if (m_storePerPairPosterior) {
            posteriors(batchIdx)(k, offset) = postValue;
          }
          if (m_storeSumOfPosterior) {
            m_decodePairsReturnStruct.sumOfPosteriors(k, offset) += postValue;
          }
        }
      }
    }
  }

  // Calculate per pair MAP
  if (m_calculatePerPairMAP) {
    MAP.setZero();
    for (long int pos = seqFrom; pos < seqTo; pos++) {
      const long int offset = pos - seqFrom;
      currentMAPValue.setZero();
      for (int k = 0; k < states; k++) {
        for (int batchIdx = 0; batchIdx < actualBatchSize; batchIdx++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + batchIdx];
          if (currentMAPValue[batchIdx] < posterior_pos_state_pair) {
            MAP(batchIdx, offset) = k;
            currentMAPValue[batchIdx] = posterior_pos_state_pair;
          }
        }
      }
    }
  }

  // Write per pair posterior mean to file
  if (m_writePerPairPosteriorMean) {
    foutPosteriorMeanPerPair << meanPost.topRows(actualBatchSize).format(m_eigenOutputFormat);
  }

  // Write per pair MAP to file
  if (m_writePerPairMAP) {
    foutMAPPerPair << MAP.topRows(actualBatchSize).format(m_eigenOutputFormat);
  }

  // Store per pair information, if required
  if (m_storePerPairMAP || m_storePerPairPosteriorMean) {
    for (int batchIdx = 0; batchIdx < actualBatchSize; ++batchIdx) {

      // Get the index: this is the row we need to write into
      const auto outIdx = static_cast<Eigen::Index>(m_decodePairsReturnStruct.getNumWritten());

      // Record the index information
      unsigned long hapA = asmc::dipToHapId(obsBatch[batchIdx].iInd, obsBatch[batchIdx].iHap);
      unsigned long hapB = asmc::dipToHapId(obsBatch[batchIdx].jInd, obsBatch[batchIdx].jHap);
      std::string indIdA =
          asmc::indPlusHapToCombinedId(data.IIDList.at(obsBatch[batchIdx].iInd), obsBatch[batchIdx].iHap);
      std::string indIdB =
          asmc::indPlusHapToCombinedId(data.IIDList.at(obsBatch[batchIdx].jInd), obsBatch[batchIdx].jHap);

      m_decodePairsReturnStruct.perPairIndices.at(outIdx) = std::make_tuple(hapA, indIdA, hapB, indIdB);

      if (m_storePerPairPosterior) {
        m_decodePairsReturnStruct.perPairPosteriors.at(outIdx) = posteriors(batchIdx);
      }

      if (m_storePerPairPosteriorMean) {
        m_decodePairsReturnStruct.perPairPosteriorMeans.row(outIdx) = meanPost.row(batchIdx);
      }

      if (m_storePerPairPosteriorMean) {
        m_decodePairsReturnStruct.perPairMAPs.row(outIdx) = MAP.row(batchIdx);
      }

      // Increment
      m_decodePairsReturnStruct.incrementNumWritten();
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksOutputPerPair += t1 - t0;
}

// *********************************************************************
// non-batched computations (for debugging and pedagogical reasons only)
// *********************************************************************

vector<vector<float>> HMM::decode(const PairObservations& observations)
{
  return decode(observations, 0, sequenceLength);
}

vector<vector<float>> HMM::decode(const PairObservations& observations, const unsigned from, const unsigned to)
{

  auto t0 = std::chrono::high_resolution_clock::now();
  vector<vector<float>> forwardOut = forward(observations, from, to);
  auto t1 = std::chrono::high_resolution_clock::now();
  ticksForward += t1 - t0;

  vector<vector<float>> backwardOut = backward(observations, from, to);
  auto t2 = std::chrono::high_resolution_clock::now();
  ticksBackward += t2 - t1;

  vector<vector<float>> posterior = asmc::elementWiseMultMatrixMatrix(forwardOut, backwardOut);
  posterior = asmc::normalizeMatrixColumns(posterior);

  auto t3 = std::chrono::high_resolution_clock::now();
  ticksCombine += t3 - t2;

  if (decodingParams.doPosteriorSums) {
    for (uint k = 0; k < m_decodingQuant.states; k++) {
      for (long int pos = 0; pos < sequenceLength; pos++) {
        m_decodingReturnValues.sumOverPairs(pos, k) += posterior[k][pos];
      }
    }
  }
  return posterior;
}

// Instead of returning a whole posterior, return the MAP and posterior mean
pair<vector<float>, vector<float>> HMM::decodeSummarize(const PairObservations& observations)
{
  vector<vector<float>> posterior = decode(observations);
  size_t num_discretizations = m_decodingQuant.expectedTimes.size();
  assert(num_discretizations == posterior.size());

  vector<float> posterior_map(posterior[0].size(), 0);
  vector<float> posterior_max_so_far(posterior[0].size(), 0);
  vector<float> posterior_mean(posterior[0].size(), 0);
  for (int i = 0; i < posterior.size(); ++i) {
    for (int j = 0; j < posterior[0].size(); ++j) {
      posterior_mean[j] += posterior[i][j] * m_decodingQuant.expectedTimes[i];
      if (posterior[i][j] > posterior_max_so_far[j]) {
        posterior_max_so_far[j] = posterior[i][j];
        posterior_map[j] = m_decodingQuant.expectedTimes[i];
      }
    }
  }
  return std::make_pair(posterior_map, posterior_mean);
}

vector<float> HMM::getEmission(int pos, int distinguished, int undistinguished, int emissionIndex)
{
  vector<float> emission;
  if (!useCSFSatThisPosition[pos]) {
    // this position is not a CSFS position, use compressed or classic
    emission = decodingParams.decodingSequence ? m_decodingQuant.classicEmissionTable[emissionIndex]
                                               : m_decodingQuant.compressedEmissionTable[emissionIndex];
  } else {
    // this position is a CSFS position
    if (decodingParams.foldData) {
      emission = decodingParams.decodingSequence
                     ? m_decodingQuant.foldedCSFSmap[undistinguished][emissionIndex]
                     : m_decodingQuant.foldedAscertainedCSFSmap[undistinguished][emissionIndex];
    } else {
      emission = decodingParams.decodingSequence ? m_decodingQuant.CSFSmap[undistinguished][distinguished]
                                                 : m_decodingQuant.ascertainedCSFSmap[undistinguished][distinguished];
    }
  }
  return emission;
}

// forward step
vector<vector<float>> HMM::forward(const PairObservations& observations, const unsigned from, const unsigned to)
{

  vector<vector<float>> alpha(states, vector<float>(sequenceLength));

  uint emissionIndex = observations.obsBits[from] ? 1 : 0;
  // if both samples are carriers, there are two distinguished, otherwise, it's the or
  // (use previous xor). This affects the number of undistinguished for the site
  uint distinguished = observations.homMinorBits[from] ? 2 : emissionIndex;
  uint undistinguished = useCSFSatThisPosition[from] ? undistinguishedCounts[from][distinguished] : -1;
  vector<float> emission = getEmission(from, distinguished, undistinguished, emissionIndex);

  vector<float> firstAlpha = asmc::elementWiseMultVectorVector(m_decodingQuant.initialStateProb, emission);

  // cumpute scaling (sum of current alpha vector)
  float m_scalingBuffer = 1.f / asmc::getSumOfVector(firstAlpha);
  // normalize current alpha vector to 1
  firstAlpha = asmc::elementWiseMultVectorScalar(firstAlpha, m_scalingBuffer);

  asmc::fillMatrixColumn(alpha, firstAlpha, from);
  // Induction Step:
  vector<float> nextAlpha(states);
  vector<float> alphaC(states + 1);
  vector<float> previousAlpha = firstAlpha;
  float lastGeneticPos = data.geneticPositions[from];
  int lastPhysicalPos = data.physicalPositions[from];

  for (long int pos = from + 1; pos < to; pos++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(data.geneticPositions[pos] - lastGeneticPos, precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos] ? 1.0f : 0.0f;

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = asmc::roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getNextAlpha(recDistFromPreviousMinusOne, alphaC, previousAlpha, nextAlpha, homozEmission, homozEmission,
                   homozEmission, 0.0f, 0.0f);
      previousAlpha = nextAlpha;
      getNextAlpha(currentRecRate, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                   emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
    } else {
      getNextAlpha(recDistFromPrevious, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos],
                   emission0minus1AtSite[pos], emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    if (pos % scalingSkip == 0) {
      // compute scaling (sum of current alpha vector)
      m_scalingBuffer = 1.f / asmc::getSumOfVector(nextAlpha);
      // normalize current alpha vector to 1
      nextAlpha = asmc::elementWiseMultVectorScalar(nextAlpha, m_scalingBuffer);
    }
    asmc::fillMatrixColumn(alpha, nextAlpha, pos);
    previousAlpha = nextAlpha;
    auto t2 = std::chrono::high_resolution_clock::now();
    t1sum += t1 - t0;
    t2sum += t2 - t1;
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
  return alpha;
}

void HMM::getNextAlpha(float recDistFromPrevious, vector<float>& alphaC, vector<float>& previousAlpha,
                       vector<float>& nextAlpha, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                       vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor)
{
  alphaC[m_decodingQuant.states - 1] = previousAlpha[m_decodingQuant.states - 1];
  for (int k = m_decodingQuant.states - 2; k >= 0; k--) {
    alphaC[k] = alphaC[k + 1] + previousAlpha[k];
  }
  const float* B = &m_decodingQuant.Bvectors.at(recDistFromPrevious)[0];
  const float* U = &m_decodingQuant.Uvectors.at(recDistFromPrevious)[0];
  const float* D = &m_decodingQuant.Dvectors.at(recDistFromPrevious)[0];
  float AUc = 0;
  for (uint k = 0; k < m_decodingQuant.states; k++) {
    if (k)
      AUc = U[k - 1] * previousAlpha[k - 1] + m_decodingQuant.columnRatios[k - 1] * AUc;
    float term = AUc + D[k] * previousAlpha[k];
    if (k < m_decodingQuant.states - 1)
      term += B[k] * alphaC[k + 1];
    float currentEmission_k =
        emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
    nextAlpha[k] = currentEmission_k * term;
  }
}

// backward step
vector<vector<float>> HMM::backward(const PairObservations& observations, const unsigned from, const unsigned to)
{
  vector<vector<float>> beta(states, vector<float>(sequenceLength));

  vector<float> lastBeta(states);
  for (uint i = 0; i < lastBeta.size(); i++) {
    lastBeta[i] = 1.f;
  }
  // normalize current alpha vector to 1
  float m_scalingBuffer = 1.f / asmc::getSumOfVector(lastBeta);
  lastBeta = asmc::elementWiseMultVectorScalar(lastBeta, m_scalingBuffer);
  asmc::fillMatrixColumn(beta, lastBeta, to - 1);
  // Induction Step:
  vector<float> currentBeta(states);
  vector<float> BL(states);
  vector<float> BU(states);
  vector<float> lastComputedBeta = lastBeta;
  float lastGeneticPos = data.geneticPositions[to - 1];
  int lastPhysicalPos = data.physicalPositions[to - 1];

  for (long int pos = to - 2; pos >= from; pos--) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(lastGeneticPos - data.geneticPositions[pos], precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos + 1] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos + 1] ? 1.0f : 0.0f;
    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getPreviousBeta(recDistFromPreviousMinusOne, lastComputedBeta, BL, BU, currentBeta, homozEmission, homozEmission,
                      homozEmission, 0.0f, 0.0f);
      lastComputedBeta = currentBeta;
      getPreviousBeta(currentRecRate, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                      emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
    } else {
      getPreviousBeta(recDistFromPrevious, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                      emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
    }
    if (pos % scalingSkip == 0) {
      m_scalingBuffer = 1.f / asmc::getSumOfVector(currentBeta);
      currentBeta = asmc::elementWiseMultVectorScalar(currentBeta, m_scalingBuffer);
    }
    asmc::fillMatrixColumn(beta, currentBeta, pos);
    lastComputedBeta = currentBeta;
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
  return beta;
}

void HMM::getPreviousBeta(float recDistFromPrevious, vector<float>& lastComputedBeta, vector<float>& BL,
                          vector<float>& BU, vector<float>& currentBeta, vector<float>& emission1AtSite,
                          vector<float>& emission0minus1AtSite, vector<float>& emission2minus0AtSite, float obsIsZero,
                          float obsIsHomMinor)
{
  vector<float> vec = vector<float>(states);
  for (int k = 0; k < states; k++) {
    float currentEmission_k =
        emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
    vec[k] = lastComputedBeta[k] * currentEmission_k;
  }
  // compute below table
  float sum = 0;
  const vector<float>& B = m_decodingQuant.Bvectors.at(recDistFromPrevious);
  for (int k = 1; k < states; k++) {
    sum += B[k - 1] * vec[k - 1];
    BL[k] = sum;
  }
  // compute above table
  const vector<float>& U = m_decodingQuant.Uvectors.at(recDistFromPrevious);
  const vector<float>& RR = m_decodingQuant.rowRatioVectors.at(recDistFromPrevious);
  for (int k = states - 2; k >= 0; k--) {
    BU[k] = vec[k + 1] * U[k] + RR[k] * BU[k + 1];
  }
  // put them together
  const vector<float>& D = m_decodingQuant.Dvectors.at(recDistFromPrevious);
  for (int k = 0; k < states; k++) {
    currentBeta[k] = BL[k] + vec[k] * D[k] + BU[k];
  }
}

const DecodingQuantities& HMM::getDecodingQuantities() const
{
  return m_decodingQuant;
}

DecodePairsReturnStruct& HMM::getDecodePairsReturnStruct()
{
  return m_decodePairsReturnStruct;
}

void HMM::updateOutputStructures() {
  m_calculatePerPairPosteriorMean =
      m_storePerPairPosteriorMean || m_writePerPairPosteriorMean || m_storePerPairPosterior;
  m_calculatePerPairMAP = m_storePerPairMAP || m_writePerPairMAP;

  const long int seqFrom = static_cast<long int>(*std::min_element(fromBatch.begin(), fromBatch.end()));
  const long int seqTo = static_cast<long int>(*std::max_element(toBatch.begin(), toBatch.end()));
  const unsigned seqLen = seqTo - seqFrom;

  if (m_calculatePerPairPosteriorMean) {
    meanPost.resize(m_batchSize, seqLen);

    if (expectedCoalTimes.empty() && !decodingParams.FastSMC) {
      if (!expectedCoalTimesFile.empty() && std::filesystem::is_regular_file(expectedCoalTimesFile)) {
        fmt::print("Reading expected coalescent times from {}\n", expectedCoalTimesFile);
        expectedCoalTimes = readExpectedTimesFromIntervalsFile(expectedCoalTimesFile.c_str());
      } else {
        fmt::print("Using expected coalescent times from {}\n", decodingParams.decodingQuantFile);
        expectedCoalTimes = m_decodingQuant.expectedTimes;
      }
    }
  }

  if (m_calculatePerPairMAP) {
    MAP.resize(m_batchSize, seqLen);
    currentMAPValue.resize(m_batchSize);
  }

  resetDecoding();
}

void HMM::setStorePerPairPosteriorMean(bool storePerPairPosteriorMean)
{
  m_storePerPairPosteriorMean = storePerPairPosteriorMean;
  updateOutputStructures();
}

void HMM::setWritePerPairPosteriorMean(bool writePerPairPosteriorMean)
{
  m_writePerPairPosteriorMean = writePerPairPosteriorMean;
  updateOutputStructures();
}

void HMM::setStorePerPairMap(bool storePerPairMAP)
{
  m_storePerPairMAP = storePerPairMAP;
  updateOutputStructures();
}

void HMM::setWritePerPairMap(bool writePerPairMAP)
{
  m_writePerPairMAP = writePerPairMAP;
  updateOutputStructures();
}

void HMM::setStorePerPairPosterior(bool storePerPairPosterior)
{
  m_storePerPairPosterior = storePerPairPosterior;
  updateOutputStructures();
}

void HMM::setStoreSumOfPosterior(bool storeSumOfPosterior)
{
  m_storeSumOfPosterior = storeSumOfPosterior;
  updateOutputStructures();
}
bool HMM::getStorePerPairPosteriorMean() const
{
  return m_storePerPairPosteriorMean;
}
bool HMM::getWritePerPairPosteriorMean() const
{
  return m_writePerPairPosteriorMean;
}
bool HMM::getStorePerPairMap() const
{
  return m_storePerPairMAP;
}
bool HMM::getWritePerPairMap() const
{
  return m_writePerPairMAP;
}
bool HMM::getStorePerPairPosterior() const
{
  return m_storePerPairPosterior;
}
bool HMM::getStoreSumOfPosterior() const
{
  return m_storeSumOfPosterior;
}
const DecodingParams& HMM::getDecodingParams() const
{
  return decodingParams;
}
