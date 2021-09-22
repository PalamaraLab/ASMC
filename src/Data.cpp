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

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Types.hpp"

#include "Data.hpp"
#include <boost/math/distributions/hypergeometric.hpp>

#include "GeneticMap.hpp"
#include "PlinkMap.hpp"

using namespace std;
namespace fs = std::filesystem;

Data::Data(const DecodingParams& params)
{
  // Copy the DecodingParams that we need
  const std::string inFileRoot = params.inFileRoot;
  const int jobID = params.jobInd;
  const int jobs = params.jobs;
  foldToMinorAlleles = params.foldData;
  decodingUsesCSFS = params.usingCSFS;
  const bool FastSMC = params.FastSMC;

  // Determine if there is jobbing based on whether the jobID and jobs are at their default values
  mJobbing = (jobID != -1) && (jobs != -1);

  sites = countHapLines(inFileRoot);
  sampleSize = countSamplesLines(inFileRoot);
  haploidSampleSize = sampleSize * 2ul;

  siteWasFlippedDuringFolding = std::vector<bool>(sites, false);

  if (params.useKnownSeed) {
    std::srand(1234u);
  } else {
    std::random_device rd;
    std::srand(rd());
  }

  if (mJobbing) {
    // the window size is the length of a square, in terms of #ind
    const auto floatSampleSize = static_cast<double>(sampleSize);
    windowSize = ceil(sqrt((2. * pow(floatSampleSize, 2) - floatSampleSize) * 2. / jobs));
    if (windowSize % 2 != 0) {
      windowSize++;
    }

    w_i = 1; // window for individual i
    int cpt_job = 1;
    int cpt_tot_job = 1;
    while (cpt_tot_job < jobID) {
      w_i++;
      cpt_job = cpt_job + 2;
      cpt_tot_job = cpt_tot_job + cpt_job;
    }
    w_j = ceil((float)(cpt_job - (cpt_tot_job - jobID)) / 2); // window for individual j
    is_j_above_diag = (cpt_job - (cpt_tot_job - jobID)) % 2 == 1;
  }

  readSamplesList(inFileRoot, jobID, jobs);

  for (auto i = 0ul; i < famAndIndNameList.size(); i++) {
    individuals.emplace_back(sites);
  }

  if (FastSMC) {
    vector<pair<unsigned long, double>> geneticMap = readMapFastSMC(inFileRoot);
    readHaps(inFileRoot, foldToMinorAlleles, jobID, jobs, geneticMap);
  } else {
    readHaps(inFileRoot, foldToMinorAlleles);
    readMap(inFileRoot);
  }
}

// *** read genetic map
vector<pair<unsigned long, double>> Data::readMapFastSMC(const string& inFileRoot)
{
  std::string mapFile;
  if (mapFile = inFileRoot + ".map.gz"; fs::is_regular_file(mapFile)) {
  } else if (mapFile = inFileRoot + ".map"; fs::is_regular_file(mapFile)) {
  } else {
    fmt::print(cerr, "ERROR. Could not find map file in {}.map.gz or {}.map\n", inFileRoot, inFileRoot);
    exit(1);
  }

  auto map = asmc::GeneticMap(mapFile);

  const std::vector<unsigned long>& physicalPositions = map.getPhysicalPositions();
  const std::vector<double>& geneticPositions = map.getGeneticPositions();

  std::vector<std::pair<unsigned long, double>> geneticMap;
  geneticMap.reserve(physicalPositions.size());

  for (auto i = 0ul; i < physicalPositions.size(); ++i) {
    geneticMap.emplace_back(physicalPositions.at(i), geneticPositions.at(i));
  }

  return geneticMap;
}

// unoptimized sampling of hypergeometric
int Data::sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize)
{
  if (numberOfSuccesses < 0 || numberOfSuccesses > populationSize) {
    return -1;
  }
  vector<unsigned short> samplingVector;
  samplingVector = vector<unsigned short>(populationSize, 0);
  for (int i = 0; i < numberOfSuccesses; i++) {
    samplingVector[i] = 1;
  }
  std::shuffle(samplingVector.begin(), samplingVector.end(), std::mt19937(std::rand()));
  int ret = 0;
  for (int i = 0; i < sampleSize; i++) {
    ret += samplingVector[i];
  }
  return ret;
}

void Data::readMap(const std::string& inFileRoot)
{
  std::string mapFile;
  if (mapFile = inFileRoot + ".map.gz"; fs::is_regular_file(mapFile)) {
  } else if (mapFile = inFileRoot + ".map"; fs::is_regular_file(mapFile)) {
  } else {
    fmt::print(cerr, "ERROR. Could not find map file in {}.map.gz or {}.map\n", inFileRoot, inFileRoot);
    exit(1);
  }
  auto map = asmc::PlinkMap(mapFile);

  if (map.getNumSites() != static_cast<unsigned long>(sites)) {
    throw std::runtime_error(
        fmt::format("ERROR. Expected map {} to contain {} sites, but found {}\n", mapFile, sites, map.getNumSites()));
  }

  if (map.getGeneticPositions().empty()) {
    throw std::runtime_error(fmt::format("ERROR. Expected map {} to contain a column of genetic positions\n", mapFile));
  }

  SNP_IDs = map.getSnpIds();
  geneticPositions.reserve(sites);
  physicalPositions.reserve(sites);
  recRateAtMarker = std::vector<float>(sites);

  for (auto i = 0ul; i < map.getNumSites(); ++i) {
    geneticPositions.emplace_back(static_cast<float>(map.getGeneticPositions().at(i)) / 100.f);
    physicalPositions.emplace_back(static_cast<int>(map.getPhysicalPositions().at(i)));
  }

  // When calculating the recombination rates, note that it is possible to have consecutive entries in the map with
  // identical physical & genetic positions. In this case we set the recombination rate to zero.
  for (auto i = 1ul; i < map.getNumSites(); ++i) {
    const float genDistFromPrevious = geneticPositions[i] - geneticPositions[i - 1];
    const float physDistFromPrevious = static_cast<float>(physicalPositions[i] - physicalPositions[i - 1]);

    if(genDistFromPrevious < std::numeric_limits<float>::epsilon()) {
      recRateAtMarker.at(i) = 0.f;
    } else {
      recRateAtMarker.at(i) = genDistFromPrevious / physDistFromPrevious;
    }
  }

  // Re-use the first rec rate in place zero
  recRateAtMarker.at(0) = recRateAtMarker.at(1);
}

void Data::readSamplesList(const string& inFileRoot, int jobID, int jobs)
{
  string line;
  // Read samples file
  FileUtils::AutoGzIfstream bufferedReader;
  if (FileUtils::fileExists(inFileRoot + ".samples")) {
    bufferedReader.openOrExit(inFileRoot + ".samples");
  } else if (FileUtils::fileExists(inFileRoot + ".sample")) {
    bufferedReader.openOrExit(inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + inFileRoot + ".sample or " + inFileRoot + ".samples" << endl;
    exit(1);
  }

  unsigned long linesProcessed = 0u; // number of current line being processed
  while (getline(bufferedReader, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing") ||
        (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }
    if (readSample(linesProcessed, jobID, jobs)) {
      string famId = splitStr[0];
      string IId = splitStr[1];
      FamIDList.push_back(famId);
      IIDList.push_back(IId);
      famAndIndNameList.push_back(famId + "\t" + IId);
    }
    linesProcessed++;
  }
  bufferedReader.close();
  cout << "Read " << linesProcessed << " samples." << endl;
}

bool Data::readSample(const unsigned linesProcessed, const int jobID, const int jobs)
{

  // If we are not doing jobbing we read all samples
  if (!mJobbing) {
    return true;
  } else {
    return (linesProcessed >= (uint)((w_i - 1) * windowSize) / 2 && linesProcessed < (uint)(w_i * windowSize) / 2) ||
           (linesProcessed >= (uint)((w_j - 1) * windowSize) / 2 && linesProcessed < (uint)(w_j * windowSize) / 2) ||
           (jobs == jobID && linesProcessed >= (uint)((w_j - 1) * windowSize) / 2);
  }
}

int Data::countHapLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".hap")) {
    hapsBr.openOrExit(inFileRoot + ".hap");
  } else if (FileUtils::fileExists(inFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".haps")) {
    hapsBr.openOrExit(inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".hap.gz, " + inFileRoot + ".hap, " +
                ".haps.gz, or " + inFileRoot + ".haps"
         << endl;
    exit(1);
  }
  string line;
  int pos = 0;
  while (getline(hapsBr, line)) {
    pos++;
  }
  return pos;
}

int Data::countSamplesLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream sampleBr;
  if (FileUtils::fileExists(inFileRoot + ".samples")) {
    sampleBr.openOrExit(inFileRoot + ".samples");
  } else if (FileUtils::fileExists(inFileRoot + ".sample")) {
    sampleBr.openOrExit(inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + inFileRoot + ".sample or " + inFileRoot + ".samples" << endl;
    exit(1);
  }
  string line;
  int cpt = 0;
  while (getline(sampleBr, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing") ||
        (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }

    cpt++;
  }
  sampleBr.close();
  return cpt;
}

void Data::readHaps(string inFileRoot, bool foldToMinorAlleles)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".hap")) {
    hapsBr.openOrExit(inFileRoot + ".hap");
  } else if (FileUtils::fileExists(inFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".haps")) {
    hapsBr.openOrExit(inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".hap.gz, " + inFileRoot + ".hap, " +
                ".haps.gz, or " + inFileRoot + ".haps"
         << endl;
    exit(1);
  }

  unsigned long pos = 0ul;
  unsigned long monomorphic = 0ul;
  totalSamplesCount = vector<int>(sites);
  derivedAlleleCounts = vector<int>(sites);
  string chr, snpID;
  unsigned long bp;
  string line;
  string alleleA;
  string alleleB;

  while (hapsBr >> chr >> snpID >> bp >> alleleA >> alleleB) {
    getline(hapsBr, line);
    int DAcount = 0;
    if (!(line.length() == 4 * sampleSize || line.length() == 4 * sampleSize + 1)) {
      cerr << "ERROR: haps line has wrong length. Length is " << line.length() << ", should be 4*"
           << famAndIndNameList.size() << " = " << 4 * famAndIndNameList.size() << "." << endl;
      cerr << "\thaps line is: " << line << endl;
      exit(1);
    }

    int totalSamples = static_cast<int>(2 * famAndIndNameList.size());
    totalSamplesCount[pos] = totalSamples;
    for (uint i = 0; i < 2 * famAndIndNameList.size(); i++) {
      if (line[2 * i + 1] == '1') {
        DAcount++;
      }
    }
    bool minorAlleleValue = foldToMinorAlleles ? (DAcount <= totalSamples - DAcount) : true;
    siteWasFlippedDuringFolding[pos] = !minorAlleleValue;
    for (uint i = 0; i < 2 * famAndIndNameList.size(); i++) {
      int indIndex = i / 2;
      Individual& ind = individuals[indIndex];
      if (line[2 * i + 1] == '1') {
        ind.setGenotype(i % 2 + 1, pos, minorAlleleValue);
      } else if (line[2 * i + 1] == '0') {
        ind.setGenotype(i % 2 + 1, pos, !minorAlleleValue);
      } else {
        cerr << "ERROR: hap is not '0' or '1'" << endl;
        exit(1);
      }
    }

    if (foldToMinorAlleles) {
      derivedAlleleCounts[pos] = std::min(DAcount, totalSamples - DAcount);
    } else {
      derivedAlleleCounts[pos] = DAcount;
    }
    if (DAcount == 0 || DAcount == totalSamples) {
      monomorphic++;
    }
    pos++;
  }

  hapsBr.close();

  cout << "Read data for " << sampleSize * 2 << " haploid samples and " << pos << " markers, " << monomorphic
       << " of which are monomorphic. This job will focus on " << FamIDList.size() * 2 << " haploid samples." << endl;
}

void Data::readHaps(string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
                    vector<pair<unsigned long int, double>>& genetic_map)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".hap")) {
    hapsBr.openOrExit(inFileRoot + ".hap");
  } else if (FileUtils::fileExists(inFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".haps")) {
    hapsBr.openOrExit(inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".hap.gz, " + inFileRoot + ".hap, " + ".haps.gz, or " +
                inFileRoot + ".haps"
         << endl;
    exit(1);
  }

  unsigned long pos = 0ul;
  unsigned long monomorphic = 0ul;
  totalSamplesCount = vector<int>(sites);
  derivedAlleleCounts = vector<int>(sites);
  string chr, snpID;
  unsigned long bp = 0ul;
  unsigned long largest_bp = bp;
  string line;
  string alleleA;
  string alleleB;

  unsigned int cur_g = 0;
  while (hapsBr >> chr >> snpID >> bp >> alleleA >> alleleB) {
    getline(hapsBr, line);
    int DAcount = 0;
    if (!(line.length() == 4 * sampleSize || line.length() == 4 * sampleSize + 1)) {
      fmt::print(stderr,
                 "ERROR: haps line has wrong length. Length is {}, but should be 4 * {} = {}.\n\tHaps line is: {}\n",
                 line.length(), famAndIndNameList.size(), 4 * famAndIndNameList.size(), line);
      exit(1);
    }

    // it is required that the rows in the haps file are sorted by increasing physical position
    if (bp > largest_bp) {
      largest_bp = bp;
    } else {
      fmt::print(stderr,
                 "ERROR: rows in haps data file must be ordered by increasing physical position, but two consecutive "
                 "values were {} and {}\n",
                 largest_bp, bp);
      exit(1);
    }

    if (cur_g == 0) {
      // splitting inputs
      std::string delimiter = ":";
      size_t chrpos = chr.find(delimiter);
      if (std::string::npos == chrpos) {
        chrNumber = std::stoi(chr);
      } else {
        chrNumber = std::stoi(chr.substr(0, chrpos));
      }
      // check if chrNumber is positive and smaller than 1260, which is the max nb of chromosomes (Ophioglossum)
      // if it is not, set a default value to 0
      if (chrNumber <= 0 || chrNumber > 1260) {
        chrNumber = 0;
      }
    }
    readGeneticMap(bp, genetic_map, cur_g, pos);
    int totalSamples = 2 * sampleSize;
    totalSamplesCount[pos] = totalSamples;
    for (uint i = 0; i < 2 * sampleSize; i++) {
      if (line[2 * i + 1] == '1') {
        DAcount++;
      }
    }
    bool minorAlleleValue = foldToMinorAlleles ? (DAcount <= totalSamples - DAcount) : true;
    siteWasFlippedDuringFolding[pos] = !minorAlleleValue;
    // cout << (DAcount <= totalSamples - DAcount) << "\t" << siteWasFlippedDuringFolding[pos] << endl;
    uint cpt = 0;
    for (uint d = 0; d < sampleSize; d++) {
      if ((d >= (uint)((w_i - 1) * windowSize) / 2 && d < (uint)(w_i * windowSize) / 2) ||
          (d >= (uint)((w_j - 1) * windowSize) / 2 && d < (uint)(w_j * windowSize) / 2) ||
          (jobs == jobID && d >= (uint)((w_j - 1) * windowSize) / 2)) {
        Individual& ind = individuals[cpt];
        cpt++;
        uint hap1 = 2 * d;
        uint hap2 = 2 * d + 1;

        if (line[2 * hap1 + 1] == '1') {
          ind.setGenotype(hap1 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap1 + 1] == '0') {
          ind.setGenotype(hap1 % 2 + 1, pos, !minorAlleleValue);
        } else {
          cerr << "ERROR: hap is not '0' or '1'" << endl;
          exit(1);
        }

        if (line[2 * hap2 + 1] == '1') {
          ind.setGenotype(hap2 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap2 + 1] == '0') {
          ind.setGenotype(hap2 % 2 + 1, pos, !minorAlleleValue);
        } else {
          cerr << "ERROR: hap is not '0' or '1'" << endl;
          exit(1);
        }
      }
    }

    if (foldToMinorAlleles) {
      derivedAlleleCounts[pos] = std::min(DAcount, totalSamples - DAcount);
    } else {
      derivedAlleleCounts[pos] = DAcount;
    }

    if (DAcount == 0 || DAcount == totalSamples) {
      monomorphic++;
    }
    pos++;
  }

  hapsBr.close();

  cout << "Read data for " << sampleSize * 2 << " haploid samples and " << pos << " markers, " << monomorphic
       << " of which are monomorphic. This job will focus on " << FamIDList.size() * 2 << " haploid samples." << endl;
}

void Data::readGeneticMap(unsigned long int bp, vector<pair<unsigned long int, double>>& genetic_map,
                          unsigned int& cur_g, unsigned int pos)
{
  double cm;

  while (bp > genetic_map[cur_g].first && cur_g < genetic_map.size() - 1) {
    cur_g++;
  }

  if (bp >= genetic_map[cur_g].first) {
    // we found this exact marker, or we reached the end of the map
    cm = genetic_map[cur_g].second;
    addMarker(bp, cm, pos);
  } else if (cur_g == 0) {
    // if we haven't hit the map yet, store first map entry
    cm = genetic_map[cur_g].second;
    addMarker(bp, cm, pos);
  } else {
    // interpolate from previous marker
    cm = genetic_map[cur_g - 1].second + (bp - genetic_map[cur_g - 1].first) *
                                             (genetic_map[cur_g].second - genetic_map[cur_g - 1].second) /
                                             (genetic_map[cur_g].first - genetic_map[cur_g - 1].first);
    addMarker(bp, cm, pos);
  }
}

void Data::addMarker(unsigned long int physicalPos, double geneticPos, unsigned int pos)
{
  geneticPositions.push_back(geneticPos / 100.f);
  physicalPositions.push_back(physicalPos);
  physicalPositionsMap[physicalPos] = pos;

  if (pos > 0) {
    double genDistFromPrevious = geneticPositions[pos] - geneticPositions[pos - 1];
    unsigned long int physDistFromPrevious = physicalPositions[pos] - physicalPositions[pos - 1];
    float recRate = genDistFromPrevious / physDistFromPrevious;
    if (pos == 1) {
      // if it's first, add it again for marker 0. Using rate to next marker instead of previous marker
      recRateAtMarker.push_back(recRate);
    }
    recRateAtMarker.push_back(recRate);
  }
}

std::vector<std::vector<int>> Data::calculateUndistinguishedCounts(const int numCsfsSamples) const {

  // This matrix is of size <derivedAlleleCounts.size() x 3>
  std::vector<std::vector<int>> undistinguished(derivedAlleleCounts.size(), std::vector<int>(3));

  for (auto i = 0ul; i < derivedAlleleCounts.size(); ++i) {
    const int derivedAlleles = derivedAlleleCounts[i];
    const int totalSamples = totalSamplesCount[i];

    if (decodingUsesCSFS && numCsfsSamples > totalSamples) {
      std::string id = SNP_IDs.size() >= i ? fmt::format("ID {}", SNP_IDs.at(i)) : fmt::format("numerical ID {}", i);
      fmt::print(std::cerr, "ERROR. SNP with {} has {} non-missing individuals, but the CSFS requires at least {}\n",
                 id, totalSamples, numCsfsSamples);
      exit(1);
    }

    const int ancestralAlleles = totalSamples - derivedAlleles;
    if (foldToMinorAlleles && derivedAlleles > ancestralAlleles) {
      cerr << "Minor alleles has frequency > 50%. Data is supposed to be folded.\n";
      exit(1);
    }

    for (int distinguished = 0; distinguished < 3; distinguished++) {
      // hypergeometric with (derivedAlleles - distinguished) derived alleles, (samples - 2) samples
      int sample = sampleHypergeometric(totalSamples - 2, derivedAlleles - distinguished, numCsfsSamples - 2);
      if (foldToMinorAlleles && (sample + distinguished > numCsfsSamples / 2)) {
        sample = (numCsfsSamples - 2 - sample);
      }
      undistinguished.at(i).at(distinguished) = sample;
    }
  }

  return undistinguished;
}

const vector<std::string>& Data::getSnpIDs() const
{
  return SNP_IDs;
}
