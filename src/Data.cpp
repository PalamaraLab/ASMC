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

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
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

  // Only FastSMC's hashing step divides samples between jobs; see mJobbing
  mJobbing = FastSMC && params.hashing;
  mJobInd = jobID;
  mJobs = jobs;

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

  readSamplesList(inFileRoot);

  for (auto i = 0ul; i < famAndIndNameList.size(); i++) {
    individuals.emplace_back(sites);
  }

  if (FastSMC) {
    vector<pair<unsigned long, double>> geneticMap = readMapFastSMC(inFileRoot, params.mapFile);
    readHaps(inFileRoot, foldToMinorAlleles, geneticMap);
  } else {
    readHaps(inFileRoot, foldToMinorAlleles);
    readMap(inFileRoot, params.mapFile);
  }
}

Data::Data(int numSites, unsigned long numDiploidSamples,
           std::vector<int> physPos,
           std::vector<float> genPos,
           const int8_t* haplotypes,
           bool usingCSFS, bool foldData_)
{
  sites = numSites;
  decodingUsesCSFS = usingCSFS;
  foldToMinorAlleles = foldData_;
  mJobbing = false;

  physicalPositions = std::move(physPos);
  geneticPositions = std::move(genPos);
  computeRecombinationRates();

  // Count the 1s at each site, folding to the minor allele if requested. haplotypes is laid out [site * numHaps + hap].
  const auto numHaps = static_cast<int>(2ul * numDiploidSamples);
  totalSamplesCount = vector<int>(sites, numHaps);
  derivedAlleleCounts = vector<int>(sites, 0);
  siteWasFlippedDuringFolding = vector<bool>(sites, false);
  for (int site = 0; site < sites; site++) {
    const int8_t* siteHaplotypes = haplotypes + static_cast<std::size_t>(site) * numHaps;
    const auto count = static_cast<int>(std::count(siteHaplotypes, siteHaplotypes + numHaps, int8_t{1}));
    const bool flip = foldToMinorAlleles && count > numHaps - count;
    siteWasFlippedDuringFolding[site] = flip;
    derivedAlleleCounts[site] = flip ? numHaps - count : count;
  }

  setGenotypes(haplotypes, numHaps);
}

// Construct from precomputed allele counts only (no genotypes).
// Emission tables are computed from these counts during HMM construction.
// Genotypes are loaded later via loadSubsetGenotypes() for each decode batch.
Data::Data(std::vector<int> derivedCounts, std::vector<int> totalCounts,
           std::vector<int> physPos, std::vector<float> genPos,
           bool usingCSFS, bool foldData_)
{
  sites = static_cast<int>(derivedCounts.size());
  sampleSize = 0;
  haploidSampleSize = 0;
  decodingUsesCSFS = usingCSFS;
  foldToMinorAlleles = foldData_;
  mJobbing = false;

  physicalPositions = std::move(physPos);
  geneticPositions = std::move(genPos);
  derivedAlleleCounts = std::move(derivedCounts);
  totalSamplesCount = std::move(totalCounts);
  computeRecombinationRates();

  if (totalSamplesCount.size() != derivedAlleleCounts.size()) {
    throw std::invalid_argument(
        fmt::format("Expected {} total sample counts, one per site, but got {}", sites, totalSamplesCount.size()));
  }

  // Fold to the minor allele if requested, as the in-memory genotype constructor does. setGenotypes flips genotypes
  // loaded later at the same sites.
  siteWasFlippedDuringFolding = std::vector<bool>(sites, false);
  for (int site = 0; site < sites; site++) {
    const int derived = derivedAlleleCounts[site];
    const int total = totalSamplesCount[site];
    if (derived < 0 || derived > total) {
      throw std::invalid_argument(fmt::format(
          "Site {} has {} derived alleles out of {}; the count must be between 0 and the total", site, derived, total));
    }
    if (foldToMinorAlleles && derived > total - derived) {
      siteWasFlippedDuringFolding[site] = true;
      derivedAlleleCounts[site] = total - derived;
    }
  }
}

// Load a subset of haplotypes for pair decoding. Called between decode batches.
// Does NOT recompute emissions — allele counts are already set from the constructor.
void Data::loadSubsetGenotypes(const int8_t* haplotypes, int numSubsetHaps)
{
  setGenotypes(haplotypes, numSubsetHaps);
}

void Data::setGenotypes(const int8_t* haplotypes, const int numHaps)
{
  const auto numDiploid = static_cast<unsigned long>(numHaps / 2);
  sampleSize = numDiploid;
  haploidSampleSize = static_cast<unsigned long>(numHaps);

  // Placeholder sample names: decoding pairs looks them up in IIDList
  FamIDList.clear();
  IIDList.clear();
  famAndIndNameList.clear();
  for (unsigned long i = 0; i < numDiploid; i++) {
    const std::string name = fmt::format("sample_{}", i);
    FamIDList.push_back(name);
    IIDList.push_back(name);
    famAndIndNameList.push_back(name + "\t" + name);
  }

  individuals.clear();
  individuals.reserve(numDiploid);
  for (unsigned long i = 0; i < numDiploid; i++) {
    individuals.emplace_back(sites);
  }

  // A 1 is the derived allele, so at a site flipped during folding it is the major allele
  for (int site = 0; site < sites; site++) {
    for (int hap = 0; hap < numHaps; hap++) {
      const int8_t genotype = haplotypes[static_cast<std::size_t>(site) * numHaps + hap];
      if (genotype != 0 && genotype != 1) {
        throw std::invalid_argument(
            fmt::format("Genotype at site {}, haplotype {} is {}, but must be 0 or 1", site, hap, genotype));
      }
      individuals[hap / 2].setGenotype((hap % 2) + 1, site, (genotype == 1) != siteWasFlippedDuringFolding[site]);
    }
  }
}

void Data::computeRecombinationRates()
{
  // Consecutive sites can have identical physical and genetic positions. In this case the recombination rate is zero.
  recRateAtMarker = std::vector<float>(sites, 0.f);
  for (auto i = 1; i < sites; ++i) {
    const float genDistFromPrevious = geneticPositions[i] - geneticPositions[i - 1];
    const float physDistFromPrevious = static_cast<float>(physicalPositions[i] - physicalPositions[i - 1]);

    if (genDistFromPrevious < std::numeric_limits<float>::epsilon()) {
      recRateAtMarker.at(i) = 0.f;
    } else {
      recRateAtMarker.at(i) = genDistFromPrevious / physDistFromPrevious;
    }
  }

  // Re-use the first rec rate in place zero
  if (sites > 1) {
    recRateAtMarker.at(0) = recRateAtMarker.at(1);
  }
}

void Data::clearGenotypes()
{
  individuals.clear();
  sampleSize = 0;
  haploidSampleSize = 0;
}

// *** read genetic map
vector<pair<unsigned long, double>> Data::readMapFastSMC(const string& inFileRoot, const string& mapFile)
{
  auto map = asmc::GeneticMap(findMapFile(inFileRoot, mapFile));

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
  if (numberOfSuccesses < 0 || numberOfSuccesses > populationSize || sampleSize > populationSize) {
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

void Data::readMap(const std::string& inFileRoot, const std::string& mapFile)
{
  const std::string checkedMapFile = findMapFile(inFileRoot, mapFile);
  auto map = asmc::PlinkMap(checkedMapFile);

  if (map.getNumSites() != static_cast<unsigned long>(sites)) {
    throw std::runtime_error(
        fmt::format("ERROR. Expected map {} to contain {} sites, but found {}\n", checkedMapFile, sites, map.getNumSites()));
  }

  if (map.getGeneticPositions().empty()) {
    throw std::runtime_error(fmt::format("ERROR. Expected map {} to contain a column of genetic positions\n", checkedMapFile));
  }

  SNP_IDs = map.getSnpIds();
  geneticPositions.reserve(sites);
  physicalPositions.reserve(sites);

  for (auto i = 0ul; i < map.getNumSites(); ++i) {
    geneticPositions.emplace_back(static_cast<float>(map.getGeneticPositions().at(i)) / 100.f);
    physicalPositions.emplace_back(static_cast<int>(map.getPhysicalPositions().at(i)));
  }

  computeRecombinationRates();
}

void Data::readSamplesList(const string& inFileRoot)
{
  string line;
  FileUtils::AutoGzIfstream bufferedReader;
  bufferedReader.openOrExit(findSamplesFile(inFileRoot));

  unsigned long linesProcessed = 0u; // number of samples read so far
  while (getline(bufferedReader, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    if (!isSampleLine(splitStr)) {
      continue;
    }
    if (isSampleInJob(linesProcessed)) {
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

bool Data::isSampleInJob(const unsigned long sampleIndex) const
{
  if (!mJobbing) {
    return true;
  }

  // Window w covers samples [(w - 1) * windowSize / 2, w * windowSize / 2). The last job also takes every sample from
  // the start of window w_j onwards, so that samples beyond the last full window are included.
  const auto windowStart = [this](const unsigned window) { return (uint)((window - 1) * windowSize) / 2; };
  const auto windowEnd = [this](const unsigned window) { return (uint)(window * windowSize) / 2; };
  return (sampleIndex >= windowStart(w_i) && sampleIndex < windowEnd(w_i)) ||
         (sampleIndex >= windowStart(w_j) && sampleIndex < windowEnd(w_j)) ||
         (mJobs == mJobInd && sampleIndex >= windowStart(w_j));
}

int Data::countHapLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  hapsBr.openOrExit(findHapsFile(inFileRoot));
  string line;
  int pos = 0;
  while (getline(hapsBr, line)) {
    // Blank lines are not sites: readHaps skips them too
    if (line.find_first_not_of(" \t\r") != string::npos) {
      pos++;
    }
  }
  return pos;
}

int Data::countSamplesLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream sampleBr;
  sampleBr.openOrExit(findSamplesFile(inFileRoot));
  string line;
  int cpt = 0;
  while (getline(sampleBr, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    if (isSampleLine(splitStr)) {
      cpt++;
    }
  }
  sampleBr.close();
  return cpt;
}

std::string Data::findSamplesFile(const std::string& inFileRoot)
{
  for (const char* extension : {".samples", ".sample"}) {
    if (FileUtils::fileExists(inFileRoot + extension)) {
      return inFileRoot + extension;
    }
  }
  throw std::runtime_error(
      fmt::format("ERROR. Could not find sample file in {}.samples or {}.sample", inFileRoot, inFileRoot));
}

std::string Data::findHapsFile(const std::string& inFileRoot)
{
  for (const char* extension : {".hap.gz", ".hap", ".haps.gz", ".haps"}) {
    if (FileUtils::fileExists(inFileRoot + extension)) {
      return inFileRoot + extension;
    }
  }
  throw std::runtime_error(
      fmt::format("ERROR. Could not find hap file in {0}.hap.gz, {0}.hap, {0}.haps.gz or {0}.haps", inFileRoot));
}

std::string Data::findMapFile(const std::string& inFileRoot, const std::string& mapFile)
{
  // A map file the user named must exist: falling back to another map would silently use the wrong one
  if (!mapFile.empty()) {
    if (!fs::is_regular_file(mapFile)) {
      throw std::runtime_error(fmt::format("ERROR. Map file {} does not exist", mapFile));
    }
    return mapFile;
  }
  for (const char* extension : {".map.gz", ".map"}) {
    if (fs::is_regular_file(inFileRoot + extension)) {
      return inFileRoot + extension;
    }
  }
  throw std::runtime_error(fmt::format("ERROR. Could not find map file in {0}.map.gz or {0}.map", inFileRoot));
}

bool Data::isSampleLine(const std::vector<std::string>& fields)
{
  if (fields.empty()) {
    return false;
  }
  const bool isHeader = fields.size() >= 3 && ((fields[0] == "ID_1" && fields[1] == "ID_2" && fields[2] == "missing") ||
                                               (fields[0] == "0" && fields[1] == "0" && fields[2] == "0"));
  if (isHeader) {
    return false;
  }
  if (fields.size() < 2) {
    throw std::runtime_error(
        fmt::format("ERROR. Samples file line \"{}\" should have a family ID and an individual ID", fields.front()));
  }
  return true;
}

int Data::parseChromosomeNumber(const std::string& chromosome)
{
  // Ignore anything after a ':' (e.g. "1:12345") and a "chr" prefix in any case (e.g. "chr1")
  std::string name = chromosome.substr(0, chromosome.find(':'));
  const auto toLower = [](const char c) { return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); };
  if (name.size() > 3 && std::equal(name.begin(), name.begin() + 3, "chr",
                                    [&toLower](const char a, const char b) { return toLower(a) == b; })) {
    name = name.substr(3);
  }

  int number = 0;
  const auto [end, error] = std::from_chars(name.data(), name.data() + name.size(), number);
  const bool isWholeNumber = error == std::errc() && end == name.data() + name.size();
  return isWholeNumber && number > 0 && number <= 1260 ? number : 0;
}

void Data::readHaps(string inFileRoot, bool foldToMinorAlleles)
{
  FileUtils::AutoGzIfstream hapsBr;
  hapsBr.openOrExit(findHapsFile(inFileRoot));

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
      throw std::runtime_error(
          fmt::format("ERROR: haps line has wrong length. Length is {}, but should be 4 * {} = {}.\n\tHaps line is: {}",
                      line.length(), famAndIndNameList.size(), 4 * famAndIndNameList.size(), line));
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
        throw std::runtime_error(fmt::format("ERROR: hap is not '0' or '1' at site {}", pos + 1));
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
  requireAllSitesRead(pos);

  cout << "Read data for " << sampleSize * 2 << " haploid samples and " << pos << " markers, " << monomorphic
       << " of which are monomorphic. This job will focus on " << FamIDList.size() * 2 << " haploid samples." << endl;
}

void Data::requireAllSitesRead(const unsigned long sitesRead) const
{
  // Reading stops at the first record that does not parse (e.g. a non-numeric position), which would otherwise leave
  // the remaining sites silently empty
  if (sitesRead != static_cast<unsigned long>(sites)) {
    throw std::runtime_error(fmt::format("ERROR: could not read site {} of {} in the haps file: check that it has a "
                                         "chromosome, SNP ID, numeric position and two alleles",
                                         sitesRead + 1, sites));
  }
}

void Data::readHaps(const string& inFileRoot, bool foldToMinorAlleles,
                    vector<pair<unsigned long int, double>>& genetic_map)
{
  FileUtils::AutoGzIfstream hapsBr;
  hapsBr.openOrExit(findHapsFile(inFileRoot));

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
      throw std::runtime_error(
          fmt::format("ERROR: haps line has wrong length. Length is {}, but should be 4 * {} = {}.\n\tHaps line is: {}",
                      line.length(), famAndIndNameList.size(), 4 * famAndIndNameList.size(), line));
    }

    // it is required that the rows in the haps file are sorted by increasing physical position
    if (bp > largest_bp) {
      largest_bp = bp;
    } else {
      throw std::runtime_error(
          fmt::format("ERROR: rows in haps data file must be ordered by increasing physical position, but two "
                      "consecutive values were {} and {}",
                      largest_bp, bp));
    }

    if (cur_g == 0) {
      chrNumber = parseChromosomeNumber(chr);
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
      if (isSampleInJob(d)) {
        Individual& ind = individuals[cpt];
        cpt++;
        uint hap1 = 2 * d;
        uint hap2 = 2 * d + 1;

        if (line[2 * hap1 + 1] == '1') {
          ind.setGenotype(hap1 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap1 + 1] == '0') {
          ind.setGenotype(hap1 % 2 + 1, pos, !minorAlleleValue);
        } else {
          throw std::runtime_error(fmt::format("ERROR: hap is not '0' or '1' at site {}", pos + 1));
        }

        if (line[2 * hap2 + 1] == '1') {
          ind.setGenotype(hap2 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap2 + 1] == '0') {
          ind.setGenotype(hap2 % 2 + 1, pos, !minorAlleleValue);
        } else {
          throw std::runtime_error(fmt::format("ERROR: hap is not '0' or '1' at site {}", pos + 1));
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
  requireAllSitesRead(pos);

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
      std::string id = SNP_IDs.size() > i ? fmt::format("ID {}", SNP_IDs.at(i)) : fmt::format("numerical ID {}", i);
      throw std::runtime_error(
          fmt::format("ERROR. SNP with {} has {} non-missing individuals, but the CSFS requires at least {}", id,
                      totalSamples, numCsfsSamples));
    }

    const int ancestralAlleles = totalSamples - derivedAlleles;
    if (foldToMinorAlleles && derivedAlleles > ancestralAlleles) {
      throw std::runtime_error("Minor alleles has frequency > 50%. Data is supposed to be folded.");
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
