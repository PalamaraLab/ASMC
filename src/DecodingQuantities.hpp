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

#ifndef DECODINGQUANTITIES_HPP
#define DECODINGQUANTITIES_HPP

#include <string>
#include <unordered_map>
#include <vector>

enum class DataType {
  TransitionType,
  States,
  CSFSSamples,
  TimeVector,
  SizeVector,
  Discretization,
  ExpectedTimes,
  CSFS,
  FoldedCSFS,
  ClassicEmission,
  AscertainedCSFS,
  FoldedAscertainedCSFS,
  CompressedAscertainedEmission,
  initialStateProb,
  ColumnRatios,
  RowRatios,
  Uvectors,
  Bvectors,
  Dvectors,
  HomozygousEmissions,
  None
};

class DecodingQuantities
{

public:
  unsigned int states = 0u;
  int CSFSSamples = 0;
  std::vector<float> initialStateProb;
  std::vector<float> expectedTimes;
  std::vector<float> discretization;
  std::vector<float> timeVector;
  std::vector<float> columnRatios;
  std::vector<std::vector<float>> classicEmissionTable;
  std::vector<std::vector<float>> compressedEmissionTable;
  std::unordered_map<float, std::vector<float>> Dvectors;
  std::unordered_map<float, std::vector<float>> Bvectors;
  std::unordered_map<float, std::vector<float>> Uvectors;
  std::unordered_map<float, std::vector<float>> rowRatioVectors;
  std::unordered_map<int, std::vector<float>> homozygousEmissionMap;
  std::vector<std::vector<std::vector<float>>> CSFSmap;
  std::vector<std::vector<std::vector<float>>> foldedCSFSmap;
  std::vector<std::vector<std::vector<float>>> ascertainedCSFSmap;
  std::vector<std::vector<std::vector<float>>> foldedAscertainedCSFSmap;

  explicit DecodingQuantities(const std::string& fileName);

private:
  // implemented, but need to update other code
  // void createFromBinary(const char *fileName);
  void createFromGzippedText(const std::string& fileName);

  /**
   * Validate that an appropriate decoding quantities file has been provided. This is achieved by:
   *
   * 1. Verifying the file exists
   * 2. Verifying the first line of the file contains exactly "TransitionType"
   *
   * @param fileName the name of the provided decoding quantities file
   */
  void validateDecodingQuantitiesFile(const std::string& fileName);

};

#endif
