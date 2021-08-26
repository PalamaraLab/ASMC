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

#ifndef DECODINGPARAMS_HPP
#define DECODINGPARAMS_HPP

#include <string>

#include <boost/program_options.hpp>

enum class DecodingMode { sequenceFolded, arrayFolded, sequence, array };

enum class DecodingModeOverall { sequence, array };

class DecodingParams
{

private:
  bool fastSmcInvokedWithProgramOptions = false;

public:
  std::string inFileRoot;
  std::string decodingQuantFile;
  std::string outFileRoot;
  int jobs = 1;
  int jobInd = 1;
  std::string decodingModeString = "array";
  DecodingModeOverall decodingModeOverall;
  DecodingMode decodingMode;
  bool decodingSequence = false;
  bool foldData = false;
  bool usingCSFS = false;
  bool compress = false;
  bool useAncestral = false;
  float skipCSFSdistance{};
  bool noBatches = false;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // New params from FastSMC that were not originally in ASMC

  int batchSize = 64;
  int recallThreshold = 3;

  float skip = 0.f;
  int gap = 1;
  int max_seeds = 0;
  float min_maf = 0;
  float min_m = 1;
  bool hashing = false;
  bool hashingOnly = false;
  bool FastSMC = false;
  bool BIN_OUT = false;
  bool useKnownSeed = false;

  /// Whether to write IBD segment length (in centimorgans)
  bool outputIbdSegmentLength = false;

  // Used by FastSCM itself
  int hashingWordSize = 64;
  int constReadAhead = 10;
  bool haploid = true;

  int time = 100; // state threshold for IBD detection


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // main tasks
  bool noConditionalAgeEstimates = false;
  bool doPosteriorSums = false;
  bool doPerPairPosteriorMean = false; // output posterior mean for each pair
  bool doPerPairMAP = false;           // output MAP for each pair
  std::string expectedCoalTimesFile;   // expected coalescence times within each interval
  bool withinOnly = false;             // only compute decoding within individuals
  bool doMajorMinorPosteriorSums = false;

  bool processOptions();
  bool processCommandLineArgs(int argc, char* argv[]);
  bool processCommandLineArgsFastSMC(int argc, char* argv[]);

  /**
   * Verify that the selected parameters are compatible. Incompatible options will cause FastSMC to exit with a message
   * explaining the incompatibility.
   *
   * @return true if the parameters are compatible
   */
  bool validateParamsFastSMC();

  /**
   * Constructor requiring only an input file root with all other parameters set to sensible defaults.
   * Decoding quantities will be generated if they are not specified and do not exist in the input file directory.
   */
  DecodingParams();
  explicit DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile = "", std::string _outFileRoot = "",
                          int _jobs = 1, int _jobInd = 1, std::string _decodingModeString = "array",
                          bool _decodingSequence = false, bool _usingCSFS = true, bool _compress = false,
                          bool _useAncestral = false, float _skipCSFSdistance = 0.f, bool _noBatches = false,
                          bool _doPosteriorSums = false, bool _doPerPairPosteriorMean = false,
                          std::string _expectedCoalTimesFile = "", bool _withinOnly = false,
                          bool _doMajorMinorPosteriorSums = false, bool _doPerPairMAP = false);

  /**
   * Minimal constructor that sets defaults for FastSMC. An er will occur if you try to use this constructor for
   * FastSMC == false.
   *
   * @param _inFileRoot the input file root
   * @param _decodingQuantFile the decoding quantities file
   * @param _outFileRoot the output file root
   * @param _fastSMC whether to run in FastSMC: if this is set to false an error will occur
   */
  DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile, std::string _outFileRoot, bool _fastSMC = true);
};

#endif
