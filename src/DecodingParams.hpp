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

#include <limits>
#include <string>

enum class DecodingMode { sequenceFolded, arrayFolded, sequence, array };

enum class DecodingModeOverall { sequence, array };

class DecodingParams
{

private:
  bool fastSmcInvokedWithProgramOptions = false;

  /**
   * Resolve skipCSFSdistance against compress. compress is shorthand for skipCSFSdistance = infinity, so it conflicts
   * with any other distance. NaN means the distance was not set: it becomes infinity with compress and 0 (CSFS at every
   * site) without. Running this more than once gives the same result, so parameters can be validated repeatedly.
   *
   * @param optionPrefix "--" when the parameters came from the command line, so that messages name the flags
   * @throws std::invalid_argument if compress is combined with useAncestral or with a finite skipCSFSdistance
   */
  void resolveCompressAndCsfsDistance(const std::string& optionPrefix);

  /**
   * Set decodingModeOverall, decodingMode, decodingSequence and foldData from decodingModeString and useAncestral.
   *
   * @throws std::invalid_argument if decodingModeString is not "sequence" or "array" (in any case)
   */
  void resolveDecodingMode();

public:
  std::string inFileRoot;
  std::string decodingQuantFile;
  std::string mapFile;
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

  /**
   * Validate the ASMC parameters and fill in the values derived from them (decoding mode, CSFS distance, jobs, output
   * root).
   *
   * @return true
   * @throws std::invalid_argument if the parameters are incompatible
   */
  bool processOptions();

  /**
   * Parse and validate ASMC command-line arguments.
   *
   * @return false, after printing a message, if the arguments are invalid
   */
  bool processCommandLineArgs(int argc, char* argv[]);

  /**
   * Parse and validate FastSMC command-line arguments.
   *
   * @return false, after printing a message, if the arguments are invalid
   */
  bool processCommandLineArgsFastSMC(int argc, char* argv[]);

  /**
   * Verify that the selected FastSMC parameters are compatible, and fill in the values derived from them.
   *
   * @return true
   * @throws std::invalid_argument if the parameters are incompatible
   */
  bool validateParamsFastSMC();

  /**
   * Print decoding properties that are currently active.
   * 
   * @return true; 
   */ 
  bool printDecodingParams();

  /**
   * Constructor requiring only an input file root with all other parameters set to sensible defaults.
   * Decoding quantities will be generated if they are not specified and do not exist in the input file directory.
   */
  DecodingParams();

  /**
   * Constructor for ASMC. Runs processOptions, so it throws std::invalid_argument for incompatible parameters.
   * _skipCSFSdistance defaults to NaN, meaning "not set": see resolveCompressAndCsfsDistance.
   */
  explicit DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile = "", std::string _outFileRoot = "",
                          int _jobs = 1, int _jobInd = 1, std::string _decodingModeString = "array",
                          bool _decodingSequence = false, bool _usingCSFS = true, bool _compress = false,
                          bool _useAncestral = false, float _skipCSFSdistance = std::numeric_limits<float>::quiet_NaN(),
                          bool _noBatches = false, bool _doPosteriorSums = false, bool _doPerPairPosteriorMean = false,
                          std::string _expectedCoalTimesFile = "", bool _withinOnly = false,
                          bool _doMajorMinorPosteriorSums = false, bool _doPerPairMAP = false,
                          std::string _mapFile = "");

  /**
   * Minimal constructor that sets defaults for FastSMC. It throws std::invalid_argument if _fastSMC is false.
   *
   * @param _inFileRoot the input file root
   * @param _decodingQuantFile the decoding quantities file
   * @param _outFileRoot the output file root
   * @param _fastSMC whether to run in FastSMC: if this is set to false an error will occur
   */
  DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile, std::string _outFileRoot, bool _fastSMC = true);


};

#endif
