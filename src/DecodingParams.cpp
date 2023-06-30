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

#include "DecodingParams.hpp"
#include "Types.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

DecodingParams::DecodingParams()
    : inFileRoot(""), decodingQuantFile(""), outFileRoot(inFileRoot), jobs(1), jobInd(1), decodingModeString("array"),
      decodingSequence(false), usingCSFS(true), compress(false), useAncestral(false), skipCSFSdistance(0.f),
      noBatches(false), doPosteriorSums(false), doPerPairPosteriorMean(false), expectedCoalTimesFile(""),
      withinOnly(false), doMajorMinorPosteriorSums(false), doPerPairMAP(false)
{
}

DecodingParams::DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile, std::string _outFileRoot,
                               int _jobs, int _jobInd, std::string _decodingModeString, bool _decodingSequence,
                               bool _usingCSFS, bool _compress, bool _useAncestral, float _skipCSFSdistance,
                               bool _noBatches, bool _doPosteriorSums, bool _doPerPairPosteriorMean,
                               std::string _expectedCoalTimesFile, bool _withinOnly, bool _doMajorMinorPosteriorSums,
                               bool _doPerPairMAP, std::string _mapFile)
    : inFileRoot(_inFileRoot), decodingQuantFile(_decodingQuantFile), mapFile(std::move(_mapFile)),
      outFileRoot(_outFileRoot), jobs(_jobs), jobInd(_jobInd), decodingModeString(_decodingModeString),
      decodingSequence(_decodingSequence), usingCSFS(_usingCSFS), compress(_compress), useAncestral(_useAncestral),
      skipCSFSdistance(_skipCSFSdistance), noBatches(_noBatches), doPosteriorSums(_doPosteriorSums),
      doPerPairPosteriorMean(_doPerPairPosteriorMean), expectedCoalTimesFile(_expectedCoalTimesFile),
      withinOnly(_withinOnly), doMajorMinorPosteriorSums(_doMajorMinorPosteriorSums), doPerPairMAP(_doPerPairMAP)
{
  if (!processOptions()) {
    throw std::exception();
  }
}

DecodingParams::DecodingParams(std::string _inFileRoot, std::string _decodingQuantFile, std::string _outFileRoot,
                               bool _fastSMC)
    : inFileRoot(std::move(_inFileRoot)), decodingQuantFile(std::move(_decodingQuantFile)),
      outFileRoot(std::move(_outFileRoot)), jobs(1), jobInd(1), decodingModeString("array"),
      decodingModeOverall(DecodingModeOverall::array), decodingMode(DecodingMode::arrayFolded), foldData(true),
      usingCSFS(true), batchSize(32), recallThreshold(3), min_m(1.5f), hashing(true), FastSMC(_fastSMC), BIN_OUT(false),
      outputIbdSegmentLength(true), time(50), noConditionalAgeEstimates(true), doPerPairPosteriorMean(true),
      doPerPairMAP(true)
{
  if (!FastSMC) {
    std::cerr << "This DecodingParams constructor sets sensible FastSMC defaults, and is only intended for use with"
                 "FastSMC. Please set the fastSMC parameter to true, or use a different constructor."
              << std::endl;
    exit(1);
  }

  validateParamsFastSMC();
}

bool DecodingParams::processCommandLineArgs(int argc, char* argv[])
{

  namespace po = boost::program_options;

  // clang-format off
  po::options_description options;
  options.add_options()
  ("inFileRoot", po::value<std::string>(&inFileRoot)->required(),
   "Prefix of hap|haps|hap.gz|haps.gz and sample|samples file")
  ("decodingQuantFile", po::value<std::string>(&decodingQuantFile),
   "Decoding quantities file")
  ("mapFile", po::value<std::string>(&mapFile),
"Map file: optional, if not in inFileRoot")
  ("jobs", po::value<int>(&jobs)->default_value(0),
   "Number of jobs being done in parallel")
  ("jobInd", po::value<int>(&jobInd)->default_value(0),
   "Job index (1..jobs)")
  ("outFileRoot", po::value<std::string>(&outFileRoot),
   "Output file for sum of posterior distribution over pairs (default: --inFileRoot argument)")
  ("mode", po::value<std::string>(&decodingModeString)->default_value("array"),
   "Decoding mode. Choose from {sequence, array}.")
  ("compress", po::bool_switch(&compress)->default_value(false),
   "Compress emission to binary (no CSFS)")
  // note: currently flipping major/minor when reading data. Instead, assume it's correctly coded to begin with
  ("useAncestral", po::bool_switch(&useAncestral)->default_value(false),
   "Assume ancestral alleles are coded as 1 in input (will assume 1 = minor otherwise)")
  // for debugging and pedagogical reasons.
 // ("noBatches", po::bool_switch(&noBatches)->default_value(false),
 //  "Do not decode in batches (for debugging, will remove)")
  ("skipCSFSdistance", po::value<float>(&skipCSFSdistance)->default_value(0.),
   "Genetic distance between two CSFS emissions")
  // main tasks
  ("majorMinorPosteriorSums", po::bool_switch(&doMajorMinorPosteriorSums)->default_value(false),
   "Output file for sum of posterior distribution over pairs, partitioned by major/minor alleles. O(sitesxstates) output")
  ("posteriorSums", po::bool_switch(&doPosteriorSums)->default_value(false),
   "Output file for sum of posterior distribution over pairs. O(sitesxstates) output");
  // ("perPairMAP", po::bool_switch(&doPerPairMAP)->default_value(false),
  //  "output per-pair MAP at each site. O(sitesxsamples^2) output")
  // ("perPairPosteriorMeans", po::value<string>(&expectedCoalTimesFile),
  //  "output per-pair posterior means at each site. O(sitesxsamples^2) output")
  // ("withinOnly", po::bool_switch(&withinOnly)->default_value(false),
  //  "Only decode pairs within unphased individuals");

  // clang-format on

  po::options_description visible("Options");
  visible.add(options);

  po::options_description all("All options");
  all.add(options);
  all.add_options()("bad-args", po::value<std::vector<std::string>>(), "bad args");

  po::positional_options_description positional_desc;
  positional_desc.add("bad-args", -1); // for error-checking command line

  po::variables_map vm;
  po::command_line_parser cmd_line(argc, argv);
  cmd_line.options(all);
  cmd_line.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
  cmd_line.positional(positional_desc);

  try {
    po::store(cmd_line.run(), vm);

    po::notify(vm); // throws an error if there are any problems

    if (vm.count("bad-args")) {
      std::cerr << "ERROR: Unknown options:";
      std::vector<std::string> bad_args = vm["bad-args"].as<std::vector<std::string>>();
      for (uint i = 0; i < bad_args.size(); i++)
        std::cerr << " " << bad_args[i];
      std::cerr << std::endl;
      return false;
    }
  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << options << std::endl;
    return false;
  }

  if (processOptions()) {
    if (!doPosteriorSums && !doPerPairMAP && !doPerPairPosteriorMean && !doMajorMinorPosteriorSums) {
      std::cerr << "ERROR: At least one of --posteriorSums, --majorMinorPosteriorSums, must be specified" << std::endl;
      return false;
    }
  } else
    return false;
  return true;
}

bool DecodingParams::processCommandLineArgsFastSMC(int argc, char* argv[])
{
  namespace po = boost::program_options;

  fastSmcInvokedWithProgramOptions = true;

  FastSMC = true;

  std::string decodingModeString;

  // clang-format off

  po::options_description options;
  options.add_options()
      ("inFileRoot", po::value<std::string>(&inFileRoot)->required(),
       "Prefix of hap|haps|hap.gz|haps.gz and sample|samples file.")
      ("outFileRoot", po::value<std::string>(&outFileRoot)->required(),
       "Output file for sum of posterior distribution over pairs.")
      ("decodingQuantFile", po::value<std::string>(&decodingQuantFile),
       "Decoding quantities file")
      ("mapFile", po::value<std::string>(&mapFile),
       "Map file: optional, if not in inFileRoot")
      ("mode", po::value<std::string>(&decodingModeString)->default_value("array"),
       "Decoding mode. Choose from {sequence, array}. [default = array]")
      ("time", po::value<int>(&time)->default_value(100),
       "Time threshold to define IBD in number of generations. [default = 100] ")
      ("jobs", po::value<int>(&jobs)->default_value(1),
       "number of jobs being done in parallel. [default = 1]")
      ("jobInd", po::value<int>(&jobInd)->default_value(1),
       "job index (1..jobs). [default = 1]")
      ("bin", po::bool_switch(&BIN_OUT)->default_value(false),
       "Binary output [default off]")
      ("batchSize", po::value<int>(&batchSize)->default_value(32),
       "Batch size [default = 32]")
      ("recall", po::value<int>(&recallThreshold)->default_value(3),
       "Recall level from 0 to 3 (higher value means higher recall). [default = 3]")

      // TASKS
      ("segmentLength", po::bool_switch(&outputIbdSegmentLength)->default_value(true),
       "Output length in centimorgans of each IBD segment. [default 1/on]")
      ("perPairMAP", po::bool_switch(&doPerPairMAP)->default_value(true),
       "Output per-pair MAP for each IBD segment. [default 1/on]")
      ("perPairPosteriorMeans", po::bool_switch(&doPerPairPosteriorMean)->default_value(true),
       "Output per-pair posterior means for each IBD segment. [default 1/on]")
      ("noConditionalAgeEstimates", po::bool_switch(&noConditionalAgeEstimates)->default_value(false),
       "Do not condition the age estimates on the TMRCA being between present time and t generations ago (where t is the time threshold). [default 0/off]")
      ("withinOnly", po::bool_switch(&withinOnly)->default_value(false),
       "Only decode pairs within unphased individuals. [default 0/off]")

      // TODO: currently flipping major/minor when reading data. Instead, assume it's correctly coded to begin with
      ("useAncestral", po::bool_switch(&useAncestral)->default_value(false),
       "Assume ancestral alleles are coded as 1 in input (will assume 1 = minor otherwise). [default 0/off]")
      ("compress", po::bool_switch(&compress)->default_value(false),
       "Compress emission to binary (no CSFS). [default 0/off]")

      // TODO: for debug. remove?
      ("noBatches", po::bool_switch(&noBatches)->default_value(false),
       "Do not decode in batches. [default 0/off]")
      ("skipCSFSdistance", po::value<float>(&skipCSFSdistance)->default_value(std::numeric_limits<float>::quiet_NaN()),
       "Genetic distance between two CSFS emissions")

      //hashing options
      ("hashing", po::bool_switch(&hashing)->default_value(true),
       "Use of hashing to pre-process IBD segments [default 1/on]")
      ("hashingOnly", po::bool_switch(&hashingOnly)->default_value(false),
       "Only perform GERMLINE2 hashing, not ASMC decoding [default 0/off]")
      ("min_m", po::value<float>(&min_m)->default_value(1.0),
       "Minimum match length (in cM). [default = 1.0]")
      ("skip", po::value<float>(&skip)->default_value(0.0),
       "Skip words with (seeds/samples) less than this value [default 0.0]")
      ("min_maf", po::value<float>(&min_maf)->default_value(0.0),
       "Minimum minor allele frequency [default 0.0]")
      ("gap", po::value<int>(&gap)->default_value(1),
       "Allowed gaps [default 1]")
      ("max_seeds", po::value<int>(&max_seeds)->default_value(0),
       "Dynamic hash seed cutoff [default 0/off]");

  // clang-format on

  po::options_description visible("Options");
  visible.add(options);

  po::options_description all("All options");
  all.add(options);
  all.add_options()("bad-args", po::value<std::vector<std::string>>(), "bad args");

  po::positional_options_description positional_desc;
  positional_desc.add("bad-args", -1); // for error-checking command line

  po::variables_map vm;
  po::command_line_parser cmd_line(argc, argv);
  cmd_line.options(all);
  cmd_line.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
  cmd_line.positional(positional_desc);

  try {
    po::store(cmd_line.run(), vm);

    po::notify(vm); // throws an error if there are any problems

    if (vm.count("bad-args")) {
      std::cerr << "ERROR: Unknown options:";
      std::vector<std::string> bad_args = vm["bad-args"].as<std::vector<std::string>>();
      for (uint i = 0; i < bad_args.size(); i++)
        std::cerr << " " << bad_args[i];
      std::cerr << std::endl;
      return false;
    }

  } catch (po::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
    std::cerr << options << std::endl;
    return false;
  }

  return validateParamsFastSMC();
}

bool DecodingParams::validateParamsFastSMC()
{
  const std::string del = fastSmcInvokedWithProgramOptions ? "--" : "";

  if (!FastSMC) {
    std::cerr
        << "Attempting to validate FastSMC parameters but FastSMC flag is false. Set DecodingParams::FastSMC to true?"
        << std::endl;
    exit(1);
  }

  if (hashingOnly) {
    if (!hashing) {
      hashing = true;
      fmt::print("Warning: {}hashingOnly set to true, but {}hashing set to false. Assuming {}hashing=true", del, del,
                 del);
    }
    if (doPerPairMAP || doPerPairPosteriorMean) {
      fmt::print("Warning: {}hashingOnly set to true which is incompatible with {}doPerPairMAP and "
                 "{}doPerPairPosteriorMean. Turning those options off.",
                 del, del, del);
      doPerPairMAP = false;
      doPerPairPosteriorMean = false;
    }
  }

  if (hashing) {
    if (withinOnly) {
      std::cerr << del << "hashing & " << del
                << "withinOnly cannot be used together. Please remove one of the two flags." << std::endl;
      exit(1);
    }
    if (time <= 0) {
      std::cerr << del << "time must be a positive integer." << std::endl;
      exit(1);
    }
  }

  if (batchSize == 0 || batchSize % 8 != 0) {
    std::cerr << del << "batchSize must be strictly positive and a multiple of 8." << std::endl;
    exit(1);
  }

  if (compress) {
    if (useAncestral) {
      std::cerr << del << "compress & " << del
                << "useAncestral cannot be used together. A compressed emission cannot use"
                << " ancestral allele information." << std::endl;
      exit(1);
    }
    if (!std::isnan(skipCSFSdistance)) {
      std::cerr << del << "compress & " << del << "skipCSFSdistance cannot be used together. " << del << "compress is a"
                << " shorthand for " << del << "skipCSFSdistance Infinity." << std::endl;
      exit(1);
    }
    skipCSFSdistance = std::numeric_limits<float>::infinity();
  } else {
    if (std::isnan(skipCSFSdistance)) {
      // default: use CSFS at all sites
      skipCSFSdistance = 0.f;
    }
  }

  if (skipCSFSdistance != std::numeric_limits<float>::infinity()) {
    usingCSFS = true;
  }

  boost::algorithm::to_lower(decodingModeString);
  if (decodingModeString == std::string("sequence")) {
    decodingSequence = true;
    if (useAncestral) {
      decodingMode = DecodingMode::sequence;
      foldData = false;
    } else {
      decodingMode = DecodingMode::sequenceFolded;
      foldData = true;
    }
  } else if (decodingModeString == std::string("array")) {
    decodingSequence = false;
    if (useAncestral) {
      decodingMode = DecodingMode::array;
      foldData = false;
    } else {
      decodingMode = DecodingMode::arrayFolded;
      foldData = true;
    }
  } else {
    std::cerr << "ERROR. Unknown decoding mode: " << decodingModeString << std::endl;
    exit(1);
  }

  if (decodingQuantFile.empty()) {
    std::cout << "Setting " << del << "decodingQuantFile to " << del << "inFileRoot + .decodingQuantities.bin"
              << std::endl;
    decodingQuantFile = inFileRoot + ".decodingQuantities.bin";
  }

  if ((jobs == 0) != (jobInd == 0)) {
    std::cerr << "ERROR: " << del << "jobs and " << del << "jobInd must either both be set or both be unset"
              << std::endl;
    exit(1);
  }

  if (jobs == 0) {
    jobs = 1;
    jobInd = 1;
  }

  if (jobInd <= 0 || jobInd > jobs || jobs <= 0) {
    std::cerr << "ERROR: " << del << "jobInd must be between 1 and " << del << "jobs inclusive" << std::endl;
    exit(1);
  }

  bool valid_job = false;
  int x = 1;
  int u = 1;
  int prev_u = u;
  for (int i = 0; i < 200; i++) {
    if (u == jobs) {
      valid_job = true;
      break;
    } else if (u > jobs) {
      break;
    }
    x = x + 2;
    prev_u = u;
    u = u + x;
  }

  if (!valid_job) {
    std::cerr << "ERROR: jobs value is incorrect. You should use either " << prev_u << " or " << u << std::endl;
    exit(1);
  }

  if (recallThreshold < 0 || recallThreshold > 3) {
    std::cerr << "ERROR: " << del << "recall must be between 0 and 3. " << std::endl;
    exit(1);
  }

  if (outFileRoot.empty()) {
    outFileRoot = inFileRoot;
    if (jobs > 0) {
      outFileRoot += "." + std::to_string(jobInd) + "-" + std::to_string(jobs);
    }
  }

  std::cout << std::boolalpha;
  std::cout << std::endl;
  std::cout << "---------------------------" << std::endl;
  std::cout << "        ASMC OPTIONS       " << std::endl;
  std::cout << "---------------------------" << std::endl;

  std::cout << "Input will have prefix : " << inFileRoot << std::endl;
  std::cout << "Decoding quantities file : " << decodingQuantFile << std::endl;
  std::cout << "Output will have prefix : " << outFileRoot << "." << jobInd << "." << jobs;

  if (hashing) {
    std::cout << ".FastSMC";
  } else {
    std::cout << ".asmc";
  }

  if (BIN_OUT) {
    std::cout << ".bibd" << std::endl;
  } else {
    std::cout << ".ibd.gz" << std::endl;
  }

  std::cout << "Binary output ? " << BIN_OUT << std::endl;
  std::cout << "Time threshold to define IBD in generations : " << time << std::endl;
  std::cout << "Use batches ? " << !noBatches << std::endl;

  if (!noBatches) {
    std::cout << "Batch size : " << batchSize << std::endl;
  }

  std::cout << "Running job " << jobInd << " of " << jobs << std::endl;
  std::cout << "Recall level " << recallThreshold << std::endl;
  std::cout << "skipCSFSdistance is " << skipCSFSdistance << std::endl;
  std::cout << "compress ? " << compress << std::endl;
  std::cout << "useAncestral ? " << useAncestral << std::endl;
  std::cout << "outputIbdSegmentLength ? " << outputIbdSegmentLength << std::endl;
  std::cout << "doPerPairPosteriorMean ? " << doPerPairPosteriorMean << std::endl;
  std::cout << "doPerPairMAP ? " << doPerPairMAP << std::endl;
  std::cout << "noConditionalAgeEstimates ? " << noConditionalAgeEstimates << std::endl;
  std::cout << "Use hashing as a preprocessing step ? " << hashing << std::endl;

  if (hashing) {
    std::cout << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "      hashing OPTIONS     " << std::endl;
    std::cout << "---------------------------" << std::endl;
    std::cout << "Minimum match length (in cM) : " << min_m << std::endl;
    std::cout << "Skipping words with (seeds/samples) less than " << skip << std::endl;
    std::cout << "Minimum minor allele frequency : " << min_maf << std::endl;
    std::cout << "Allowed gaps " << gap << std::endl;
    std::cout << "Dynamic hash seed cutoff : " << max_seeds << std::endl;
  }
  std::cout << std::noboolalpha;

  return true;
}

bool DecodingParams::processOptions()
{

  if (compress) {
    if (useAncestral) {
      std::cerr
          << "--compress & --useAncestral cannot be used together. A compressed emission cannot use ancestral allele "
             "information."
          << std::endl;
      exit(1);
    }
    if (!std::isnan(skipCSFSdistance)) {
      std::cerr << "--compress & --skipCSFSdistance cannot be used together. --compress is a shorthand for "
                   "--skipCSFSdistance Infinity."
                << std::endl;
      exit(1);
    }
    skipCSFSdistance = std::numeric_limits<float>::infinity();
  } else {
    if (std::isnan(skipCSFSdistance)) {
      // default: use CSFS at all sites
      skipCSFSdistance = 0.f;
    }
  }

  if (!expectedCoalTimesFile.empty()) {
    doPerPairPosteriorMean = true;
  }

  if (skipCSFSdistance != std::numeric_limits<float>::infinity()) {
    usingCSFS = true;
  }

  boost::algorithm::to_lower(decodingModeString);
  if (decodingModeString == std::string("sequence"))
    decodingModeOverall = DecodingModeOverall::sequence;
  else if (decodingModeString == std::string("array"))
    decodingModeOverall = DecodingModeOverall::array;
  else {
    std::cerr << "Decoding mode should be one of {sequence, array}";
    return false;
  }

  if (decodingModeOverall == DecodingModeOverall::sequence) {
    decodingSequence = true;
    if (useAncestral) {
      decodingMode = DecodingMode::sequence;
      foldData = false;
    } else {
      decodingMode = DecodingMode::sequenceFolded;
      foldData = true;
    }
  } else if (decodingModeOverall == DecodingModeOverall::array) {
    decodingSequence = false;
    if (useAncestral) {
      decodingMode = DecodingMode::array;
      foldData = false;
    } else {
      decodingMode = DecodingMode::arrayFolded;
      foldData = true;
    }
  } else {
    std::cerr << "ERROR. Unknown decoding mode: " << decodingModeString << std::endl;
    exit(1);
  }

  if (decodingQuantFile.empty()) {
    std::cout << "Setting --decodingQuantFile to --inFileRoot + .decodingQuantities.bin" << std::endl;
    decodingQuantFile = inFileRoot + ".decodingQuantities.bin";
  }

  if ((jobs == 0) != (jobInd == 0)) {
    std::cerr << "ERROR: --jobs and --jobInd must either both be set or both be unset" << std::endl;
    return false;
  }

  if (jobs == 0) {
    jobs = 1;
    jobInd = 1;
  }

  if (jobInd <= 0 || jobInd > jobs) {
    std::cerr << "ERROR: --jobInd must be between 1 and --jobs inclusive" << std::endl;
    return false;
  }
  if (outFileRoot.empty()) {
    outFileRoot = inFileRoot;
    if (jobs > 0) {
      outFileRoot += "." + std::to_string(jobInd) + "-" + std::to_string(jobs);
    }
  }
  return true;
}
