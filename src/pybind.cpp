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

#include "ASMC.hpp"
#include "BinaryDataReader.hpp"
#include "Data.hpp"
#include "DecodePairsReturnStruct.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "HMM.hpp"
#include "Individual.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <vector>

PYBIND11_MAKE_OPAQUE(std::vector<bool>)
PYBIND11_MAKE_OPAQUE(std::vector<float>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float>>)
PYBIND11_MAKE_OPAQUE(std::vector<Individual>)
PYBIND11_MAKE_OPAQUE(std::vector<PairObservations>)
PYBIND11_MAKE_OPAQUE(std::unordered_map<float, std::vector<float>>)
PYBIND11_MAKE_OPAQUE(std::unordered_map<int, std::vector<float>>)
PYBIND11_MAKE_OPAQUE(DecodingQuantities)
PYBIND11_MAKE_OPAQUE(DecodePairsReturnStruct)
PYBIND11_MAKE_OPAQUE(DecodingReturnValues)
PYBIND11_MAKE_OPAQUE(PairObservations)
PYBIND11_MAKE_OPAQUE(Individual)
PYBIND11_MAKE_OPAQUE(IbdPairDataLine)
PYBIND11_MAKE_OPAQUE(BinaryDataReader)
PYBIND11_MAKE_OPAQUE(Data)
PYBIND11_MAKE_OPAQUE(HMM)
PYBIND11_MAKE_OPAQUE(ASMC::FastSMC)

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MODULE(asmc_python_bindings, m)
{
  py::enum_<DecodingModeOverall>(m, "DecodingModeOverall", py::arithmetic())
      .value("sequence", DecodingModeOverall::sequence)
      .value("array", DecodingModeOverall::array);
  py::enum_<DecodingMode>(m, "DecodingMode", py::arithmetic())
      .value("sequenceFolded", DecodingMode::sequenceFolded)
      .value("arrayFolded", DecodingMode::arrayFolded)
      .value("sequence", DecodingMode::sequence)
      .value("array", DecodingMode::array);
  py::bind_vector<std::vector<bool>>(m, "VectorBool");
  py::bind_vector<std::vector<float>>(m, "VectorFloat");
  py::bind_vector<std::vector<Individual>>(m, "VectorIndividual");
  py::bind_vector<std::vector<uint>>(m, "VectorUInt");
  py::bind_vector<std::vector<PairObservations>>(m, "VectorPairObservations");
  py::bind_vector<std::vector<std::vector<float>>>(m, "Matrix");
  py::bind_map<std::unordered_map<float, std::vector<float>>>(m, "UMapFloatToVectorFloat");
  py::bind_map<std::unordered_map<int, std::vector<float>>>(m, "UMapIntToVectorFloat");
  py::class_<DecodingReturnValues>(m, "DecodingReturnValues")
      .def_readwrite("sumOverPairs", &DecodingReturnValues::sumOverPairs)
      .def_readwrite("sumOverPairs00", &DecodingReturnValues::sumOverPairs00)
      .def_readwrite("sumOverPairs01", &DecodingReturnValues::sumOverPairs01)
      .def_readwrite("sumOverPairs11", &DecodingReturnValues::sumOverPairs11)
      .def_readwrite("sites", &DecodingReturnValues::sites)
      .def_readwrite("states", &DecodingReturnValues::states)
      .def_readwrite("siteWasFlippedDuringFolding", &DecodingReturnValues::siteWasFlippedDuringFolding);
  py::class_<DecodePairsReturnStruct>(m, "DecodePairsReturnStruct")
      .def_readwrite("per_pair_indices", &DecodePairsReturnStruct::perPairIndices,
                     py::return_value_policy::reference_internal)
      .def_readwrite("per_pair_posteriors", &DecodePairsReturnStruct::perPairPosteriors,
                     py::return_value_policy::reference_internal)
      .def_readwrite("sum_of_posteriors", &DecodePairsReturnStruct::sumOfPosteriors,
                     py::return_value_policy::reference_internal)
      .def_readwrite("per_pair_posterior_means", &DecodePairsReturnStruct::perPairPosteriorMeans,
                     py::return_value_policy::reference_internal)
      .def_readwrite("min_posterior_means", &DecodePairsReturnStruct::minPosteriorMeans,
                     py::return_value_policy::reference_internal)
      .def_readwrite("argmin_posterior_means", &DecodePairsReturnStruct::argminPosteriorMeans,
                     py::return_value_policy::reference_internal)
      .def_readwrite("per_pair_MAPs", &DecodePairsReturnStruct::perPairMAPs,
                     py::return_value_policy::reference_internal)
      .def_readwrite("min_MAPs", &DecodePairsReturnStruct::minMAPs, py::return_value_policy::reference_internal)
      .def_readwrite("argmin_MAPs", &DecodePairsReturnStruct::argminMAPs, py::return_value_policy::reference_internal);
  py::class_<Individual>(m, "Individual")
      .def(py::init<int>(), "numOfSites"_a = 0)
      .def("setGenotype", &Individual::setGenotype, "hap"_a, "pos"_a, "val"_a)
      .def_readwrite("genotype1", &Individual::genotype1)
      .def_readwrite("genotype2", &Individual::genotype2);
  py::class_<PairObservations>(m, "PairObservations")
      .def_readwrite("obsBits", &PairObservations::obsBits)
      .def_readwrite("homMinorBits", &PairObservations::homMinorBits);
  py::class_<DecodingQuantities>(m, "DecodingQuantities")
      .def(py::init<const char*>())
      .def_readwrite("CSFSSamples", &DecodingQuantities::CSFSSamples)
      .def_readwrite("states", &DecodingQuantities::states)
      .def_readwrite("initialStateProb", &DecodingQuantities::initialStateProb)
      .def_readwrite("expectedTimes", &DecodingQuantities::expectedTimes)
      .def_readwrite("discretization", &DecodingQuantities::discretization)
      .def_readwrite("timeVector", &DecodingQuantities::timeVector)
      .def_readwrite("columnRatios", &DecodingQuantities::columnRatios)
      .def_readwrite("classicEmissionTable", &DecodingQuantities::classicEmissionTable)
      .def_readwrite("compressedEmissionTable", &DecodingQuantities::compressedEmissionTable)
      .def_readwrite("Dvectors", &DecodingQuantities::Dvectors)
      .def_readwrite("Bvectors", &DecodingQuantities::Bvectors)
      .def_readwrite("Uvectors", &DecodingQuantities::Uvectors)
      .def_readwrite("rowRatioVectors", &DecodingQuantities::rowRatioVectors)
      .def_readwrite("homozygousEmissionMap", &DecodingQuantities::homozygousEmissionMap)
      .def_readwrite("CSFSmap", &DecodingQuantities::CSFSmap)
      .def_readwrite("foldedCSFSmap", &DecodingQuantities::foldedCSFSmap)
      .def_readwrite("ascertainedCSFSmap", &DecodingQuantities::ascertainedCSFSmap)
      .def_readwrite("foldedAscertainedCSFSmap", &DecodingQuantities::foldedAscertainedCSFSmap);
  py::class_<DecodingParams>(m, "DecodingParams")
      .def(py::init<std::string, std::string, std::string, int, int, std::string, bool, bool, bool, bool, float, bool,
                    bool, bool, std::string, bool, bool>(),
           "inFileRoot"_a, "decodingQuantFile"_a, "outFileRoot"_a = "", "jobs"_a = 1, "jobInd"_a = 1,
           "decodingModeString"_a = "array", "decodingSequence"_a = false, "usingCSFS"_a = true, "compress"_a = false,
           "useAncestral"_a = false, "skipCSFSdistance"_a = 0.f, "noBatches"_a = false, "doPosteriorSums"_a = false,
           "doPerPairPosteriorMean"_a = false, "expectedCoalTimesFile"_a = "", "withinOnly"_a = false,
           "doMajorMinorPosteriorSums"_a = false)
      .def(py::init<>())
      .def(py::init<std::string, std::string, std::string, bool>(), "in_dir"_a, "decoding_quants"_a, "out_dir"_a,
           "FastSMC"_a = true)
      .def("validateParamsFastSMC", &DecodingParams::validateParamsFastSMC)
      .def_readwrite("inFileRoot", &DecodingParams::inFileRoot)
      .def_readwrite("decodingQuantFile", &DecodingParams::decodingQuantFile)
      .def_readwrite("outFileRoot", &DecodingParams::outFileRoot)
      .def_readwrite("jobs", &DecodingParams::jobs)
      .def_readwrite("jobInd", &DecodingParams::jobInd)
      .def_readwrite("decodingModeString", &DecodingParams::decodingModeString)
      .def_readwrite("decodingMode", &DecodingParams::decodingMode)
      .def_readwrite("decodingSequence", &DecodingParams::decodingSequence)
      .def_readwrite("foldData", &DecodingParams::foldData)
      .def_readwrite("usingCSFS", &DecodingParams::usingCSFS)
      .def_readwrite("compress", &DecodingParams::compress)
      .def_readwrite("useAncestral", &DecodingParams::useAncestral)
      .def_readwrite("skipCSFSdistance", &DecodingParams::skipCSFSdistance)
      .def_readwrite("noBatches", &DecodingParams::noBatches)
      .def_readwrite("batchSize", &DecodingParams::batchSize)
      .def_readwrite("recallThreshold", &DecodingParams::recallThreshold)
      .def_readwrite("min_m", &DecodingParams::min_m)
      .def_readwrite("hashing", &DecodingParams::hashing)
      .def_readwrite("FastSMC", &DecodingParams::FastSMC)
      .def_readwrite("BIN_OUT", &DecodingParams::BIN_OUT)
      .def_readwrite("useKnownSeed", &DecodingParams::useKnownSeed)
      .def_readwrite("outputIbdSegmentLength", &DecodingParams::outputIbdSegmentLength)
      .def_readwrite("hashingWordSize", &DecodingParams::hashingWordSize)
      .def_readwrite("constReadAhead", &DecodingParams::constReadAhead)
      .def_readwrite("haploid", &DecodingParams::haploid)
      .def_readwrite("time", &DecodingParams::time)
      .def_readwrite("noConditionalAgeEstimates", &DecodingParams::noConditionalAgeEstimates)
      .def_readwrite("doPosteriorSums", &DecodingParams::doPosteriorSums)
      .def_readwrite("doPerPairMAP", &DecodingParams::doPerPairMAP)
      .def_readwrite("doPerPairPosteriorMean", &DecodingParams::doPerPairPosteriorMean)
      .def_readwrite("expectedCoalTimesFile", &DecodingParams::expectedCoalTimesFile)
      .def_readwrite("withinOnly", &DecodingParams::withinOnly)
      .def_readwrite("doMajorMinorPosteriorSums", &DecodingParams::doMajorMinorPosteriorSums);

  py::class_<IbdPairDataLine>(m, "IbdPairDataLine")
      .def_readwrite("ind1FamId", &IbdPairDataLine::ind1FamId)
      .def_readwrite("ind1Id", &IbdPairDataLine::ind1Id)
      .def_readwrite("ind1Hap", &IbdPairDataLine::ind1Hap)
      .def_readwrite("ind2FamId", &IbdPairDataLine::ind2FamId)
      .def_readwrite("ind2Id", &IbdPairDataLine::ind2Id)
      .def_readwrite("ind2Hap", &IbdPairDataLine::ind2Hap)
      .def_readwrite("chromosome", &IbdPairDataLine::chromosome)
      .def_readwrite("ibdStart", &IbdPairDataLine::ibdStart)
      .def_readwrite("ibdEnd", &IbdPairDataLine::ibdEnd)
      .def_readwrite("lengthInCentimorgans", &IbdPairDataLine::lengthInCentimorgans)
      .def_readwrite("ibdScore", &IbdPairDataLine::ibdScore)
      .def_readwrite("postEst", &IbdPairDataLine::postEst)
      .def_readwrite("mapEst", &IbdPairDataLine::mapEst)
      .def("toString", &IbdPairDataLine::toString);

  py::class_<BinaryDataReader>(m, "BinaryDataReader")
      .def(py::init<const std::string&>(), "binaryFile"_a)
      .def("getNextLine", &BinaryDataReader::getNextLine)
      .def("moreLinesInFile", &BinaryDataReader::moreLinesInFile);

  py::class_<Data>(m, "Data")
      .def(py::init<const DecodingParams&>(), "params"_a)
      .def_static("countHapLines", &Data::countHapLines)
      .def_readwrite("FamIDList", &Data::FamIDList)
      .def_readwrite("IIDList", &Data::IIDList)
      .def_readwrite("famAndIndNameList", &Data::famAndIndNameList)
      .def_readwrite("individuals", &Data::individuals)
      .def_readwrite("sampleSize", &Data::sampleSize)
      .def_readwrite("haploidSampleSize", &Data::haploidSampleSize)
      .def_readwrite("sites", &Data::sites)
      .def_readwrite("decodingUsesCSFS", &Data::decodingUsesCSFS)
      .def_readwrite("geneticPositions", &Data::geneticPositions)
      .def_readwrite("physicalPositions", &Data::physicalPositions)
      .def_readwrite("siteWasFlippedDuringFolding", &Data::siteWasFlippedDuringFolding)
      .def_readwrite("recRateAtMarker", &Data::recRateAtMarker);
  py::class_<HMM>(m, "HMM")
      .def(py::init<Data&, DecodingParams&, int>(), "data"_a, "params"_a, "scalingSkip"_a = 1)
      .def("decode", py::overload_cast<const PairObservations&>(&HMM::decode))
      .def("decode", py::overload_cast<const PairObservations&, unsigned, unsigned>(&HMM::decode))
      .def("decodeAll", &HMM::decodeAll)
      .def("decodeSummarize", &HMM::decodeSummarize)
      .def("getDecodingReturnValues", &HMM::getDecodingReturnValues)
      .def("decodePair", &HMM::decodePair)
      .def("decodePairs", &HMM::decodePairs)
      .def("getBatchBuffer", &HMM::getBatchBuffer)
      .def("finishDecoding", &HMM::finishDecoding)
      .def("getDecodingQuantities", &HMM::getDecodingQuantities)
      .def("makePairObs", &HMM::makePairObs, "iHap"_a, "ind1"_a, "jHap"_a, "ind2"_a);
  py::class_<ASMC::FastSMC>(m, "FastSMC")
      .def(py::init<DecodingParams>(), "decodingParams"_a)
      .def(py::init<const std::string&, const std::string&>(), "in_dir"_a, "out_dir"_a)
      .def("run", &ASMC::FastSMC::run);
  py::class_<ASMC::ASMC>(m, "ASMC")
      .def(py::init<DecodingParams>(), "decodingParams"_a)
      .def(py::init<const std::string&, const std::string&>(), "in_dir"_a, "out_dir"_a)
      .def("decode_all_in_job", &ASMC::ASMC::decodeAllInJob)
      .def("decode_pairs", py::overload_cast<>(&ASMC::ASMC::decodePairs))
      .def("decode_pairs",
           py::overload_cast<const std::vector<unsigned long>&, const std::vector<unsigned long>&>(
               &ASMC::ASMC::decodePairs),
           "hap_indices_a"_a, "hap_indices_b"_a)
      .def(
          "decode_pairs",
          py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&>(&ASMC::ASMC::decodePairs),
          "hap_ids_a"_a, "hap_ids_b"_a)
      .def("get_copy_of_results", &ASMC::ASMC::getCopyOfResults, py::return_value_policy::copy)
      .def("get_ref_of_results", &ASMC::ASMC::getRefOfResults, py::return_value_policy::reference_internal)
      .def("set_store_per_pair_posterior_mean", &ASMC::ASMC::setStorePerPairPosteriorMean,
           "store_per_pair_posterior_mean"_a = false)
      .def("set_write_per_pair_posterior_mean", &ASMC::ASMC::setWritePerPairPosteriorMean,
           "write_per_pair_posterior_mean"_a = false)
      .def("set_store_per_pair_map", &ASMC::ASMC::setStorePerPairMap, "store_per_pair_map"_a = false)
      .def("set_write_per_pair_map", &ASMC::ASMC::setWritePerPairMap, "write_per_pair_map"_a = false)
      .def("set_store_per_pair_posterior", &ASMC::ASMC::setStorePerPairPosterior, "store_per_pair_posterior"_a = false)
      .def("set_store_sum_of_posterior", &ASMC::ASMC::setStoreSumOfPosterior, "store_sum_of_posterior"_a = false);
}
