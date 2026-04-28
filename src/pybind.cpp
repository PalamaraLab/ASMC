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
#include "Cluster.hpp"
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

#include <limits>
#include <vector>

namespace {
// Helpers shared by from_arrays/from_counts pybind factories

// Parameters for in-memory input. The validating constructor checks the decoding mode (in any case) and sets the
// folding and CSFS options from use_ancestral and skip_csfs_distance, as for file input. It reads no input files.
// usingCSFS starts false and is turned on unless skip_csfs_distance is infinite.
DecodingParams makeDecodingParams(const std::string& dq_file, bool use_ancestral, float skip_csfs_distance,
                                  const std::string& decoding_mode)
{
  return DecodingParams(/*_inFileRoot=*/"", dq_file, /*_outFileRoot=*/"", /*_jobs=*/1, /*_jobInd=*/1, decoding_mode,
                        /*_decodingSequence=*/false, /*_usingCSFS=*/false, /*_compress=*/false, use_ancestral,
                        skip_csfs_distance);
}

std::pair<std::vector<int>, std::vector<float>> copyPositions(
    pybind11::array_t<int, pybind11::array::c_style | pybind11::array::forcecast> physical_positions,
    pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast> genetic_positions,
    int num_sites) {
  auto phys = physical_positions.unchecked<1>();
  auto gen = genetic_positions.unchecked<1>();
  std::vector<int> physPos(num_sites);
  std::vector<float> genPos(num_sites);
  for (int i = 0; i < num_sites; i++) { physPos[i] = phys(i); genPos[i] = gen(i); }
  return {std::move(physPos), std::move(genPos)};
}
} // namespace

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
                    bool, bool, std::string, bool, bool, bool, std::string>(),
           "in_file_root"_a, "dq_file"_a, "out_file_root"_a = "", "jobs"_a = 1, "job_ind"_a = 1,
           "decoding_mode_string"_a = "array", "decoding_sequence"_a = false, "using_CSFS"_a = true,
           "compress"_a = false, "use_ancestral"_a = false,
           "skip_CSFS_distance"_a = std::numeric_limits<float>::quiet_NaN(), "no_batches"_a = false,
           "do_posterior_sums"_a = false, "do_per_pair_posterior_mean"_a = false, "expected_coal_times_file"_a = "",
           "within_only"_a = false, "do_major_minor_posterior_sums"_a = false, "do_per_pair_MAP"_a = false,
           "map_file"_a = "")
      .def(py::init<>())
      .def(py::init<std::string, std::string, std::string, bool>(), "in_file_root"_a, "dq_file"_a,
           "out_file_root"_a, "FastSMC"_a = true)
      .def("validateParamsFastSMC", &DecodingParams::validateParamsFastSMC)
      .def("print_decoding_params", &DecodingParams::printDecodingParams)
      .def_readwrite("inFileRoot", &DecodingParams::inFileRoot)
      .def_readwrite("decodingQuantFile", &DecodingParams::decodingQuantFile)
      .def_readwrite("mapFile", &DecodingParams::mapFile)
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
      .def(py::init<const std::string&, const std::string&, const std::string&>(), "in_dir"_a, "dq_file"_a, "out_dir"_a)
      .def("run", &ASMC::FastSMC::run);
  py::class_<ASMC::ASMC>(m, "ASMC")
      .def(py::init<DecodingParams>(), "decodingParams"_a)
      .def(py::init<const std::string&, const std::string&, const std::string&, const std::string&>(), "in_dir"_a,
           "dq_file"_a, "out_dir"_a = "", "decoding_mode"_a = "array")
      .def("get_decoding_params", &ASMC::ASMC::getDecodingParams)
      .def("get_diploid_sample_size", &ASMC::ASMC::getDiploidSampleSize)
      .def("get_haploid_sample_size", &ASMC::ASMC::getHaploidSampleSize)
      .def("get_num_sites", &ASMC::ASMC::getNumSites)
      .def("get_physical_positions", &ASMC::ASMC::getPhysicalPositions)
      .def("get_genetic_positions", &ASMC::ASMC::getGeneticPositions)
      .def("decode_all_in_job", &ASMC::ASMC::decodeAllInJob)
      .def("decode_pairs", py::overload_cast<unsigned, unsigned, float>(&ASMC::ASMC::decodePairs), "from"_a = 0, "to"_a = 0, "cm_burn_in"_a = 0.5)
      .def("decode_pairs",
           py::overload_cast<const std::vector<unsigned long>&, const std::vector<unsigned long>&, unsigned, unsigned, float>(
               &ASMC::ASMC::decodePairs),
           "hap_indices_a"_a, "hap_indices_b"_a, "from"_a = 0, "to"_a = 0, "cm_burn_in"_a = 0.5)
      .def(
          "decode_pairs",
          py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&, unsigned, unsigned, float>(&ASMC::ASMC::decodePairs),
          "hap_ids_a"_a, "hap_ids_b"_a, "from"_a = 0, "to"_a = 0, "cm_burn_in"_a = 0.5)
      .def("decode_pairs_parallel", &ASMC::ASMC::decodePairsParallel,
           "hap_indices_a"_a, "hap_indices_b"_a, "from"_a = 0, "to"_a = 0, "cm_burn_in"_a = 0.5, "n_threads"_a = 4)
      .def("get_copy_of_results", &ASMC::ASMC::getCopyOfResults, py::return_value_policy::copy)
      .def("get_ref_of_results", &ASMC::ASMC::getRefOfResults, py::return_value_policy::reference_internal)
      .def("get_expected_times", &ASMC::ASMC::getExpectedTimes, py::return_value_policy::reference_internal)
      .def("set_store_per_pair_posterior_mean", &ASMC::ASMC::setStorePerPairPosteriorMean,
           "store_per_pair_posterior_mean"_a = true)
      .def("set_write_per_pair_posterior_mean", &ASMC::ASMC::setWritePerPairPosteriorMean,
           "write_per_pair_posterior_mean"_a = true)
      .def("set_store_per_pair_map", &ASMC::ASMC::setStorePerPairMap, "store_per_pair_map"_a = true)
      .def("set_write_per_pair_map", &ASMC::ASMC::setWritePerPairMap, "write_per_pair_map"_a = true)
      .def("set_store_per_pair_posterior", &ASMC::ASMC::setStorePerPairPosterior, "store_per_pair_posterior"_a = true)
      .def("set_store_sum_of_posterior", &ASMC::ASMC::setStoreSumOfPosterior, "store_sum_of_posterior"_a = true)
      // Construct from in-memory genotype array, bypassing all file I/O.
      .def_static("from_arrays", [](
            py::array_t<int8_t, py::array::c_style | py::array::forcecast> genotypes,
            py::array_t<int, py::array::c_style | py::array::forcecast> physical_positions,
            py::array_t<float, py::array::c_style | py::array::forcecast> genetic_positions,
            const std::string& dq_file, bool use_ancestral,
            float skip_csfs_distance, const std::string& decoding_mode
          ) -> std::unique_ptr<ASMC::ASMC> {
            auto geno = genotypes.unchecked<2>();
            int num_sites = static_cast<int>(geno.shape(0));
            int num_haps = static_cast<int>(geno.shape(1));
            if (num_haps % 2 != 0)
              throw std::invalid_argument("genotypes num_haps (axis 1) must be even (diploid)");
            if (physical_positions.size() != num_sites || genetic_positions.size() != num_sites)
              throw std::invalid_argument("position arrays must have length == num_sites");
            auto num_diploid = static_cast<unsigned long>(num_haps / 2);
            auto [physPos, genPos] = copyPositions(physical_positions, genetic_positions, num_sites);
            auto params = makeDecodingParams(dq_file, use_ancestral, skip_csfs_distance, decoding_mode);
            Data data(num_sites, num_diploid, std::move(physPos), std::move(genPos), genotypes.data(), params.usingCSFS,
                      params.foldData);
            return std::make_unique<ASMC::ASMC>(std::move(params), std::move(data));
          },
          "genotypes"_a, "physical_positions"_a, "genetic_positions"_a,
          "dq_file"_a, "use_ancestral"_a = false, "skip_csfs_distance"_a = 0.f,
          "decoding_mode"_a = "array",
          "Construct ASMC from in-memory genotype array (num_sites x num_haps), bypassing file I/O.")
      // Construct from precomputed allele counts only (no genotypes needed).
      // Genotypes are loaded per-batch via load_genotypes() for streaming decoding.
      .def_static("from_counts", [](
            py::array_t<int, py::array::c_style | py::array::forcecast> derived_allele_counts,
            py::array_t<int, py::array::c_style | py::array::forcecast> total_sample_counts,
            py::array_t<int, py::array::c_style | py::array::forcecast> physical_positions,
            py::array_t<float, py::array::c_style | py::array::forcecast> genetic_positions,
            const std::string& dq_file, bool use_ancestral,
            float skip_csfs_distance, const std::string& decoding_mode
          ) -> std::unique_ptr<ASMC::ASMC> {
            auto dac = derived_allele_counts.unchecked<1>();
            auto tsc = total_sample_counts.unchecked<1>();
            int num_sites = static_cast<int>(dac.shape(0));
            if (tsc.shape(0) != num_sites)
              throw std::invalid_argument("total_sample_counts must have same length as derived_allele_counts");
            if (physical_positions.size() != num_sites || genetic_positions.size() != num_sites)
              throw std::invalid_argument("position arrays must have length == num_sites");
            auto [physPos, genPos] = copyPositions(physical_positions, genetic_positions, num_sites);
            std::vector<int> dacVec(num_sites), tscVec(num_sites);
            for (int i = 0; i < num_sites; i++) {
              dacVec[i] = dac(i);
              tscVec[i] = tsc(i);
            }
            auto params = makeDecodingParams(dq_file, use_ancestral, skip_csfs_distance, decoding_mode);
            Data data(std::move(dacVec), std::move(tscVec), std::move(physPos), std::move(genPos), params.usingCSFS,
                      params.foldData);
            return std::make_unique<ASMC::ASMC>(std::move(params), std::move(data));
          },
          "derived_allele_counts"_a, "total_sample_counts"_a,
          "physical_positions"_a, "genetic_positions"_a,
          "dq_file"_a, "use_ancestral"_a = false, "skip_csfs_distance"_a = 0.f,
          "decoding_mode"_a = "array",
          "Construct ASMC from allele counts only. Use load_genotypes()/clear_genotypes() to stream haplotype subsets.")
      // Load/clear genotypes for streaming pair decoding (used with from_counts).
      .def("load_genotypes", [](ASMC::ASMC& self,
            py::array_t<int8_t, py::array::c_style | py::array::forcecast> genotypes) {
            if (genotypes.ndim() != 2)
              throw std::invalid_argument("genotypes must be 2D (num_sites x num_haps)");
            if (static_cast<int>(genotypes.shape(0)) != self.getNumSites())
              throw std::invalid_argument("genotypes num_sites (axis 0) must match ASMC num_sites");
            if (genotypes.shape(1) % 2 != 0)
              throw std::invalid_argument("genotypes num_haps (axis 1) must be even (diploid)");
            self.loadSubsetGenotypes(genotypes.data(), static_cast<int>(genotypes.shape(1)));
          }, "genotypes"_a, "Load genotypes (num_sites x num_subset_haps) for a batch of pair decoding.")
      .def("clear_genotypes", &ASMC::ASMC::clearGenotypes, "Clear loaded genotypes to free memory between batches.")
      // ASMC-clust: per-site UPGMA over posterior-mean TMRCAs. Returns a tuple of
      // numpy arrays so callers can ingest it into tskit / arg-needle-lib without
      // copying. See src/Cluster.hpp for layout details.
      .def("cluster_chunk", [](
            ASMC::ASMC& self,
            unsigned from, unsigned to,
            unsigned site_block, int n_threads,
            float cm_burn_in, int arg_offset, double arg_length,
            unsigned anchor_stride_k
          ) -> py::tuple {
            ASMC::ClusterChunkResult r =
                ASMC::clusterChunk(self, from, to, site_block, n_threads,
                                   cm_burn_in, arg_offset, arg_length,
                                   anchor_stride_k);
            // Move the vectors into NumPy arrays without copying, by handing
            // ownership to capsule-managed allocations.
            auto take_vec = [](auto&& vec) {
              using T = typename std::decay_t<decltype(vec)>::value_type;
              auto* heap = new std::vector<T>(std::move(vec));
              py::capsule owner(heap, [](void* p) {
                delete reinterpret_cast<std::vector<T>*>(p);
              });
              return py::array_t<T>(
                  static_cast<py::ssize_t>(heap->size()),
                  heap->data(), owner);
            };
            py::dict result;
            result["num_haps"]        = r.num_haps;
            result["num_sites"]       = r.num_sites;
            result["site_boundaries"] = take_vec(std::move(r.site_boundaries));
            result["node_time"]       = take_vec(std::move(r.node_time));
            result["edge_left"]       = take_vec(std::move(r.edge_left));
            result["edge_right"]      = take_vec(std::move(r.edge_right));
            result["edge_parent"]     = take_vec(std::move(r.edge_parent));
            result["edge_child"]      = take_vec(std::move(r.edge_child));
            return py::tuple(py::make_tuple(result));
          },
          "from_"_a = 0u, "to"_a = 0u,
          "site_block"_a = 0u, "n_threads"_a = 0,
          "cm_burn_in"_a = 0.5f, "arg_offset"_a = 0,
          "arg_length"_a = 0.0,
          "anchor_stride_k"_a = 1u,
          "ASMC-clust: per-site UPGMA over posterior-mean TMRCAs in [from, to). "
          "anchor_stride_k>1 runs UPGMA only at every K-th site (anchors) and "
          "copies topology+heights to non-anchor sites; K=1 is full per-site UPGMA.");
}
