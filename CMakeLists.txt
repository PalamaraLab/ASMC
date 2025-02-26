# This file is part of https://github.com/PalamaraLab/ASMC_new_dev which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

cmake_minimum_required(VERSION 3.15)
message(STATUS "CMake version: ${CMAKE_VERSION}")

include(cmake/CheckDataModule.cmake)

project(asmc LANGUAGES CXX VERSION 1.3.1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    message(STATUS "  You can specify with the flag -DCMAKE_BUILD_TYPE=<Debug|Release|MinSizeRel|RelWithDebInfo>")
    message(STATUS "  'Release' will build optimised binaries, but 'Debug' may be better while developing FastSMC.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()

message(STATUS "CMake build type is set to ${CMAKE_BUILD_TYPE}")

# ensure C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#export compiler flags for code completion engines
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

include(FetchContent)

option(ASMC_PYTHON_BINDINGS "Whether to build the python bindings" OFF)
if (ASMC_PYTHON_BINDINGS)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11
            GIT_TAG 7c33cdc2d39c7b99a122579f53bc94c8eb3332ff # Version 2.13.5
    )
    FetchContent_MakeAvailable(pybind11)
endif ()

# Set global compiler warnings for debug builds
if(CMAKE_BUILD_TYPE MATCHES "Deb")
    if (MSVC)
        add_compile_options(/W3 /WX)
    else ()
        add_compile_options(-Wall -pedantic -fPIC)
    endif ()
endif()

# Set SIMD compiler options
include(cmake/SIMD.cmake)

set(ASMC_exe_dir ${CMAKE_CURRENT_SOURCE_DIR}/exe)
set(ASMC_source_dir ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(ASMC_test_dir ${CMAKE_CURRENT_SOURCE_DIR}/test)
set(ASMC_data_module_dir ${CMAKE_CURRENT_SOURCE_DIR}/DataModule)
set(ASMC_data_dir ${CMAKE_CURRENT_SOURCE_DIR}/ASMC_data)
add_definitions(-DASMC_TEST_DIR=\"${ASMC_test_dir}\")
add_definitions(-DASMC_DATA_DIR=\"${ASMC_data_dir}\")

option(ASMC_TESTING "Whether to generate testing targets" ON)
if (NOT EXISTS ${ASMC_data_dir}/README.md)
    message(STATUS "The data directory ${ASMC_data_dir} does not exist, so testing will not work.")
    message(STATUS "For testing functionality, please check out all submodules.")
    set(ASMC_TESTING OFF)
endif ()

set(ASMC_data_module_source
        ${ASMC_data_module_dir}/src/GeneticMap.cpp
        ${ASMC_data_module_dir}/src/PlinkMap.cpp
        ${ASMC_data_module_dir}/src/utils/FileUtils.cpp
        ${ASMC_data_module_dir}/src/utils/StringUtils.cpp
        )

set(ASMC_data_module_headers
        ${ASMC_data_module_dir}/src/GeneticMap.hpp
        ${ASMC_data_module_dir}/src/PlinkMap.hpp
        ${ASMC_data_module_dir}/src/utils/FileUtils.hpp
        ${ASMC_data_module_dir}/src/utils/StringUtils.hpp
        ${ASMC_data_module_dir}/src/utils/VectorUtils.hpp
        )

# The ASMC Data Module library
add_library(ASMC_data_module ${ASMC_data_module_source} ${ASMC_data_module_headers})
target_include_directories(ASMC_data_module PUBLIC ${ASMC_data_module_dir}/src)

set(ASMC_source
        ${ASMC_source_dir}/HMM.cpp
        ${ASMC_source_dir}/HmmUtils.cpp
        ${ASMC_source_dir}/BinaryDataReader.cpp
        ${ASMC_source_dir}/Data.cpp
        ${ASMC_source_dir}/DecodingParams.cpp
        ${ASMC_source_dir}/DecodingQuantities.cpp
        ${ASMC_source_dir}/FastSMC.cpp
        ${ASMC_source_dir}/FileUtils.cpp
        ${ASMC_source_dir}/Individual.cpp
        ${ASMC_source_dir}/MemoryUtils.cpp
        ${ASMC_source_dir}/Simd.cpp
        ${ASMC_source_dir}/StringUtils.cpp
        ${ASMC_source_dir}/Timer.cpp
        ${ASMC_source_dir}/ASMC.cpp
        ${ASMC_source_dir}/hashing/Utils.cpp
        )

set(ASMC_headers
        ${ASMC_source_dir}/AvxDefinitions.hpp
        ${ASMC_source_dir}/BinaryDataReader.hpp
        ${ASMC_source_dir}/HMM.hpp
        ${ASMC_source_dir}/HmmUtils.hpp
        ${ASMC_source_dir}/Data.hpp
        ${ASMC_source_dir}/DecodePairsReturnStruct.hpp
        ${ASMC_source_dir}/DecodingParams.hpp
        ${ASMC_source_dir}/DecodingQuantities.hpp
        ${ASMC_source_dir}/FastSMC.hpp
        ${ASMC_source_dir}/FileUtils.hpp
        ${ASMC_source_dir}/IbdPairDataLine.hpp
        ${ASMC_source_dir}/Individual.hpp
        ${ASMC_source_dir}/MemoryUtils.hpp
        ${ASMC_source_dir}/Simd.hpp
        ${ASMC_source_dir}/StringUtils.hpp
        ${ASMC_source_dir}/Timer.hpp
        ${ASMC_source_dir}/ASMC.hpp
        ${ASMC_source_dir}/hashing/ExtendHash.hpp
        ${ASMC_source_dir}/hashing/Individuals.hpp
        ${ASMC_source_dir}/hashing/Match.hpp
        ${ASMC_source_dir}/hashing/SeedHash.hpp
        ${ASMC_source_dir}/hashing/Utils.hpp
        )

set(ASMC_main
        ${ASMC_exe_dir}/main.cpp
        )

set(FastSMC_main
        ${ASMC_exe_dir}/main_fastsmc.cpp
        )

set(convertBinary_main
        ${ASMC_exe_dir}/main_convertBinary.cpp
        )

set(ASMC_unit_tests
        ${ASMC_test_dir}/unit_tests.cpp
        ${ASMC_test_dir}/test_ASMC.cpp
        ${ASMC_test_dir}/test_binary_data_reader.cpp
        ${ASMC_test_dir}/test_hashing.cpp
        ${ASMC_test_dir}/test_HMM.cpp
        ${ASMC_test_dir}/test_hmm_utils.cpp
        ${ASMC_test_dir}/test_decoding_params.cpp
        ${ASMC_test_dir}/test_decoding_quantities.cpp
        ${ASMC_test_dir}/test_SIMD.cpp
        )

set(ASMC_regression
        ${ASMC_test_dir}/test_regression.cpp
        ${ASMC_test_dir}/test_fastsmc_regression.cpp
        )

if (ASMC_PYTHON_BINDINGS)
    pybind11_add_module(asmc_python_bindings ${ASMC_source_dir}/pybind.cpp)
    target_link_libraries(asmc_python_bindings PRIVATE ASMC)
endif ()

# The main ASMC library
add_library(ASMC ${ASMC_source} ${ASMC_headers})
target_include_directories(ASMC PUBLIC ${ASMC_source_dir})
target_link_libraries(ASMC PRIVATE ASMC_data_module)

# The main ASMC executable
add_executable(ASMC_exe ${ASMC_main})
target_link_libraries(ASMC_exe PRIVATE ASMC)

# The main FastSMC executable
add_executable(FastSMC_exe ${FastSMC_main})
target_link_libraries(FastSMC_exe PRIVATE ASMC)

# The main convertBinary executable
add_executable(convertBinary_exe ${convertBinary_main})
target_link_libraries(convertBinary_exe PRIVATE ASMC)

# boost as a required dependency for ASMC
find_package(Boost CONFIG REQUIRED COMPONENTS iostreams program_options)
target_link_libraries(ASMC PRIVATE Boost::iostreams Boost::program_options)
target_link_libraries(ASMC PUBLIC Boost::boost)

# zlib is required (at least on Linux) for boost iostreams
find_package(ZLIB REQUIRED)
target_link_libraries(ASMC PRIVATE ZLIB::ZLIB)

include(FetchContent)
set(FETCHCONTENT_QUIET ON)

# Fetch Eigen3
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "Disable range-v3 tests")
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        68f4e58cfacc686583d16cff90361f0b43bc2c1b # Latest commit on 3.4 branch
)
FetchContent_MakeAvailable(eigen)
target_link_libraries(ASMC PUBLIC Eigen3::Eigen)

# Fetch fmt
set(FMT_TEST OFF CACHE BOOL "Disable fmt tests")
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281 # v10.2.1
)
FetchContent_MakeAvailable(fmt)
target_link_libraries(ASMC PRIVATE fmt::fmt)

# Fetch range-v3
set(RANGE_V3_TESTS OFF CACHE BOOL "Disable range-v3 tests")
FetchContent_Declare(
        range-v3
        GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
        GIT_TAG        a81477931a8aa2ad025c6bda0609f38e09e4d7ec  # v0.12.0
)
FetchContent_MakeAvailable(range-v3)

# Fetch Highway
set(HWY_ENABLE_TESTS OFF CACHE BOOL "Disable Highway tests")
FetchContent_Declare(
        highway
        GIT_REPOSITORY https://github.com/google/highway.git
        GIT_TAG        457c891775a7397bdb0376bb1031e6e027af1c48 # v1.2.0
)
FetchContent_MakeAvailable(highway)
target_link_libraries(ASMC PRIVATE hwy)
target_compile_definitions(ASMC PRIVATE HWY_DYNAMIC_DISPATCH=1)

target_link_libraries(ASMC_data_module PRIVATE Eigen3::Eigen fmt::fmt range-v3 ZLIB::ZLIB)


if (ASMC_TESTING)
    enable_testing()

    add_executable(ASMC_unit_tests ${ASMC_unit_tests})
    target_include_directories(ASMC_unit_tests PRIVATE ${ASMC_test_dir})
    target_link_libraries(ASMC_unit_tests PRIVATE ASMC fmt::fmt)
    add_test(Asmc_unit_tests ASMC_unit_tests)

    add_executable(ASMC_regression ${ASMC_regression})
    target_include_directories(ASMC_regression PRIVATE ${ASMC_test_dir})
    target_link_libraries(ASMC_regression PRIVATE ASMC ASMC_data_module)
    add_test(regression ASMC_regression)
endif ()


# Optionally run regression tests to ensure core functionality remains unchanged
if (NOT ASMC_NO_PYTHON)
    option(ASMC_REGRESSION_TESTING "Whether to generate regression testing target" OFF)
    if (ASMC_REGRESSION_TESTING)

        # Python is required for running regression tests
        find_package(Python3 3.9 COMPONENTS Interpreter REQUIRED)

        add_custom_target(
                Regression
                COMMAND ${Python3_EXECUTABLE} regression_test.py $<TARGET_FILE:ASMC_exe>
                DEPENDS ASMC_exe
                COMMENT "Running ASMC regression test"
                WORKING_DIRECTORY ${ASMC_test_dir}
        )

    endif ()
endif ()

# Need position independent code for Python bindings
#if (ASMC_PYTHON_BINDINGS)
#    set_target_properties(ASMC PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#    set_target_properties(ASMC_data_module PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
#endif ()

# Clang tidy as optional static analyzer
option(Template_USE_CLANG_TIDY "Use clang tidy for static analysis" OFF)
if (Template_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-6.0" "clang-tidy-5.0" "clang-tidy-4.0"
            DOC "Path to clang-tidy executable")
    if (NOT CLANG_TIDY_EXE)
        message(STATUS "clang-tidy not found.")
    else ()
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set_target_properties(ASMC PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif ()
endif ()

# LLVM AddressSanitizer to detecting memory errors
# Note that there are many other sanitizers in LLVM to help detect errors, see
# http://travistoptips.blogspot.co.uk/2015/11/sanitize-all-things.html
option(ASMC_MEMCHECK "Use LLVM AddressSanitizer to detecting memory errors" OFF)
if (ASMC_MEMCHECK)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Configuring with LLVM AddressSanitizer")
        set(ASMC_MEMCHECK_FLAGS -fno-optimize-sibling-calls
                -fsanitize=address
                -fsanitize-address-use-after-scope
                )
        target_compile_options(ASMC PUBLIC -O1 -g -fno-omit-frame-pointer ${ASMC_MEMCHECK_FLAGS})
        target_link_libraries(ASMC PUBLIC -g ${ASMC_MEMCHECK_FLAGS})
    else ()
        message(FATAL_ERROR "clang compiler required with ASMC_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()

# Setup coverage testing for GCC or Clang
option(ASMC_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (ASMC_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message(STATUS "Configuring with coverage")
        target_compile_options(ASMC PUBLIC --coverage -O0)
        target_link_libraries(ASMC PUBLIC --coverage)
    else ()
        message(FATAL_ERROR "GCC or Clang required with ASMC_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()
