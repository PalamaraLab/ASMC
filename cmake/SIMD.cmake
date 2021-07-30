# fast-math and SIMD instruction settings below has been copied and modified from 
# th GLM library CMakeLists.txt (MIT license)
#
# https://github.com/g-truc/glm/blob/master/CMakeLists.txt

option(ASMC_ENABLE_FAST_MATH "Enable fast math optimizations" OFF)
if(ASMC_ENABLE_FAST_MATH)
  message(STATUS "Build with fast math optimizations")

  if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
    add_compile_options(-ffast-math)

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/fp:fast)
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/fp:precise)
  endif()
endif()

option(ASMC_ENABLE_SIMD_SSE3 "Enable SSE 1, 2 & 3 optimizations" OFF)
option(ASMC_ENABLE_SIMD_AVX "Enable AVX optimizations" ON)
option(ASMC_ENABLE_SIMD_AVX512 "Enable AVX512 optimizations" OFF)
option(ASMC_FORCE_PURE "Force 'pure' instructions" OFF)

if(ASMC_FORCE_PURE)
  add_definitions(-DNO_SSE)

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-mfpmath=387)
  endif()
  message(STATUS "No SIMD instruction set")

elseif(ASMC_ENABLE_SIMD_AVX)
  add_definitions(-DAVX)
  add_compile_definitions(EIGEN_MAX_ALIGN_BYTES=64)

  if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    add_compile_options(-mavx)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(/QxAVX)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/arch:AVX)
  endif()
  message(STATUS "AVX instruction set")

elseif(ASMC_ENABLE_SIMD_AVX512)
  add_definitions(-DAVX)

  if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    add_compile_options(-mavx512f)
    add_compile_options(-mavx512cd)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-xCOMMON-AVX512)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/arch:AVX512)
  endif()
  message(STATUS "AVX-512 instruction set")

elseif(ASMC_ENABLE_SIMD_SSE)
  add_definitions(-DSSE)

  if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    add_compile_options(-msse3)
    add_compile_options(-msse2)
    add_compile_options(-msse)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(/QxSSE3)
    add_compile_options(/QxSSE2)
    add_compile_options(/QxSSE)
  elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC"))
    add_compile_options(/arch:SSE2) # VC doesn't support SSE3
    add_compile_options(/arch:SSE) 
  endif()
  message(STATUS "SSE2 & SSE3 instruction set")

elseif(ASMC_ENABLE_SIMD_SSE2)
  add_definitions(-DGLM_FORCE_INTRINSICS)

  if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    add_compile_options(-msse2)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(/QxSSE2)
  elseif((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND NOT CMAKE_CL_64)
    add_compile_options(/arch:SSE2)
  endif()
  message(STATUS "SSE2 instruction set")
endif()

