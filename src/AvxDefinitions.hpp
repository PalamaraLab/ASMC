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


#ifndef AVXDEFINITIONS_HPP
#define AVXDEFINITIONS_HPP

#include <immintrin.h>

#ifdef NO_SSE
#define MODE "NO_SSE"
#define VECX 4
#endif

// SSE vectorization (block size = 4)
#ifdef SSE
#define MODE "SSE"
#define VECX 4
#define FLOAT __m128
#define LOAD _mm_load_ps
#define STORE _mm_store_ps
#define MULT _mm_mul_ps
#define ADD _mm_add_ps
#define RECIPROCAL _mm_rcp_ps
#define LOAD1 _mm_load1_ps
#endif

// AVX vectorization (block size = 8)
#ifdef AVX
#define MODE "AVX"
#include <immintrin.h>
#define VECX 8
#define FLOAT __m256
#define LOAD _mm256_load_ps
#define STORE _mm256_store_ps
#define MULT _mm256_mul_ps
#define ADD _mm256_add_ps
#define RECIPROCAL _mm256_rcp_ps
#define LOAD1 _mm256_broadcast_ss
#endif

// AVX512 vectorization (block size = 16)
#ifdef AVX512
#define MODE "AVX512"
#include <immintrin.h>
#define VECX 16
#define FLOAT __m512
#define LOAD _mm512_load_ps
#define STORE _mm512_store_ps
#define MULT _mm512_mul_ps
#define ADD _mm512_add_ps
#define RECIPROCAL _mm512_rcp14_ps
#define LOAD1 _mm512_set1_ps
#endif


#endif // AVXDEFINITIONS_HPP
