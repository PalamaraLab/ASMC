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


#ifndef TYPES_HPP
#define TYPES_HPP

#include <cstdint>
#include <cinttypes>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef uint64_t uint64;
typedef int64_t int64;
typedef uint64_t hash_size;

struct uint64_masks {
	uint64 is0, is2, is9;
};

#endif
