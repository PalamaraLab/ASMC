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


#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <vector>
#include <string>

namespace StringUtils {

/**
 * Convert string to float, taking account of the fact that inputs may be given
 * at a precision too great to be representable as a float.
 *
 * This function converts first to long double, and then explicitly performs a
 * static cast to narrow the output to a float.
 *
 * @param str the string to convert to a float
 * @return the closest float representation of the string
 */
float stof(const std::string &str);

/**
 * Convert string to double, taking account of the fact that inputs may be given
 * at a precision too great to be representable as a double.
 *
 * This function converts first to long double, and then explicitly performs a
 * static cast to narrow the output to a double.
 *
 * @param str the string to convert to a double
 * @return the closest double representation of the string
 */
double stod(const std::string &str);

const std::string RANGE_DELIMS = "{:}";

std::string findDelimiters(const std::string &s, const std::string &c);

std::vector <std::string> tokenizeMultipleDelimiters(const std::string &s, const std::string &c);
std::vector <std::string> expandRangeTemplate(const std::string &str);
std::vector <std::string> expandRangeTemplates(const std::vector <std::string> &rangeTemplates);
}

#endif
