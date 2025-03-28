// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef DATA_MODULE_STRING_UTILS_HPP
#define DATA_MODULE_STRING_UTILS_HPP

#include <string>
#include <string_view>
#include <vector>

namespace asmc {

/**
 * Split a string of text into a vector of strings, splitting by a given delimiter.
 *
 * @param text the text to split
 * @param del the delimiter to split by
 * @return a vector of substrings of text, split by del
 */
std::vector<std::string> splitTextByDelimiter(std::string_view text, std::string_view del);

/**
 * Remove whitespace characters '\n', ' ', '\t', '\r' from the end of a string
 *
 * @param s the string to strip whitespace from the back of
 * @return the string with whitespace stripped from the back
 */
std::string stripBack(std::string s);

/**
 * Convert a string to unsigned long, verifying the number was an integer rather than floating point type.
 * A std::runtime_error will be thrown if the string is not representable as an unsigned long.
 *
 * @param s the string to convert to unsigned long
 * @return unsigned long representation of the string
 */
unsigned long ulFromString(const std::string& s);

/**
 * Convert a string to double. A std::runtime_error will be thrown if the string is not representable as a double.
 *
 * @param s the string to convert to double
 * @return double representation of the string
 */
double dblFromString(const std::string& s);

} // namespace asmc

#endif // DATA_MODULE_STRING_UTILS_HPP
