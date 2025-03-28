// This file is part of https://github.com/PalamaraLab/DataModule which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "StringUtils.hpp"

#include <exception>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/core.h>

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/split.hpp>

namespace asmc {

std::vector<std::string> splitTextByDelimiter(std::string_view text, std::string_view del) {
  return text | ranges::views::split(del) | ranges::to<std::vector<std::string>>();
}

std::string stripBack(std::string s) {
  while (!s.empty() && (s.back() == '\n' || s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) {
    s.pop_back();
  }
  return s;
}

unsigned long ulFromString(const std::string& s) {

  bool valid = true;
  long double ld{};
  unsigned long ul{};

  try {
    ld = std::stold(s);
  } catch (const std::exception&) {
    valid = false;
  }

  try {
    ul = std::stoul(s);
  } catch (const std::exception&) {
    valid = false;
  }

  if (!valid || static_cast<long double>(ul) != ld) {
    throw std::runtime_error(fmt::format("String {} not representable as an unsigned integer\n", s));
  }

  return ul;
}

double dblFromString(const std::string& s) {
  try {
    return static_cast<double>(std::stold(s));
  } catch (const std::exception&) {
    throw std::runtime_error(fmt::format("String {} not representable as a double\n", s));
  }
}

} // namespace asmc
