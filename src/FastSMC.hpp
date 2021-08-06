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

#ifndef ASMC_FASTSMC_HPP
#define ASMC_FASTSMC_HPP

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "HMM.hpp"

namespace ASMC
{

class FastSMC
{

private:

  DecodingParams mParams;
  Data mData;
  HMM mHmm;

public:

  /**
   * FastSMC constructor with full control over parameters, by manually specifying a DecodingParams object.
   *
   * @param params the decoding parameters
   */
  explicit FastSMC(DecodingParams params);

  /**
   * FastSMC constructor that will set sensible defaults. If you wish to fine-tune parameters, use the constructor that
   * takes a DecodingParams object, which you can configure manually.
   *
   * @param inFileRoot the input file root
   * @param dqFile the decoding quantities file
   * @param outFileRoot the output file root
   */
  FastSMC(const std::string& inFileRoot, const std::string& dqFile, const std::string& outFileRoot);

  void run();

};

} // namespace ASMC

#endif // ASMC_FASTSMC_HPP
