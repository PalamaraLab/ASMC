# This file is part of ASMC, developed by Pier Francesco Palamara.

# ASMC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ASMC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ASMC.  If not, see <https://www.gnu.org/licenses/>.


from .asmc_python_bindings import BinaryDataReader
from .asmc_python_bindings import DecodingModeOverall
from .asmc_python_bindings import DecodingMode
from .asmc_python_bindings import DecodingReturnValues
from .asmc_python_bindings import DecodePairsReturnStruct
from .asmc_python_bindings import IbdPairDataLine
from .asmc_python_bindings import Individual
from .asmc_python_bindings import PairObservations
from .asmc_python_bindings import DecodingQuantities
from .asmc_python_bindings import DecodingParams
from .asmc_python_bindings import Data
from .asmc_python_bindings import HMM
from .asmc_python_bindings import FastSMC
from .asmc_python_bindings import ASMC


#
# ASMCReturnValues = collections.namedtuple(
#     "ASMCReturnValues",
#     "sumOverPairs sumOverPairs00 sumOverPairs01 sumOverPairs11")


# def to_array(x):
#     a = list(x)
#     if a:
#         return np.array(a)
#     else:
#         return None
#
#
# def flip_rows(a1, a2, flips):
#     # Swap rows according to boolean flips vector
#     if a1 is None or a2 is None:
#         return None, None
#     a1[flips], a2[flips] = a2[flips], a1[flips]
#     return a1, a2


# def run(in_file_root, decoding_quant_file, out_file_root="",
#         mode=DecodingModeOverall.array, jobs=0,
#         job_index=0, skip_csfs_distance=0,
#         compress=False, use_ancestral=False,
#         posterior_sums=False, major_minor_posterior_sums=False):
#     ret = asmc(in_file_root=in_file_root,
#                decoding_quant_file=decoding_quant_file,
#                mode=mode, jobs=jobs, job_index=job_index,
#                skip_csfs_distance=skip_csfs_distance,
#                compress=compress, use_ancestral=use_ancestral,
#                posterior_sums=posterior_sums,
#                major_minor_posterior_sums=major_minor_posterior_sums)
#     sumOverPairs00, sumOverPairs11 = flip_rows(
#         to_array(ret.sumOverPairs00), to_array(ret.sumOverPairs11),
#         ret.siteWasFlippedDuringFolding)
#     return ASMCReturnValues(
#         sumOverPairs=to_array(ret.sumOverPairs),
#         sumOverPairs00=sumOverPairs00,
#         sumOverPairs01=to_array(ret.sumOverPairs01),
#         sumOverPairs11=sumOverPairs11)
