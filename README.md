[![Unit tests: Ubuntu](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-unit.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-unit.yml)
[![Unit tests: macOS](https://github.com/PalamaraLab/ASMC/actions/workflows/macos-unit.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/macos-unit.yml)
[![Python 3.8 3.11](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-python.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-python.yml)
[![Regression test](https://github.com/PalamaraLab/ASMC/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![Ubuntu asan](https://github.com/PalamaraLab/ASMC/workflows/Ubuntu%20asan/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![Ubuntu no sse/avx](https://github.com/PalamaraLab/ASMC/workflows/Ubuntu%20no%20sse/avx/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/ASMC/branch/main/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/ASMC)

# ASMC and FastSMC

This repository contains ASMC and an extension, FastSMC, together with python bindings for both.

The following pages of documentation contains specific information:
- [ASMC](./docs/asmc.md)
- [ASMC python bindings](./docs/asmc_python.md)
- [FastSMC](./docs/fastsmc.md)
- [FastSMC python bindings](./docs/fastsmc_python.md)

## Installation

ASMC and FastSMC are regularly built and tested on Ubuntu and macOS.
They consist of a C++ library, C++ executables, and optional Python bindings.

The C++ libraries and executables require:

- A C++ compiler (C++17 or later)
- CMake (3.15 or later)
- Boost (1.62 or later)
- Eigen (3.3.4 or later)
- {fmt}
- range-v3
- OpenMP
- zlib

We recommend installing dependencies using vcpkg, distributed with this repository as a submodule.
Information below.

Building the optional Python bindings additionally requires:

- Python (3.6 or later) with development files
- PyBind11 (distributed with this repository as a submodule)

## Quickstart guides

- [For users](./docs/quickstart_user.md)
- [For developers](./docs/quickstart_developer.md)

## Decoding Quantities

Decoding quantities files are required in order to run ASMC and FastSMC.
These can be generated directly from a Python module, and instructions can be found [here](https://github.com/PalamaraLab/PrepareDecoding).
Input and output file formats for the tool used to create decoding quantities are described [here](https://github.com/PalamaraLab/PrepareDecoding/blob/master/docs/file_formats.md).

Note: the CEU.demo demographic model and the decoding quantities for CEU+UKBB previously provided in [this repository](https://github.com/PalamaraLab/FastSMC) and [this repository](https://github.com/PalamaraLab/ASMC_legacy) were mistakenly encoded as diploid rather than haploid.
The file [CEU.demo](https://github.com/PalamaraLab/ASMC_data/tree/main/demographies) and CEU+UKBB decoding quantities [here](https://github.com/PalamaraLab/ASMC_data/tree/main/decoding_quantities) have now been fixed.
They were generated using [v2.2.1](https://github.com/PalamaraLab/PrepareDecoding/releases/tag/v2.2.1) of the [PrepareDecoding tool](https://github.com/PalamaraLab/PrepareDecoding), which also provides a simpler interface for computing decoding quantities as well as support for additional demographic models.
Using these new decoding quantities with v1.2 of ASMC will tend to produce more recent estimates for TMRCAs compared to the decoding quantities distributed with v1.0 and v1.1.
This should not have a substantial impact on most downstream analyses.

## For developers: making a release

- Bump the version number in [setup.py](setup.py), and [CMakeLists.txt](CMakeLists.txt)
- Update [RELEASE_NOTES.md](RELEASE_NOTES.md)
- Push changes and check that all [GitHub workflows](https://github.com/PalamaraLab/ASMC/actions) pass
- Tag the commit in Git using syntax `vX.Y`
- Make a release on GitHub, which should trigger a new build that will upload Python wheels to PyPI

## License

ASMC and FastSMC are distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Reference

If you use this software, please cite the appropriate reference(s) below.

The ASMC algorithm and software were developed in
- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

The FastSMC algorithm and software were developed in
- J. Nait Saada, G. Kalantzis, D. Shyr, F. Cooper, M. Robinson, A. Gusev, P. F. Palamara. Identity-by-descent detection across 487,409 British samples reveals fine-scale evolutionary history and trait associations. *Nature Communications*, 2020.
