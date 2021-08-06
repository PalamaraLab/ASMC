[![Ubuntu unit tests](https://github.com/PalamaraLab/FastSMC/workflows/Ubuntu%20unit/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![macOS unit tests](https://github.com/PalamaraLab/FastSMC/workflows/macOS%20unit/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![Python tests](https://github.com/PalamaraLab/FastSMC/workflows/Python%203.5%203.8/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![Regression test](https://github.com/PalamaraLab/FastSMC/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![Ubuntu asan](https://github.com/PalamaraLab/FastSMC/workflows/Ubuntu%20asan/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![Ubuntu no sse/avx](https://github.com/PalamaraLab/FastSMC/workflows/Ubuntu%20no%20sse/avx/badge.svg)](https://github.com/PalamaraLab/FastSMC/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/FastSMC/branch/master/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/FastSMC)

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
- CMake (3.12 or later)
- Boost (1.62 or later)
- Eigen (3.3.4 or later)
- {fmt}
- zlib

We recommend installing these dependencies using vcpkg, distributed with this repository as a submodule.
Information below.

To build the Python bindings you additionally require:

- Python (3.6 or later)
- PyBind11 (distributed with this repository as a submodule)

## Getting the code

**For code users**
```bash
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC
cd ASMC
```

**For code developers**
```bash
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC_dev
cd ASMC
```

## Install dependencies

### Using vcpkg (recommended)

Dependencies will be automatically built and installed during the CMake configuration step.
No additional steps are required to install dependencies.

### Using your system package manager

**Ubuntu**
```bash
sudo apt install g++ cmake libboost-all-dev libeigen3-dev
```

**macOS (using homebrew and assuming Xcode is installed)**
```bash
brew install cmake boost eigen
```

## Build the project

These instructions will build the library, executables and Python bindings.

```bash
mkdir build && cd build
cmake ..
cmake --build . --parallel 4
```

Note: you can locate the build directory anywhere you like: just run `cmake /path/to/ASMC`.

## Install the Python module

For code development, install the python module locally with `pip` from the top-level ASMC directory:

```bash
pip install .
```

If you are a code user, instead install directly from `PyPI` following instructions in [ASMC python bindings](./docs/asmc_python.md) or [FastSMC python bindings](./docs/fastsmc_python.md).

## Decoding Quantities

Decoding quantities files are required in order to run ASMC and FastSMC.
These can be generated directly from a Python module, and instructions can be found [here](https://github.com/PalamaraLab/PrepareDecoding).

## License

ASMC and FastSMC are distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Reference

If you use this software, please cite the appropriate reference(s) below.

The ASMC algorithm and software were developed in
- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

The FastSMC algorithm and software were developed in
- J. Nait Saada, G. Kalantzis, D. Shyr, F. Cooper, M. Robinson, A. Gusev, P. F. Palamara. Identity-by-descent detection across 487,409 British samples reveals fine-scale evolutionary history and trait associations. *Nature Communications*, in press.
