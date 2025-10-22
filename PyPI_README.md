[![Unit tests: Ubuntu](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-unit.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-unit.yml)
[![Unit tests: macOS](https://github.com/PalamaraLab/ASMC/actions/workflows/macos-unit.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/macos-unit.yml)
[![Python 3.8 3.11](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-python.yml/badge.svg)](https://github.com/PalamaraLab/ASMC/actions/workflows/ubuntu-python.yml)
[![Regression test](https://github.com/PalamaraLab/ASMC/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![Ubuntu asan](https://github.com/PalamaraLab/ASMC/workflows/Ubuntu%20asan/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![Ubuntu no sse/avx](https://github.com/PalamaraLab/ASMC/workflows/Ubuntu%20no%20sse/avx/badge.svg)](https://github.com/PalamaraLab/ASMC/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/ASMC/branch/main/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/ASMC)

# ASMC and FastSMC

This repository provides ASMC and its extension FastSMC, implemented in C++ with Python bindings.
Prebuilt CPython wheels are available for Linux (compatible with glibc ≥ 2.28) and macOS (built on macOS 15 for x86_64 and macOS 14 for arm64).

| Platform \ CPython          | ≤3.8 | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|-----------------------------| ---- | --- | ---- | ---- | ---- | ---- | ---- |
| Linux x86_64                | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| Linux aarch64               | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Intel (x86_64)        | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Apple Silicon (arm64) | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |

## Quickstart

### Install the Python module from PyPI

Most functionality is available through a Python module which can be installed with:

```bash
pip install asmc-asmc
```

### Documentation

The following pages of documentation contains specific information:
- [Quickstart guide for users](https://github.com/PalamaraLab/ASMC/blob/main/docs/quickstart_user.md)
- [ASMC python docs](https://github.com/PalamaraLab/ASMC/blob/main/docs/asmc_python.md)
- [FastSMC python docs](https://github.com/PalamaraLab/ASMC/blob/main/docs/fastsmc_python.md)

This Python module is currently available on Linux and macOS.

Example Jupyter notebooks showcasing basic functionality can be found here:
- [Example notebooks](https://github.com/PalamaraLab/ASMC/tree/main/notebooks)

## License

ASMC and FastSMC are distributed under the GNU General Public License v3.0 (GPLv3). For any questions or comments on ASMC, please contact Pier Palamara using `<lastname>@stats.ox.ac.uk`.

## Reference

If you use this software, please cite the appropriate reference(s) below.

The ASMC algorithm and software were developed in
- P. Palamara, J. Terhorst, Y. Song, A. Price. High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability. *Nature Genetics*, 2018.

The FastSMC algorithm and software were developed in
- J. Nait Saada, G. Kalantzis, D. Shyr, F. Cooper, M. Robinson, A. Gusev, P. F. Palamara. Identity-by-descent detection across 487,409 British samples reveals fine-scale evolutionary history and trait associations. *Nature Communications*, 2020.
