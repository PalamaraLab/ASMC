# Quickstart guide for developers

- [Linux](#linux)
- [macOS](#macos)
- [ResComp (oxford research computing)](#rescomp-oxford-research-computing)
- [Python bindings](#python-bindings)

## Linux

This guide assumes you have a C++17 compatible compiler (e.g. gcc >= 8.3 or clang >= 7) and [CMake >= 3.15](https://cmake.org/install/).
Additionally, to compile the Python bindings you need Python with development files:

```bash
sudo apt install python3-dev
```

Then, follow these steps:

```bash
# Get the source
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC_dev
cd ASMC_dev

# Create a build directory
mkdir build && cd build

# Configure and build
# On first run, CMake will build the required dependencies
cmake ..
cmake --build . --parallel 4
```

## macOS

This guide assumes you have a recent version of the [Xcode command line tools](https://developer.apple.com/xcode/features/) and [Homebrew](https://brew.sh/).
Install the following dependencies:

```bash
brew install cmake
brew install libomp
brew install python  # for python bindings, if required
```

Then, follow these steps:

```bash
# Get the source
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC_dev
cd ASMC_dev

# Create a build directory
mkdir build && cd build

# Configure and build
# On first run, CMake will build the required dependencies
cmake ..
cmake --build . --parallel 4
```

## ResComp (oxford research computing)

All necessary dependencies are already installed on ResComp. Simply follow these steps:

```bash
# Load required modules
module load GCC/10.2.0
module load CMake/3.18.4-GCCcore-10.2.0
module load git/2.28.0-GCCcore-10.2.0-nodocs
module load Python/3.8.6-GCCcore-10.2.0

# Get the source
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC_dev
cd ASMC_dev

# Create a build directory
mkdir build && cd build

# Configure and build
# On first run, CMake will build the required dependencies
cmake ..
cmake --build . --parallel 4
```

## Python bindings

These instructions are platform independent, assuming you have installed all dependencies (excluding those from vcpkg) according to the instructions above.
From the `ASMC_dev` directory:

```bash
python3 -m venv venv
source venv/bin/activate

pip install --upgrade pip setuptools wheel ninja
pip install .
```
