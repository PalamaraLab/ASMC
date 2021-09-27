# Quickstart guide for users

- [Python bindings](#python-bindings)
- [Linux](#linux)
- [macOS](#macos)
- [Without vcpkg](#without-vcpkg)

## Python bindings

If you want to use ASMC or FastSMC via their Python interface, you can simply install ASMC using pip:

```
pip install asmc-asmc
```

For examples, see the [ASMC python documentation](./asmc_python.md) and [FastSMC python documentation](./fastsmc_python.md).

If you want to compile the C++ executables, read on.

## Linux

This guide assumes you have a C++17 compatible compiler (e.g. gcc >= 8.3 or clang >= 7) and [CMake >= 3.15](https://cmake.org/install/).
Then, follow these steps to build the ASMC, FastSMC and binary conversion executables:

```bash
# Get the source
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC
cd ASMC

# Create a build directory
mkdir build && cd build

# Configure and build
# On first run, CMake will build the required dependencies
cmake -DASMC_NO_PYTHON=TRUE ..
cmake --build . --parallel 4
```

## macOS

This guide assumes you have a recent version of the [Xcode command line tools](https://developer.apple.com/xcode/features/) and [Homebrew](https://brew.sh/).
Install the following dependencies:

```bash
brew install cmake
brew install libomp
```

Then, follow these steps to build the ASMC, FastSMC and binary conversion executables:

```bash
# Get the source
git clone --recurse-submodules https://github.com/PalamaraLab/ASMC
cd ASMC

# Create a build directory
mkdir build && cd build

# Configure and build
# On first run, CMake will build the required dependencies
cmake -DASMC_NO_PYTHON=TRUE ..
cmake --build . --parallel 4
```

## Without vcpkg

If you would like to compile ASMC without using [vcpkg](https://github.com/microsoft/vcpkg/) to handle dependencies, you should first ensure all dependencies are installed:

**Ubuntu**
```bash
sudo apt install libboost-iostreams-dev libboost-math-dev libboost-program-options-dev libeigen3-dev libfmt-dev librange-v3-dev zlib1g-dev
```

**macOS**
```bash
brew install boost eigen fmt range-v3 zlib
````

Then, when you run CMake, add the following definition:

```bash
cmake -DASMC_AVOID_VCPKG=true ..
```

You may additionally choose to not recursively clone all submodules, as long as you still obtain the `DataModule` submodule.
From the ASMC directory:

```bash
git clone https://github.com/PalamaraLab/ASMC
cd ASMC
git submodule update --init DataModule
```
