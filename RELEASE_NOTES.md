# ASMC Release Notes

## v1.4.0 (2025-10-22)

### Breaking changes

- Dependencies are no longer managed with vcpkg. Boost and zlib should now be obtained from a system package manager, while other dependencies are fetched using CMake's FetchContent during configuration.

### Other changes

- Added support for cross-platform SIMD vectorization using the Google Highway library. This improves performance by dynamically dispatching to the most powerful supported instruction set at runtime (e.g., AVX-512, AVX2, or NEON), ensuring optimal performance on different CPU architectures.

- Python wheels are now available for Linux and macOS on both x86_64 and arm64/AArch64 architectures, for CPython versions 3.9 to 3.14 inclusive.

## v1.3.1 (2023-06-30)

### Breaking changes

None

### Other changes

- The location of a `.map` or `.map.gz` file can now be optionally specified explicitly: previously it was assumed to be at the `inFileRoot`.


## v1.3 (2023-03-03)

### Breaking changes

None

### Other changes

- Decoding a batch can now be done in a selected subregion with from / to indices.
  A `cm_burn_in` parameter takes into account additional variants on either side of the subregion for HMM burn-in.
- Allow the user to access selected attributes of the DecodingParams and Data from the ASMC object.
- Python continuous integration now uses Python 3.8 and 3.11 (previously 3.6 and 3.9)
- Update Catch to v2.13.


## v1.2 (2021-09-28)

All functionality for ASMC and FastSMC is now in this repository ([link](https://github.com/PalamaraLab/ASMC)).

### Breaking changes

- Fixed an issue with demographic models.
  The `CEU.demo` demographic model and the decoding quantities for CEU+UKBB previously provided in the repository were mistakenly encoded as diploid rather than haploid. 
  CEU.demo and CEU+UKBB decoding quantities have now been updated and can be found in [this repository](https://github.com/PalamaraLab/ASMC_data).
  Also see the manual for a note on how this affects analyses.

### Other changes

- New API for decoding pairs with ASMC.
  In addition to running full analyses as described in the ASMC paper, users can now decode specific pairs and get back a variety of summary statistics.
  See the [ASMC python documentation](https://github.com/PalamaraLab/ASMC/blob/main/docs/asmc_python.md) for details.
- New, more extensive, [documentation](https://github.com/PalamaraLab/ASMC/blob/main/docs/) is available.


## v1.1 (2021-01-20)

[Legacy repository](https://github.com/PalamaraLab/FastSMC/releases/tag/v1.1)

Improvements to documentation and default use.
No changes to any core functionality.

### Breaking changes

- The hashing functionality, previously named `GERMLINE`, has been renamed to `hashing`.
  This includes the command line flag for turning this behaviour on/off, which is now `--hashing`.

### Other changes

- `--hashing` is now ON by default when running the FastSMC executable: previously, `--GERMLINE` was OFF by default.
- Extra output, including the IBD segment length, posterior mean, and MAP, are now on by default.
  This behaviour can be toggled with the flags `--segmentLength`, `--perPairPosteriorMeans`, `--perPairMAP`.
- An example script has been added to `cpp_example/FastSMC_example_multiple_jobs.sh` that demonstrates how to run FastSMC with multiple jobs simultaneously.
- The README has been updated to focus on FastSMC functionality.
- More robust checking is now used to verify the decoding quantities file is correct before reading it.
- CMake will now, by default, build in Release mode (giving 03 optimisation on Linux).
  Previously, Debug was used by default.


## v1.0 (2020-09-18)

[Legacy repository](https://github.com/PalamaraLab/FastSMC/releases/tag/v1.0)

First public release of FastSMC, with functionality as described and used in [this paper](https://doi.org/10.1038/s41467-020-19588-x).
