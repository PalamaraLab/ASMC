# ASMC Release Notes

## v1.2 (2021-08-??)

All functionality for ASMC and FastSMC is now in [this repository](link???).




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