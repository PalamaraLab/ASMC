# ASMC

- [Summary (TL;DR)](#summary-tldr)
- [Running ASMC](#running-asmc)
- [A complete example](#a-complete-example)
- [Detailed command line options](#detailed-command-line-options)
- [Input/output file formats](#inputoutput-file-formats)

This page describes compiling and running the ASMC program described in [this paper](https://doi.org/10.1038/s41588-018-0177-x):

- P. Palamara, J. Terhorst, Y. Song, A. Price. **High-throughput inference of pairwise coalescence times identifies signals of selection and enriched disease heritability.** *Nature Genetics*, 2018.


A page with data and annotations from the paper can be found [here](https://palamaralab.github.io/software/asmc/data/).

## Summary (TL;DR)

The Ascertained Sequentially Markovian Coalescent is a method to efficiently estimate pairwise coalescence times along the genome.
It can be run using SNP array or whole-genome sequencing (WGS) data.

To run ASMC, download it [here](https://github.com/PalamaraLab/ASMC), and follow the [quickstart guide](./quickstart_user.md).
Compile the `ASMC_exe` target to build the ASMC executable.
To compute pairwise coalescence times for the following files containing SNP array data for 150 phased diploid samples:

- `ASMC_data/examples/asmc/exampleFile.n300.array.hap.gz`
- `ASMC_data/examples/asmc/exampleFile.n300.array.samples`
- `ASMC_data/examples/asmc/exampleFile.n300.array.map.gz`

(see [here](#inputoutput-file-formats) for file formats), you can run the following ASMC command:
```bash
./build/ASMC_exe \
        --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
        --inFileRoot ASMC_data/examples/asmc/exampleFile.n300.array \
        --posteriorSums
```

This will generate `ASMC_data/examples/asmc/exampleFile.n300.array.1-1.sumOverPairs.gz`, which contains a matrix of size SxD, where S is the number of sites in the data, and D is the number of discrete coalescence time intervals defined in `FILES/DISC/30-100-2000.disc`.
The [*i*,*j*]-th entry of this matrix contains the sum of posterior coalescence probabilities for all samples at SNP *i* and time *j*.

### Running ASMC

Once you have computed or downloaded the decoding quantities corresponding to your demographic model, time discretization, and SNP allele frequencies, you can analyze SNP or WGS data using the `ASMC` program:

```
./ASMC \
        --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
        --inFileRoot ASMC_data/examples/asmc/exampleFile.n300.array \
        --majorMinorPosteriorSums
```

For WGS data, add the `--mode sequence` option:

```
./ASMC \
        --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
        --inFileRoot ASMC_data/examples/asmc/exampleFile.n300 \
        --majorMinorPosteriorSums \
        --mode sequence
```

Using the `--majorMinorPosteriorSums` flag, ASMC will output the sum of posterior coalescence probabilities for all analyzed pairs of individuals.
These will be written in `ASMC_data/examples/asmc/exampleFile.n300.array.{00,01,11}.sumOverPairs.gz` for the SNP array example, and `ASMC_data/examples/asmc/exampleFile.n300.{00,01,11}.sumOverPairs.gz` for the WGS example.

If you are decoding a large number of samples, you can break down the computation in several *jobs* using the `--jobs int` and `--jobInd int` flags, which take the total number of jobs to be performed and the current job index as arguments.
If you don't specify a name for the output files, a job index will be automatically added to the default output path.

The full set of command line options for `ASMC` is as follows:

```
Mandatory:
  --inFileRoot arg              Prefix of hap|haps|hap.gz|haps.gz and sample|samples file
  --decodingQuantFile arg       Decoding quantities file

Choose one of:
  --posteriorSums               Output file for sum of posterior distribution
                                over pairs.
  --majorMinorPosteriorSums     Output file for sum of posterior distribution 
                                over pairs, partitioned by major/minor alleles.                                

Optional:
  --outFileRoot arg             Output file for sum of posterior distribution
                                over pairs (default: --hapsFileRoot argument)
  --jobs int (=0)               Number of jobs being done in parallel
  --jobInd int (=0)             Job index (1..jobs)
  --mode string (=array)        Decoding mode. Choose from {sequence, array}.
  --compress (=false)           Compress emission to binary (no CSFS)
  --useAncestral (=false)       Assume ancestral alleles are coded as 1 in
                                input (will assume 1 = minor otherwise)
  --skipCSFSdistance int (=0)   Genetic distance between two CSFS emissions
```

In addition to the arguments described above, `ASMC` options include:
- `--compress` is a shorthand for `--skipCSFSdistance Infinity` (see below).
- `--useAncestral` can be used to specify that a `1` in the data specifies an ancestral allele. This will cause the CSFS to be used without folding. This is mostly not needed.
- `--skipCSFSdistance int`, which takes an integer argument, specifies the minimum distance for a CSFS emission to be used. The default is `0` (always use CSFS). Setting `--skipCSFSdistance Infinity` (which is the same as `--compress`) leads to never using the CSFS (i.e. the classic PSMC emission if decoding WGS data, or a binary emission which controls for ascertainment if decoding SNP array data).

## Input/output file formats

You may want to look at files in FILES/\* for examples of the file formats described below.

### Phased haplotypes in Oxford haps/sample format (\*.hap/hap.gz, samples)

[File format specification.](https://www.cog-genomics.org/plink/2.0/formats#haps)

These files are provided in input to `ASMC` and optionally `ASMCprepareDecoding`. The file format explained [here](https://www.cog-genomics.org/plink/2.0/formats#haps). These files are output by phasing programs like Eagle and Shapeit.

### Genetic map (\*.map/map.gz)

The genetic map provided in input to `ASMC` is in [Plink map](https://www.cog-genomics.org/plink2/formats#map) format, in which each line has four columns with format "Chromosome SNPName GeneticPosition PhysicalPosition". Genetic positions are in centimorgans, physical positions are in bp. The map can be optionally compressed using gzip.

### Demographic history (\*.demo)

The demographic history provided in input to `ASMCprepareDecoding` represents a piece-wise constant history of past effective population sizes, with format
```
TimeStart   PopulationSize
```

Where TimeStart is the first generation where the population has size PopulationSize. Note that population size is *haploid*, and that the demographic model is usually built assuming a specific mutation rate, which is passed as an argument to the `ASMCprepareDecoding` program. The first line should contain generation `0`. You can obtain this model using e.g. PSMC/MSMC/SMC++. If your model is not piecewise constant, you will need to approximate it as piecewise constant. The last provided interval is assumed to last until time=Infinity (and is usually remote enough to have negligible effects on the results).

### Time discretization (\*.disc)

The list of discrete time intervals provided in input to `ASMC` contains a single number per line, representing time measured in (continuous) generations, and starting at generation `0.0`. For instance, the list `FILES/DISC/10.disc` contains 10 time intervals:
```
0.0
1118.2
1472.2
1849.7
2497.0
3963.8
9120.8
15832.9
24139.9
34891.6
```

The intervals defined in this file are: `{0.0-1118.2, 1118.2-1472.2, 1472.2-1849.7, 1849.7-2497.0, 2497.0-3963.8, 3963.8-9120.8, 9120.8-15832.9, 15832.9-24139.9, 24139.9-34891.6, 34891.6-Infinity}`.

### Decoding quantities (\*.decodingQuantities.gz)

The `*.decodingQuantities.gz` file is generated by `ASMCprepareDecoding` and input into `ASMC`. It is used to perform efficient inference of pairwise coalescence times. There is no need to understand the content of this file.

### Time discretization intervals (\*.intervalsInfo)

The `*.intervalsInfo` file is generated by the `ASMCprepareDecoding` and input into `ASMC`. It contains some useful information about the time discretization and the demographic model. It contains a number of lines corresponding to the number of discrete time intervals used in the analysis. Each line has format:

```
IntervalStart   ExpectedCoalescenceTime IntervalEnd
```

IntervalStart and IntervalEnd represent the start/end of each discrete time interval, ExpectedCoalescenceTime is the expected coalescence time for a pair of individuals who have been inferred to coalesce within this time interval, and depends on the demographic model.

### Sum of pairwise posterior coalescence probabilities `*.sumOverPairs.gz`

The output of the `ASMC` analysis is written in `*.{00,01,11}.sumOverPairs.gz` files. Each file contains a matrix of size SxD, where S is the number of sites in the data, and D is the number of discrete time intervals used in the analysis. The [*i*,*j*]-th entry of each matrix contains the sum of posterior coalescence probabilities for all samples at SNP *i* and discrete coalescence time *j*. The output breaks down coalescence events of samples carrying different alleles at each site, using the `{00,01,11}` suffixes. Specifically:

- The *i*-th row of the matrix in `*.00.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that are homozygous `0` at site *i*.
- The *i*-th row of the matrix in `*.01.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that heterozygous at site *i*.
- The *i*-th row of the matrix in `*.11.sumOverPairs.gz` contains the sum of posterior probabilities for all pairs of samples that are homozygous `1` at site *i*.
