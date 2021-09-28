# FastSMC

- [Running FastSMC](#running-fastsmc)
	- [Detailed command line options](#detailed-command-line-options)
	- [Input file formats](#input-file-formats)
	- [Output format](#output-format)
	- [Binary output](#binary-output)
	- [Example using the compiled FastSMC executable (C++)](#example-using-the-compiled-fastsmc-executable-c)
- [Relationship to ASMC](#relationship-to-asmc)

The Fast Sequentially Markovian Coalescent (FastSMC) algorithm is an extension to the ASMC algorithm, adding an identification step by hashing (currently using an improved version of the GERMLINE algorithm).
FastSMC is an accurate method to detect Identical-By-Descent segments which enables estimating the time to most recent common ancestor for IBD individuals, and provides an estimate of uncertainty for detected IBD regions.

**This document is not intended as an extensive guide, a more detailed user manual is under development, data and annotations from the FastSMC paper can be found [here](https://palamaralab.github.io/software/fastsmc/).**

## Running FastSMC

You can run FastSMC as a C++ compiled executable or using [Python bindings](./fastsmc_python.md).
Follow this [quickstart guide](./quickstart_user.md) for information on compiling the FastSMC targets.

### Detailed command line options
See the [ASMC documentation](./asmc.md) for parameters related to the validation step.
Additional parameters related to the identification step are listed below.
Note: default parameter values are likely to change in future versions.

```
  --inFileRoot                	Prefix of input files (.hap, .samples, .map).
                              	[mandatory]
  --decodingQuantFile         	Decoding quantities file.
                              	[mandatory]
  --outFileRoot               	Prefix of output file.
                              	[mandatory]
  --hashing                  	Use of hashing to pre-process IBD segments. If off, no identification step will be performed.
                              	[default 1/on]
  --min_m arg (=1)		        Minimum match length (in cM).
				                [default = 1.0]
  --time arg (=100)		        Time threshold to define IBD in number of generations.
				                [default = 100]
  --skip arg (=0)		        Skip words with (seeds/samples) less than this value
				                [default 0.0]
  --min_maf arg (=0)		    Minimum minor allele frequency
				                [default 0.0]
  --gap arg (=1)		        Allowed gaps
                                [default 1]
  --max_seeds arg (=0)		    Dynamic hash seed cutoff
				                [default 0/off]
  --recall arg (=3)		        Recall level from 0 to 3 (higher value means higher recall).
				                [default = 3]
  --segmentLength		        Output length in centimorgans of each IBD segment.
				                [default 1/on]
  --perPairMAP			        Output MAP age estimate for each IBD segment.
				                [default 1/on]
  --perPairPosteriorMeans	    Output posterior mean age estimate for each IBD segment.
				                [default 1/on]
  --noConditionalAgeEstimates	Do not condition the age estimates on the TMRCA being between present time and t 
				                generations ago (where t is the time threshold).
				                [default 0/off]
  --bin				            Binary output
				                [default off]
  --batchSize			        Size of batches to be decoded.
				                [default = 32]
  --hashingOnly                 Only perform GERMLINE2 hashing, not ASMC decoding
  				                [default off]
```

> Note: the `hashingOnly` flag has not been extensively tested.
You may also want to look into [this repository](https://github.com/gusevlab/germline2) for a standalone version.

Suggested optimal parameters for IBD detection within the past 25, 50, 100, 150 and 200 generations are provided in the FastSMC paper.

### Input file formats

Input files are provided to FastSMC with the --inFileRoot option. You may want to look at files in `ASMC_data` for examples of the file formats described below.

#### Phased haplotypes in Oxford haps/sample format (.hap/.hap.gz, .samples)
These files are provided in input to FastSMC. The file format explained [here](https://www.cog-genomics.org/plink/2.0/formats#haps). These files are output by phasing programs like Eagle and Shapeit.

#### Genetic map (.map)
The genetic map file needs to provide physical positions (in base pairs) and genetic positions (in centimorgans).
FastSMC expects a tab separated file containing at least 3 columns, with physical positions in the first column and genetic positions in the third column.
The second column may contain any values.
The file may be optionally compressed using gzip.

#### Decoding quantities (.decodingQuantities.gz)
See the instructions above to generate decoding quantities files and the [ASMC manual](./asmc.md#decoding-quantities-decodingquantitiesgz) for more details.

Note: the CEU.demo demographic model and the decoding quantities for CEU+UKBB previously provided in [this repository](https://github.com/PalamaraLab/FastSMC) and [this repository](https://github.com/PalamaraLab/ASMC_legacy) were mistakenly encoded as diploid rather than haploid.
The file [CEU.demo](https://github.com/PalamaraLab/ASMC_data/tree/main/demographies) and CEU+UKBB decoding quantities [here](https://github.com/PalamaraLab/ASMC_data/tree/main/decoding_quantities) have now been fixed.
They were generated using v2.2.1 of the [PrepareDecoding tool](https://github.com/PalamaraLab/PrepareDecoding/releases/tag/v2.2.1), which also provides a simpler interface for computing decoding quantities as well as support for additional demographic models.
Using these new decoding quantities with v1.2 of ASMC will tend to produce more recent estimates for TMRCAs compared to the decoding quantities distributed with v1.0 and v1.1.
This should not have a substantial impact on most downstream analyses.

### Output format

FastSMC generates an .ibd.gz (or .bibd.gz if binary output) file in the specified location.
Each line corresponds to a pairwise shared segment, with the following fields:

	0. First individual's family identifier
	1. First individual identifier
	2. First individual haplotype identifier (1 or 2)
	3. Second individual's family identifier
	4. Second individual identifier
	5. Second individual haplotype identifier (1 or 2)
	6. Chromosome number
	7. Starting position of the IBD segment (inclusive)
	8. Ending position of the IBD segment (inclusive)
	9. (optional) Length in centimorgans of IBD segment
	10. IBD score
	11. (optional) Average mean posterior age estimate of the IBD segment
	12. (optional) Average MAP age estimate of the IBD segment

### Binary output

If you use the --bin option, FastSMC will generate a compressed binary (.bib.gz) output. This can be then converted to text format using the BinaryDataReader class in Python (see notebooks for an example) and using the convertBinary executable in C++ (see the C++ example below).

### Example using the compiled FastSMC executable (C++)

Following the compilation instructions in the [quickstart guide](./quickstart_user.md) will create an executable

```
build/FastSMC_exe
```

which can be used by providing command line arguments summarised above.
For an example of IBD detection within the past 50 generations, please run the following command line:

```bash
sh cpp_example/FastSMC_example.sh
```

To parallelise this run over 4 independent jobs, you can run the following command line instead:

```bash
sh cpp_example/FastSMC_example_multiple_jobs.sh
```

This example will run multiple jobs in different threads on the same machine.
If you are running FastSMC on a cluster then it may be more appropriate to instead use the job scheduler such as `qsub`.

A binary output file will be generated and then converted to text format using the convertBinary executable. The first 10 lines will be printed.

Either way of running FastSMC (Python bindings or C++) will run it on a simulated dataset as described in the FastSMC paper.
An output file with IBD segments will be generated (in notebooks/ or c++\_example/ respectively), and run time should be less than 4s.

## Relationship to ASMC

The Ascertained Sequentially Markovian Coalescent is a method to efficiently estimate pairwise coalescence time along the genome.
It can be run using SNP array or whole-genome sequencing (WGS) data.

FastSMC builds on ASMC, and this repository can be used to run ASMC analysis.
A user manual can be found [here](./asmc.md) and data and annotations from the ASMC paper can be found [here](https://palamaralab.github.io/software/asmc/data/).

We don't currently provide scripts to automate some of the analyses described in the FastSMC paper.
Some thoughts on the selection analysis may be found [here](./asmc.md#density-of-recent-coalescence-drc-statistic).
