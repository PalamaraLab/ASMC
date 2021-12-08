ASMC
====

-  `Summary (TL;DR) <#summary-tldr>`__

   -  `Running ASMC <#running-asmc>`__

-  `Input/output file formats <#inputoutput-file-formats>`__

   -  `Phased haplotypes in Oxford haps/sample format (*.hap/hap.gz,
      samples) <#phased-haplotypes-in-oxford-hapssample-format-haphapgz-samples>`__
   -  `Genetic map (*.map/map.gz) <#genetic-map-mapmapgz>`__
   -  `Demographic history (*.demo) <#demographic-history-demo>`__
   -  `Time discretization (*.disc) <#time-discretization-disc>`__
   -  `Decoding quantities
      (*.decodingQuantities.gz) <#decoding-quantities-decodingquantitiesgz>`__
   -  `Time discretization intervals
      (*.intervalsInfo) <#time-discretization-intervals-intervalsinfo>`__
   -  `Sum of pairwise posterior coalescence probabilities
      ``*.sumOverPairs.gz`` <#sum-of-pairwise-posterior-coalescence-probabilities-sumoverpairsgz>`__

-  `Tools, scripts, and analyses <#tools-scripts-and-analyses>`__

   -  `Tool to merge output of parallel ASMC
      jobs <#tool-to-merge-output-of-parallel-asmc-jobs>`__
   -  `Script to visualize average coalescence density in a
      region <#script-to-visualize-average-coalescence-density-in-a-region>`__
   -  `Density of recent coalescence (DRC)
      statistic <#density-of-recent-coalescence-drc-statistic>`__

-  `Precomputed decoding
   quantities <#precomputed-decoding-quantities>`__

This page describes compiling and running the ASMC program described in
`this paper <https://doi.org/10.1038/s41588-018-0177-x>`__:

-  P. Palamara, J. Terhorst, Y. Song, A. Price. **High-throughput
   inference of pairwise coalescence times identifies signals of
   selection and enriched disease heritability.** *Nature Genetics*,
   2018.

A page with data and annotations from the paper can be found
`here <https://palamaralab.github.io/software/asmc/data/>`__.

Summary (TL;DR)
---------------

The Ascertained Sequentially Markovian Coalescent is a method to
efficiently estimate pairwise coalescence times along the genome. It can
be run using SNP array or whole-genome sequencing (WGS) data.

To run ASMC, download it `here <https://github.com/PalamaraLab/ASMC>`__,
and follow the `quickstart guide <./quickstart_user.md>`__. Compile the
``ASMC_exe`` target to build the ASMC executable. To compute pairwise
coalescence times for the following files containing SNP array data for
150 phased diploid samples:

-  ``ASMC_data/examples/asmc/exampleFile.n300.array.hap.gz``
-  ``ASMC_data/examples/asmc/exampleFile.n300.array.samples``
-  ``ASMC_data/examples/asmc/exampleFile.n300.array.map.gz``

(see `here <#inputoutput-file-formats>`__ for file formats), you can run
the following ASMC command:

.. code:: bash

   ./build/ASMC_exe \
           --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
           --inFileRoot ASMC_data/examples/asmc/exampleFile.n300.array \
           --posteriorSums

This will generate
``ASMC_data/examples/asmc/exampleFile.n300.array.1-1.sumOverPairs.gz``,
which contains a matrix of size SxD, where S is the number of sites in
the data, and D is the number of discrete coalescence time intervals
defined in ``ASMC_data/discretizations/30-100-2000_CEU.disc``. The
[*i*,\ *j*]-th entry of this matrix contains the sum of posterior
coalescence probabilities for all samples at SNP *i* and time *j*.

Running ASMC
~~~~~~~~~~~~

Once you have computed or downloaded the decoding quantities
corresponding to your demographic model, time discretization, and SNP
allele frequencies, you can analyze SNP or WGS data using the ``ASMC``
program:

::

   ./build/ASMC_exe \
           --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
           --inFileRoot ASMC_data/examples/asmc/exampleFile.n300.array \
           --majorMinorPosteriorSums

For WGS data, add the ``--mode sequence`` option:

::

   ./build/ASMC_exe \
           --decodingQuantFile ASMC_data/decoding_quantities/30-100-2000_CEU.decodingQuantities.gz \
           --inFileRoot ASMC_data/examples/asmc/exampleFile.n300 \
           --majorMinorPosteriorSums \
           --mode sequence

Using the ``--majorMinorPosteriorSums`` flag, ASMC will output the sum
of posterior coalescence probabilities for all analyzed pairs of
individuals. These will be written in
``ASMC_data/examples/asmc/exampleFile.n300.array.{00,01,11}.sumOverPairs.gz``
for the SNP array example, and
``ASMC_data/examples/asmc/exampleFile.n300.{00,01,11}.sumOverPairs.gz``
for the WGS example.

If you are decoding a large number of samples, you can break down the
computation in several *jobs* using the ``--jobs int`` and
``--jobInd int`` flags, which take the total number of jobs to be
performed and the current job index as arguments. If you don't specify a
name for the output files, a job index will be automatically added to
the default output path.

The full set of command line options for ``ASMC`` is as follows:

::

   Mandatory:
     --inFileRoot arg                Prefix of hap|haps|hap.gz|haps.gz and sample|samples file
     --decodingQuantFile arg         Decoding quantities file

   Choose one of:
     --posteriorSums                 Output file for sum of posterior distribution
                                     over pairs.
     --majorMinorPosteriorSums       Output file for sum of posterior distribution 
                                     over pairs, partitioned by major/minor alleles.                                

   Optional:
     --outFileRoot arg               Output file for sum of posterior distribution
                                     over pairs (default: --hapsFileRoot argument)
     --jobs int (=0)                 Number of jobs being done in parallel
     --jobInd int (=0)               Job index (1..jobs)
     --mode string (=array)          Decoding mode. Choose from {sequence, array}.
     --compress (=false)             Compress emission to binary (no CSFS)
     --useAncestral (=false)         Assume ancestral alleles are coded as 1 in
                                     input (will assume 1 = minor otherwise)
     --skipCSFSdistance float (=0.0) Genetic distance (in cM) between two CSFS emissions

In addition to the arguments described above, ``ASMC`` options include:

-  ``--compress`` is a shorthand for ``--skipCSFSdistance Infinity``
   (see below).
-  ``--useAncestral`` can be used to specify that a ``1`` in the data
   specifies an ancestral allele. This will cause the CSFS to be used
   without folding. This is mostly not needed.
-  ``--skipCSFSdistance float``, which takes a floating point argument,
   specifies the minimum distance for a CSFS emission to be used. The
   default is ``0.0`` (always use CSFS). Setting
   ``--skipCSFSdistance Infinity`` (which is the same as ``--compress``)
   leads to never using the CSFS (i.e. the classic PSMC emission if
   decoding WGS data, or a binary emission which controls for
   ascertainment if decoding SNP array data).

Input/output file formats
-------------------------

You may want to look at files in ``ASMC_data/`` for examples of the file
formats described below.

.. _phased-haplotypes-in-oxford-hapssample-format-haphapgz-samples:

Phased haplotypes in Oxford haps/sample format (*.hap/hap.gz, samples)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`File format
specification. <https://www.cog-genomics.org/plink/2.0/formats#haps>`__

These files are provided in input to ``ASMC`` and optionally
``ASMCprepareDecoding``. The file format explained
`here <https://www.cog-genomics.org/plink/2.0/formats#haps>`__. These
files are output by phasing programs like Eagle and Shapeit.

.. _genetic-map-mapmapgz:

Genetic map (*.map/map.gz)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The genetic map provided in input to ``ASMC`` is in `Plink
map <https://www.cog-genomics.org/plink2/formats#map>`__ format, in
which each line has four columns with format "Chromosome SNPName
GeneticPosition PhysicalPosition". Genetic positions are in
centimorgans, physical positions are in bp. The map can be optionally
compressed using gzip.

.. _demographic-history-demo:

Demographic history (*.demo)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The demographic history provided in input to ``ASMCprepareDecoding``
represents a piece-wise constant history of past effective population
sizes, with format

::

   TimeStart   PopulationSize

Where TimeStart is the first generation where the population has size
PopulationSize. Note that population size is *haploid*, and that the
demographic model is usually built assuming a specific mutation rate,
which is passed as an argument to the ``ASMCprepareDecoding`` program.
The first line should contain generation ``0``. You can obtain this
model using e.g. PSMC/MSMC/SMC++. If your model is not piecewise
constant, you will need to approximate it as piecewise constant. The
last provided interval is assumed to last until time=Infinity (and is
usually remote enough to have negligible effects on the results).

The demographic models used with ASMC can be found
`here <https://github.com/PalamaraLab/ASMC_data/tree/main/demographies>`__
and were inferred using smc++ in the following paper:

   Spence, J.P. and Song, Y.S. **Inference and analysis of
   population-specific fine-scale recombination maps across 26 diverse
   human populations.** *Science Advances*, Vol. 5, No. 10, eaaw9206
   (2019), [`doi <https://doi.org/10.1126/sciadv.aaw9206>`__].

They correspond to `these population
sizes <https://github.com/popgenmethods/pyrho/blob/master/smcpp_popsizes_1kg.csv>`__,
but rescaled to assume mutation rate of 1.65e-8.

.. _time-discretization-disc:

Time discretization (*.disc)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The list of discrete time intervals provided in input to ``ASMC``
contains a single number per line, representing time measured in
(continuous) generations, and starting at generation ``0.0``. For
instance, the list ``ASMC_data/discretizations/30-100-2000_CEU.disc``
contains time intervals:

::

   0.0
   30.0
   60.0
   90.0
   120.0
   150.0
   180.0
   210.0
   240.0
   270.0
   300.0
   330.0
   360.0
   ... <lines omitted>
   79855.6
   96263.0
   124311.7

The intervals defined in this file are:
``{0.0-30.0, 30.0-60.0, ..., 96263.0-124311.7, 124311.7-Infinity}``.

.. _decoding-quantities-decodingquantitiesgz:

Decoding quantities (*.decodingQuantities.gz)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``*.decodingQuantities.gz`` file is generated by
``ASMCprepareDecoding`` and input into ``ASMC``. It is used to perform
efficient inference of pairwise coalescence times. There is no need to
understand the content of this file.

Note: the CEU.demo demographic model and the decoding quantities for
CEU+UKBB previously provided in `this
repository <https://github.com/PalamaraLab/FastSMC>`__ and `this
repository <https://github.com/PalamaraLab/ASMC_legacy>`__ were
mistakenly encoded as diploid rather than haploid. The file
`CEU.demo <https://github.com/PalamaraLab/ASMC_data/tree/main/demographies>`__
and CEU+UKBB decoding quantities
`here <https://github.com/PalamaraLab/ASMC_data/tree/main/decoding_quantities>`__
have now been fixed. They were generated using v2.2.1 of the
`PrepareDecoding
tool <https://github.com/PalamaraLab/PrepareDecoding/releases/tag/v2.2.1>`__,
which also provides a simpler interface for computing decoding
quantities as well as support for additional demographic models. Using
these new decoding quantities with v1.2 of ASMC will tend to produce
more recent estimates for TMRCAs compared to the decoding quantities
distributed with v1.0 and v1.1. This should not have a substantial
impact on most downstream analyses.

.. _time-discretization-intervals-intervalsinfo:

Time discretization intervals (*.intervalsInfo)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``*.intervalsInfo`` file is generated by the ``ASMCprepareDecoding``
and input into ``ASMC``. It contains some useful information about the
time discretization and the demographic model. It contains a number of
lines corresponding to the number of discrete time intervals used in the
analysis. Each line has format:

::

   IntervalStart   ExpectedCoalescenceTime IntervalEnd

IntervalStart and IntervalEnd represent the start/end of each discrete
time interval, ExpectedCoalescenceTime is the expected coalescence time
for a pair of individuals who have been inferred to coalesce within this
time interval, and depends on the demographic model.

.. _sum-of-pairwise-posterior-coalescence-probabilities-sumoverpairsgz:

Sum of pairwise posterior coalescence probabilities ``*.sumOverPairs.gz``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output of the ``ASMC`` analysis is written in
``*.{00,01,11}.sumOverPairs.gz`` files. Each file contains a matrix of
size SxD, where S is the number of sites in the data, and D is the
number of discrete time intervals used in the analysis. The
[*i*,\ *j*]-th entry of each matrix contains the sum of posterior
coalescence probabilities for all samples at SNP *i* and discrete
coalescence time *j*. The output breaks down coalescence events of
samples carrying different alleles at each site, using the
``{00,01,11}`` suffixes. Specifically:

-  The *i*-th row of the matrix in ``*.00.sumOverPairs.gz`` contains the
   sum of posterior probabilities for all pairs of samples that are
   homozygous ``0`` at site *i*.
-  The *i*-th row of the matrix in ``*.01.sumOverPairs.gz`` contains the
   sum of posterior probabilities for all pairs of samples that
   heterozygous at site *i*.
-  The *i*-th row of the matrix in ``*.11.sumOverPairs.gz`` contains the
   sum of posterior probabilities for all pairs of samples that are
   homozygous ``1`` at site *i*.

Tools, scripts, and analyses
----------------------------

These are some useful tools and scripts to be used with ASMC. They have
not yet been relocated to this repository, but are standalone tools that
are available at the indicated locations.

Tool to merge output of parallel ASMC jobs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The folder
```MERGE_POSTERIORS`` <https://github.com/PalamaraLab/ASMC_legacy/tree/master/TOOLS/MERGE_POSTERIORS>`__
contains the ``ASMCmergePosteriorSums.jar`` program, which may be used
to merge the output of different ASMC jobs. You can type

::

   java -jar TOOLS/MERGE_POSTERIORS/ASMCmergePosteriorSums.jar  -h

for a list of command line options. Also see the ``merge.sh`` file. This
tool assumes the decoding has been done using
``--majorMinorPosteriorSums``. You may normalize the output so that the
posterior sums to ``1`` for each site using the ``--norm`` flag. If you
used the ``--posteriorSums`` flag and you want to simply normalize the
output, you can run:

.. code:: bash

   zcat path/to/output/exampleFile.n100.array.merged.sumOverPairs.gz | \
       awk '{ c=0; for (i=1; i<=NF; i++) c+=$i; for (i=1; i<=NF; i++) $i/=c; print; }' | \
       gzip -c -v - > path/to/output/exampleFile.n100.array.merged.norm.sumOverPairs.gz

Script to visualize average coalescence density in a region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The folder
```PLOT_POSTERIORS`` <https://github.com/PalamaraLab/ASMC_legacy/tree/master/TOOLS/PLOT_POSTERIORS>`__
contains the ``plotPosteriorHeatMap.py`` program, which can be used to
visualize the coalescence posterior density in specific regions (e.g.
figures 3.b, 3.c, and S7 in the ASMC paper). For an example, see the
`data <https://palamaralab.github.io/software/asmc/data/>`__ page, where
the
`UKBB_posteriors <https://www.stats.ox.ac.uk/~palamara/ASMC/data/UKBB_posteriors.180813.tar>`__
files can be downloaded to generate the figures from the paper.

Density of recent coalescence (DRC) statistic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DRC statistic can be obtained by summing entries of the
``sumOverPairs`` matrices that correspond to the desired time interval,
and averaging across SNPs within a window. Please refer to the paper for
details.

A few notes on analyzing the DRC statistic: it is important to remember
that using a Gamma distribution as a null model for the DRC statistic is
an approximation (see paper), and further work is needed to derive an
exact null model. In particular the DRC cannot be larger than 1, so
approximate p-values obtained under a Gamma model for very large DRC
values will be artificially small. We found the Gamma approximation to
be reasonable under some conditions, such as a large recent effective
population size, which makes the DRC small in neutral regions, but this
approximation will not work well in all populations (e.g. isolated
populations) or for large time intervals. If you decide to use this
approach to detect candidate regions for selection, we recommend
carefully checking that the Gamma approximation is reasonable. Because
particularly large DRC values may lead to artificially small p-values,
you may also want to use the Gamma approximation only to determine an
approximate significance threshold, then report raw DRC values larger
than this threshold, rather than reporting approximate p-values.

Finally, we also suggest that you check that candidate regions detected
this way do not overlap problematic regions of the genome, such as
regions of very high/low recombination rate or LD, as well as large
structural variants, such as inversions, which may affect recombination.

Precomputed decoding quantities
-------------------------------

`This repository <https://github.com/PalamaraLab/ASMC_data>`__ contains
various data files related to ASMC that users might find helpful. `This
directory <https://github.com/PalamaraLab/ASMC_data/tree/main/decoding_quantities>`__
contains two sets of decoding quantities that have been precomputed
using the `PrepareDecoding
tool <https://github.com/PalamaraLab/PrepareDecoding>`__.

They are built using the European demographic model
`CEU.demo <https://github.com/PalamaraLab/ASMC_data/blob/main/demographies/CEU.demo>`__,
SNP allele frequencies from the UK Biobank in
`UKBB.frq <https://github.com/PalamaraLab/ASMC_data/blob/main/frequencies/UKBB.frq>`__,
and the time discretizations that can be found in `this
directory <https://github.com/PalamaraLab/ASMC_data/tree/main/discretizations>`__.
The discretizations contain several short time intervals in the recent
generations, followed by intervals calculated from the mutation age
distribution from generation ``2,000`` on. This enables getting more
fine-grained information for recent generation, though note that smaller
time intervals will contain fewer coalescent events on average.
