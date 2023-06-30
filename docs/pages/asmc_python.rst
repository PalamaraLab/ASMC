ASMC Python API
===============

-  `Examples using the Python
   bindings <#examples-using-the-python-bindings>`__
-  `API <#api>`__

   -  `ASMC <#asmc>`__
   -  `DecodePairsReturnStruct <#decodepairsreturnstruct>`__

ASMC includes Python bindings which can be installed using pip:

::

   pip install asmc-asmc

Before reading further you may wish to read the `ASMC
docs <./asmc.md>`__. In particular, these sections are directly
relevant:

-  `Summary (TL;DR) <./asmc.md#summary-tldr>`__
-  `Input/output file formats <./asmc.md#inputoutput-file-formats>`__
-  `Tools, scripts, and
   analyses <./asmc.md#tools-scripts-and-analyses>`__
-  `Precomputed decoding
   quantities <./asmc.md#precomputed-decoding-quantities>`__

Examples using the Python bindings
----------------------------------

See the ``notebooks`` directory for examples. There are two Jupyter
notebooks:

-  a `minimal working example <../notebooks/asmc-minimal.ipynb>`__,
   where sensible defaults for parameters are chosen automatically
-  a `more detailed example <../notebooks/asmc.ipynb>`__ that
   demonstrates how to customise parameters

API
---

The core Python API for ASMC consists of the following classes:

-  ``ASMC``
-  ``DecodingParams``
-  ``DecodePairsReturnStruct``

ASMC
~~~~

The main ``ASMC`` object can be constructed minimally with an input file
root and a decoding quantities file. Optional parameters are the output
file root and the decoding mode. The full signature (with defaults
indicated) is as follows:

.. code:: python

   asmc = ASMC(
       in_dir=input_files_root,          # path to the 
       dq_file=dq_file,                  # path to the decoding quantities file
       out_dir="",                       # location of output files (default is the input file root)
       decoding_mode="array"             # one of "array" or "sequence"
   )

This creates an ASMC object with sensible defaults. To fine-tune
parameters you can instead create the ASMC object with an instance of
``DecodingParams``:

.. code:: python

   # These are the arguments (with defaults indicated) for constructing a decoding paramters object
   params = DecodingParams(
       in_file_root=input_files_root,
       dq_file=dq_file,
       map_file="",                          # Optional override for map|map.gz file, if not in in_file_root
       out_file_root="",
       jobs=1,                               # Number of jobs being done in total
       job_ind=1,                            # Job index (0, ..., jobs)
       decoding_mode_string="array",         # One of {"squence", "array"}
       decoding_sequence=False,
       using_CSFS=True,                      # Whether to use CSFS
       compress=False,                       # Compress emission to binary (no CSFS)
       use_ancestral=False,                  # Assume ancestral alleles are coded as 1 in input (will assume 1 = minor otherwise)
       skip_CSFS_distance=0.0,               # Genetic distance between two CSFS emissions
       no_batches=False,                     # Decode with no vectorization (do not use without good reason)
       do_posterior_sums=False,
       do_per_pair_posterior_mean=False,
       expected_coal_times_file="",
       within_only=False,
       do_major_minor_posterior_sums=False,  # 
       do_per_pair_MAP=False    
   )

   asmc = ASMC(params)

You can specify the outputs that will be calculated with the following
methods (with defaults indicated):

.. code:: python

   # Per pair posterior mean, MAP and full posteriors, as well as the sum of posteriors can be stored in matrices
   asmc.set_store_per_pair_posterior_mean(True)  # <-- true by default; others false by default
   asmc.set_store_per_pair_map(False)
   asmc.set_store_per_pair_posterior(False)
   asmc.set_store_sum_of_posterior(False)

   # Per pair posterior mean and MAP can be written to file. This is typically slow.
   asmc.set_write_per_pair_posterior_mean(False)
   asmc.set_write_per_pair_map(False)

Finally, the ASMC method ``decode_pairs`` will run the analysis. There
are three different signatures available:

.. code:: python

   a = [1, 2, 3]
   b = [4, 5, 6]

   a_str = [f"1_{x}_1" for x in range(1,149)]
   b_str = [f"1_{x}_2" for x in range(1,149)]

   asmc.decode_pairs(a, b)           # two lists of haplotype indices
   asmc.decode_pairs(a_str, b_str)   # two lists of haplotype IDs, with _1 and _2 indicating the haplotype
   asmc.decode_pairs()               # <-- decode all pairs in the dataset

The results can then be accessed either by copy or reference:

.. code:: python

   return_vals = asmc.get_copy_of_results()
   return_vals_ref = asmc.get_ref_of_results()

Getting the values by reference is safe if you are only planning to call
``decode_pairs`` once, or if you are performing calculations that do not
require the results to persist after the first call to ``decode_pairs``.
If you call ``decode_pairs`` multiple times, the results will be
overwritten, so you should ensure you get results by copy.

DecodePairsReturnStruct
~~~~~~~~~~~~~~~~~~~~~~~

The return structure will contain results based on the options selected
on the ASMC object before calling ``decode_pairs``.

.. code:: python

   # The index information for the pairs decoded
   return_vals.per_pair_indices

   # The `per_pair_posteriors` option gives the largest amount of information: a list of 2D numpy arrays
   # The list has length numPairs, and each 2D array has size (numStates x numSites)
   return_vals.per_pair_posteriors

   # The sum of posteriors is a single 2D numpy array of size (numStates x numSites)
   return_vals.sum_of_posteriors

   # Turning on the per_pair_posteriors flag gives you the the following:
   # A 2D numpy array with posterior means, of size (numPairs x numSites)
   return_vals.per_pair_posterior_means
   # Two 1D numpy arrays with the column-wise min and argmin of this array:
   return_vals.min_posterior_means
   return_vals.argmin_posterior_means

   # Turning on the per_pair_MAPs flag gives you the the following:
   # A 2D numpy array with posterior MAPs, of size (numPairs x numSites)
   return_vals.per_pair_MAPs
   # Two 1D numpy arrays with the column-wise min and argmin of this array:
   return_vals.min_MAPs
   return_vals.argmin_MAPs

Finally, the ASMC object can also return the list of expected coalescent
times from the decoding quantities file: asmc.get_expected_times()
