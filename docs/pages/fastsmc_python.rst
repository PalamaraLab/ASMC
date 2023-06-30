FastSMC Python API
==================

-  `Examples using the Python
   bindings <#examples-using-the-python-bindings>`__
-  `API <#api>`__

   -  `FastSMC <#fastsmc>`__
   -  `DecodingParams <#decodingparams>`__
   -  `BinaryDataReader <#binarydatareader>`__

FastSMC includes Python bindings which can be installed using pip:

::

   pip install asmc-asmc

Before reading further you may wish to read the `FastSMC
docs <./fastsmc.md>`__. In particular, these sections are directly
relevant:

-  `Summary (TL;DR) <./fastsmc.md#input-file-formats>`__
-  `Input/output file formats <./fastsmc.md#output-format>`__
-  `Tools, scripts, and analyses <./fastsmc.md#binary-output>`__
-  `Precomputed decoding
   quantities <./fastsmc.md#relationship-to-asmc>`__

And, from the `ASMC docs <./asmc.md>`__:

-  `Precomputed decoding
   quantities <./asmc.md#precomputed-decoding-quantities>`__

Examples using the Python bindings
----------------------------------

See the ``notebooks`` directory for examples. There are two Jupyter
notebooks:

-  a `minimal working example <../notebooks/fastsmc-minimal.ipynb>`__,
   where sensible defaults for parameters are chosen automatically
-  a `more detailed example <../notebooks/fastsmc.ipynb>`__ that
   demonstrates how to customise parameters, how to convert the binary
   file to text format, and how to analyse the output if it is too large
   to fit in memory.

API
---

The core Python API for FastSMC consists of the following classes:

-  ``FastSMC``
-  ``DecodingParams``
-  ``BinaryDataReader``

FastSMC
~~~~~~~

The main ``FastSMC`` object can be constructed minimally with an input
file root, a decoding quantities file, and an output directory. Simply
construct a FastSMC object and call ``run()`` to generate output in the
output file root:

.. code:: python

   fast_smc = FastSMC(in_dir=input_files_root, dq_file=dq_file, out_dir=output_files_root)
   fast_smc.run()

This creates a FastSMC object with sensible defaults. To fine-tune
parameters you can instead create the FastSMC object with an instance of
``DecodingParams``.

DecodingParams
~~~~~~~~~~~~~~

Create an empty ``DecodingParams`` object:

.. code:: python

   params = DecodingParams()

The following parameters can be set:

.. code:: python

   params.decodingQuantFile = dq_file
   params.inFileRoot = input_files_root
   params.map_file = map_file  # Optional override for .map file, if not in input_files_root
   params.outFileRoot = output_files_root
   params.decodingModeString = 'array'
   params.usingCSFS = True
   params.batchSize = 32
   params.recallThreshold = 3
   params.min_m = 1.5
   params.hashing = True
   params.FastSMC = True
   params.BIN_OUT = True
   params.outputIbdSegmentLength = True
   params.time = 50
   params.noConditionalAgeEstimates = True
   params.doPerPairMAP = True
   params.doPerPairPosteriorMean = True
   params.hashingOnly = False

..

   Note: the ``hashingOnly`` flag has not been extensively tested. You
   may also want to look into `this
   repository <https://github.com/gusevlab/germline2>`__ for a
   standalone version.

Finally, you can validate that the parameters are consistent for running
FastSMC:

.. code:: python

   assert params.validateParamsFastSMC()

Then, construct and run a ``FastSMC`` object using these parameters:

.. code:: python

   fast_smc = FastSMC(params)
   fast_smc.run()

BinaryDataReader
~~~~~~~~~~~~~~~~

If you turn on ``BIN_OUT`` in the decoding parameters, the
``BinaryDataReader`` class can read sequential lines in a file. This is
useful particularly if the output is too large to process entirely in
memory.

.. code:: python

   binary_data_reader = BinaryDataReader(output_files_root + '.1.1.FastSMC.bibd.gz')

   while binary_data_reader.moreLinesInFile():
       line = binary_data_reader.getNextLine()

For each line, the following attributes and methods are available:

.. code:: python

   line.ind1FamId
   line.ind1Id
   line.ind1Hap
   line.ind2FamId
   line.ind2Id
   line.ind2Hap
   line.chromosome
   line.ibdStart
   line.ibdEnd
   line.lengthInCentimorgans
   line.ibdScore
   line.postEst
   line.mapEst

   line.toString()
