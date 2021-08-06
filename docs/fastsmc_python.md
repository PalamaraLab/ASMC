# FastSMC Python API

FastSMC includes Python bindings which can be installed using pip:

```
pip install asmc-asmc
```

## Examples using the Python bindings

See the `notebooks` directory for examples.
There are two Jupyter notebooks:
- a [minimal working example](notebooks/fastsmc-minimal.ipynb), where sensible defaults for parameters are chosen automatically
- a [more detailed example](notebooks/fastsmc.ipynb) that demonstrates how to customise parameters, how to convert the binary file to text format, and how to analyse the output if it is too large to fit in memory.

## API

The core Python API for FastSMC consists of the following classes:
- `FastSMC`
- `DecodingParams`
- `BinaryDataReader`

### FastSMC

The main `FastSMC` object can be constructed minimally with an input file root, a decoding quantities file, and an output directory.
Simply construct a FastSMC object and call `run()` to generate output in the output file root:

```python
fast_smc = FastSMC(in_dir=input_files_root, dq_file=dq_file, out_dir=output_files_root)
fast_smc.run()
```

This creates a FastSMC object with sensible defaults.
To fine-tune parameters you can instead create the FastSMC object with an instance of `DecodingParams`.

### DecodingParams

Create an empty `DecodingParams` object:

```python
params = DecodingParams()
```

The following parameters can be set:

```python
params.decodingQuantFile = dq_file
params.inFileRoot = input_files_root
params.outFileRoot = output_files_root
params.decodingModeString = 'array'
params.decodingMode = DecodingMode.arrayFolded
params.foldData = True
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
```

Finally, you can validate that the parameters are consistent for running FastSMC:

```python
assert params.validateParamsFastSMC()
```

Then, construct and run a `FastSMC` object using these parameters:

```python
fast_smc = FastSMC(params)
fast_smc.run()
```

### BinaryDataReader

If you turn on `BIN_OUT` in the decoding parameters, the `BinaryDataReader` class can read sequential lines in a file.
This is useful particularly if the output is too large to process entirely in memory.

```python
binary_data_reader = BinaryDataReader(output_files_root + '.1.1.FastSMC.bibd.gz')

while binary_data_reader.moreLinesInFile():
    line = binary_data_reader.getNextLine()
```

For each line, the following attributes and methods are available:

```python
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
```
