# ASMC Python API

FastSMC includes Python bindings which can be installed using pip:

```
pip install asmc-asmc
```

## Examples using the Python bindings

See the `notebooks` directory for an example:
- a [asmc example](notebooks/asmc.ipynb), where sensible defaults for parameters are chosen automatically

## API

The core Python API for FastSMC consists of the following classes:
- `ASMC`
- `DecodePairsReturnStruct`

### ASMC

The main `ASMC` object can be constructed minimally with an input file root, a decoding quantities file, and an output directory.
Simply construct a FastSMC object and call `run()` to generate output in the output file root:

```python
asmc = ASMC(in_dir=input_files_root, dq_file=dq_file, out_dir=output_files_root)
```

This creates an ASMC object with sensible defaults.
To fine-tune parameters you can instead create the ASMC object with an instance of `DecodingParams`.

## Finish this...
