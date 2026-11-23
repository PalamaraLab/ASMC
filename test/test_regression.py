import pathlib
import sys
import tempfile
import unittest

import numpy as np

from asmc.asmc import *

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
if not data_dir.exists():
    print(f'ERROR. {data_dir} does not exist. Did you clone ASMC recursively with submodules?')
    sys.exit()


# These tests check that the Python bindings give the same results as the C++ regression tests, against the same
# snapshots. The C++ tests (test_regression.cpp, test_fastsmc_regression.cpp) cover many more configurations.
snapshot_dir = data_dir / 'testing' / 'regression'


class TestASMCRegression(unittest.TestCase):

    def setUp(self):
        in_file_root = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        decoding_quant_file = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')

        params = DecodingParams(in_file_root, decoding_quant_file, skip_CSFS_distance=float('nan'))

        self.asmc = ASMC(params)
        self.asmc.set_store_per_pair_posterior_mean(True)
        self.asmc.set_store_per_pair_map(True)

    def test_regression(self):
        # Same pairs as "ASMC regression: array data" in test_regression.cpp
        self.asmc.decode_pairs([0, 2, 17, 150], [1, 5, 230, 151])
        res = self.asmc.get_ref_of_results()

        self.assertEqual(res.per_pair_posterior_means.shape, (4, 6760))

        existing_post = np.loadtxt(str(snapshot_dir / 'asmc_array.perPairPosteriorMeans.gz'))
        self.assertTrue(np.allclose(res.per_pair_posterior_means, existing_post, rtol=0.001))

        existing_map = np.loadtxt(str(snapshot_dir / 'asmc_array.perPairMAPs.gz'))
        self.assertTrue(np.array_equal(res.per_pair_MAPs, existing_map))


class TestFastSMCRegression(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

        # Same configuration as "FastSMC regression: hashing" in test_fastsmc_regression.cpp
        self.params = DecodingParams()
        self.params.decodingQuantFile = str(data_dir / 'decoding_quantities' / '10-20-2000_CEU.decodingQuantities.gz')
        self.params.inFileRoot = str(data_dir / 'examples' / 'fastsmc' / 'example')
        self.params.outFileRoot = str(pathlib.Path(self.temp_dir.name) / 'out')
        self.params.decodingModeString = 'array'
        self.params.usingCSFS = True
        self.params.batchSize = 32
        self.params.recallThreshold = 3
        self.params.min_m = 1.5
        self.params.hashing = True
        self.params.FastSMC = True
        self.params.BIN_OUT = False
        self.params.outputIbdSegmentLength = True
        self.params.time = 50
        self.params.noConditionalAgeEstimates = True
        self.params.doPerPairMAP = True
        self.params.doPerPairPosteriorMean = True
        self.params.useKnownSeed = True

        assert self.params.validateParamsFastSMC()

        fast_smc = FastSMC(self.params)
        fast_smc.run()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_regression(self):
        original_text = np.loadtxt(str(snapshot_dir / 'fastsmc_hashing.ibd.gz'), usecols=(7, 8, 9, 10, 11, 12))
        generated_text = np.loadtxt(self.params.outFileRoot + ".1.1.FastSMC.ibd.gz", usecols=(7, 8, 9, 10, 11, 12))

        self.assertEqual(original_text.shape, generated_text.shape)
        self.assertTrue(np.allclose(original_text, generated_text, rtol=0.001))
