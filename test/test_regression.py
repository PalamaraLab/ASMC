import pathlib
import sys
import unittest

import numpy as np

from asmc.asmc import *

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
if not data_dir.exists():
    print(f'ERROR. {data_dir} does not exist. Did you clone ASMC recursively with submodules?')
    sys.exit()


class TestASMCRegression(unittest.TestCase):

    def setUp(self):
        in_file_root = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        decoding_quant_file = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')

        params = DecodingParams(in_file_root, decoding_quant_file, do_posterior_sums=True)

        self.asmc = ASMC(params)
        self.asmc.set_store_per_pair_posterior_mean(True)
        self.asmc.set_store_per_pair_map(True)

    def test_regression(self):
        self.asmc.decode_pairs([1, 2, 3], [2, 3, 4])
        res = self.asmc.get_ref_of_results()

        self.assertEqual(res.per_pair_posterior_means.shape[0], 3)
        self.assertEqual(res.per_pair_posterior_means.shape[1], 6760)

        existing_post = np.loadtxt(
            str(data_dir / 'testing' / 'asmc' / 'regression' / 'regression.perPairPosteriorMeans.gz'))
        self.assertEqual(np.allclose(res.per_pair_posterior_means, existing_post), True)

        existing_map = np.loadtxt(
            str(data_dir / 'testing' / 'asmc' / 'regression' / 'regression.perPairMAP.gz'))
        self.assertEqual(np.allclose(res.per_pair_MAPs, existing_map), True)


class TestFastSMCRegression(unittest.TestCase):

    def setUp(self):
        # Create decoding params object with required options
        self.params = DecodingParams()
        self.params.decodingQuantFile = str(data_dir / 'decoding_quantities' / '10-20-2000_CEU.decodingQuantities.gz')
        self.params.inFileRoot = str(data_dir / 'examples' / 'fastsmc' / 'example')
        self.params.outFileRoot = '/tmp/FastSMCresults'
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

    def test_regression(self):
        original_text = np.loadtxt(str(data_dir / 'testing' / 'fastsmc' / 'regression' / 'regression_output.ibd.gz'),
                                   usecols=(7, 8, 9, 10, 11))
        generated_text = np.loadtxt(self.params.outFileRoot + ".1.1.FastSMC.ibd.gz", usecols=(7, 8, 9, 10, 11))

        self.assertEqual(original_text.shape, generated_text.shape)
        self.assertEqual(np.allclose(original_text, generated_text), True)


class TestFastSMCRegressionWithoutHashing(unittest.TestCase):

    def setUp(self):
        # Create decoding params object with required options
        self.params = DecodingParams()
        self.params.decodingQuantFile = str(data_dir / 'decoding_quantities' / '10-20-2000_CEU.decodingQuantities.gz')
        self.params.inFileRoot = str(data_dir / 'examples' / 'fastsmc' / 'example')
        self.params.outFileRoot = '/tmp/FastSMCresults'
        self.params.decodingModeString = 'array'
        self.params.usingCSFS = True
        self.params.batchSize = 32
        self.params.recallThreshold = 3
        self.params.min_m = 1.5
        self.params.hashing = False
        self.params.FastSMC = True
        self.params.BIN_OUT = False
        self.params.outputIbdSegmentLength = True
        self.params.time = 50
        self.params.noConditionalAgeEstimates = True
        self.params.doPerPairMAP = True
        self.params.doPerPairPosteriorMean = True
        self.params.jobInd = 7
        self.params.jobs = 25
        self.params.useKnownSeed = True

        assert self.params.validateParamsFastSMC()

        fast_smc = FastSMC(self.params)
        fast_smc.run()

    def test_regression(self):
        original_text = np.loadtxt(
            str(data_dir / 'testing' / 'fastsmc' / 'regression' / 'regression_output_no_hashing.ibd.gz'),
            usecols=(7, 8, 9, 10, 11))
        generated_text = np.loadtxt(self.params.outFileRoot + ".7.25.FastSMC.ibd.gz", usecols=(7, 8, 9, 10, 11))

        self.assertEqual(original_text.shape, generated_text.shape)
        self.assertEqual(np.allclose(original_text, generated_text), True)
