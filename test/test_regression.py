import os
import unittest
import numpy as np
from pathlib import Path
from asmc import (
    HMM,
    DecodingQuantities,
    DecodingParams,
    Data,
    DecodingMode,
    FastSMC,
)


class TestASMCRegression(unittest.TestCase):

    def setUp(self):
        inFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        decodingQuantFile = "FILES/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz"
        self.params = DecodingParams(inFileRoot, decodingQuantFile, doPosteriorSums=True)
        self.data = Data(self.params)
        self.hmm = HMM(self.data, self.params)

    def test_regression(self):
        oldSumOverPairs = np.loadtxt(Path(__file__).parent / 'data' /
                                      'regression_test_original.gz')
        self.hmm.decodeAll(self.params.jobs, self.params.jobInd)
        ret = self.hmm.getDecodingReturnValues()
        self.assertEqual(np.allclose(ret.sumOverPairs, oldSumOverPairs), True)


class TestFastSMCRegression(unittest.TestCase):

    def setUp(self):
        self.file_dir = os.path.join(os.getcwd(), 'FILES', 'FASTSMC_EXAMPLE')
        self.name_prefix = 'example'

        # Create decoding params object with required options
        self.params = DecodingParams()
        self.params.decodingQuantFile = os.path.join(self.file_dir, '{}.decodingQuantities.gz'.format(self.name_prefix))
        self.params.inFileRoot = os.path.join(self.file_dir, self.name_prefix)
        self.params.outFileRoot = os.path.join('/tmp/FastSMCresults')
        self.params.decodingModeString = 'array'
        self.params.foldData = True
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

        original_text = np.loadtxt(os.path.join(self.file_dir, 'regression_output.ibd.gz'), usecols=(7, 8, 9, 10, 11))
        generated_text = np.loadtxt(self.params.outFileRoot + ".1.1.FastSMC.ibd.gz", usecols=(7, 8, 9, 10, 11))

        self.assertEqual(original_text.shape, generated_text.shape)
        self.assertEqual(np.allclose(original_text, generated_text), True)


class TestFastSMCRegressionWithoutHashing(unittest.TestCase):

    def setUp(self):
        self.file_dir = os.path.join(os.getcwd(), 'FILES', 'FASTSMC_EXAMPLE')
        self.name_prefix = 'example'

        # Create decoding params object with required options
        self.params = DecodingParams()
        self.params.decodingQuantFile = os.path.join(self.file_dir, '{}.decodingQuantities.gz'.format(self.name_prefix))
        self.params.inFileRoot = os.path.join(self.file_dir, self.name_prefix)
        self.params.outFileRoot = os.path.join('/tmp/FastSMCresults')
        self.params.decodingModeString = 'array'
        self.params.foldData = True
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
        self.params.jobs = 9
        self.params.useKnownSeed = True

        assert self.params.validateParamsFastSMC()

        fast_smc = FastSMC(self.params)
        fast_smc.run()

    def test_regression(self):
        original_text = np.loadtxt(os.path.join(self.file_dir, 'regression_output_no_hashing.ibd.gz'),
                                   usecols=(7, 8, 9, 10, 11))
        generated_text = np.loadtxt(self.params.outFileRoot + ".7.9.FastSMC.ibd.gz", usecols=(7, 8, 9, 10, 11))

        self.assertEqual(original_text.shape, generated_text.shape)
        self.assertEqual(np.allclose(original_text, generated_text), True)
