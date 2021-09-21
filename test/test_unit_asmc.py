import numpy as np

import pathlib
import sys
import unittest

from asmc.asmc import *

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
if not data_dir.exists():
    print(f'ERROR. {data_dir} does not exist. Did you clone ASMC recursively with submodules?')
    sys.exit()


class TestASMC(unittest.TestCase):

    def setUp(self):
        inFileRoot = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        decodingQuantFile = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')
        self.sequenceLength = Data.countHapLines(inFileRoot)
        params = DecodingParams(inFileRoot, decodingQuantFile)
        self.data = Data(params)
        self.hmm = HMM(self.data, params)

    def test_initialization(self):
        self.assertGreater(len(self.data.individuals), 20)

    def test_sum_over_pairs_shape(self):
        ret = self.hmm.getDecodingReturnValues()
        self.assertEqual(ret.sumOverPairs.shape,
                         (self.sequenceLength, self.hmm.getDecodingQuantities().states))

    def test_decode_pair(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePair(0, 9)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 4)
        self.hmm.decodePair(1, 1)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 5)

    def test_decode_pairs(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePairs([0, 1], [9, 1])
        self.assertEqual(len(self.hmm.getBatchBuffer()), 5)

    def test_decode_pair_observation(self):
        self.assertEqual(len(self.hmm.getDecodingQuantities().discretization),
                         len(self.hmm.getDecodingQuantities().expectedTimes) + 1)
        self.assertEqual(self.data.sites, self.sequenceLength)

        for p in [
            self.hmm.makePairObs(1, 0, 2, 0),
            self.hmm.makePairObs(1, 0, 1, 0),
            self.hmm.makePairObs(2, 0, 2, 0)]:
            d = self.hmm.decode(p)
            self.assertEqual(len(d), len(self.hmm.getDecodingQuantities().expectedTimes))
            for i in range(len(d)):
                self.assertEqual(len(d[i]), self.data.sites)

    def test_finish_decoding(self):
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)
        self.hmm.decodePair(0, 9)
        self.assertEqual(len(self.hmm.getBatchBuffer()), 4)
        self.hmm.finishDecoding()
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)

    def test_fill_up_buffer(self):
        for i in range(1, (64 // 4) + 1):
            self.hmm.decodePair(0, i)
        # buffer should be empty now
        self.assertEqual(len(self.hmm.getBatchBuffer()), 0)


class TestASMCDecodingParams(unittest.TestCase):
    def test_no_compress(self):
        inFileRoot = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        decodingQuantFile = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')
        params = DecodingParams(inFileRoot, decodingQuantFile, compress=True,
                                skip_CSFS_distance=float('nan'))

        self.assertEqual(params.compress, True)
        self.assertEqual(params.skipCSFSdistance, float('inf'))

        data = Data(params)
        hmm = HMM(data, params)

        p = hmm.makePairObs(1, 0, 2, 0)
        d = hmm.decode(p)
        self.assertEqual(len(d), len(hmm.getDecodingQuantities().expectedTimes))
        for i in range(len(d)):
            self.assertEqual(len(d[i]), data.sites)


class TestASMCDecodePairsAPI(unittest.TestCase):
    def test_decode_pairs(self):

        input_files_root = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        dq_file = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')

        asmc = ASMC(input_files_root, dq_file)

        asmc.set_store_per_pair_posterior_mean(True)  # <-- true by default; others false by default
        asmc.set_store_per_pair_map(True)
        asmc.set_store_per_pair_posterior(True)
        asmc.set_store_sum_of_posterior(True)

        a = [1, 2, 3]
        b = [2, 3, 4]

        asmc.decode_pairs(a, b)

        return_vals = asmc.get_copy_of_results()

        for posterior in return_vals.per_pair_posteriors:
            self.assertTrue(np.allclose(np.sum(posterior, axis=0), 1.0, 1e-2))


if __name__ == "__main__":
    unittest.main()
