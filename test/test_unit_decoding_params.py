import pathlib
import sys
import unittest

from asmc.asmc import DecodingParams, DecodingMode

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
if not data_dir.exists():
    print(f'ERROR. {data_dir} does not exist. Did you clone ASMC recursively with submodules?')
    sys.exit()


class TestDecodingParams(unittest.TestCase):

    def setUp(self):
        self.inFileRoot = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
        self.decodingQuantFile = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')

    def test_array_folded(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile)
        self.assertEqual(params.decodingMode, DecodingMode.arrayFolded)
        self.assertEqual(params.compress, False)

    def test_sequence_folded(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile,
                                compress=True, skip_CSFS_distance=float('nan'),
                                decoding_mode_string="sequence")
        self.assertEqual(params.decodingMode, DecodingMode.sequenceFolded)
        self.assertEqual(params.compress, True)

    def test_sequence(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile,
                                decoding_mode_string="sequence",
                                use_ancestral=True)
        self.assertEqual(params.decodingMode, DecodingMode.sequence)


if __name__ == "__main__":
    unittest.main()
