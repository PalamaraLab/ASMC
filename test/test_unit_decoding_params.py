import unittest
from asmc import DecodingParams, DecodingMode


class TestDecodingParams(unittest.TestCase):

    def setUp(self):
        self.inFileRoot = "FILES/EXAMPLE/exampleFile.n300.array"
        self.decodingQuantFile = "FILES/DECODING_QUANTITIES" \
            "/30-100-2000.decodingQuantities.gz"

    def test_array_folded(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile)
        self.assertEqual(params.decodingMode, DecodingMode.arrayFolded)
        self.assertEqual(params.compress, False)

    def test_sequence_folded(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile,
                                compress=True, skipCSFSdistance=float('nan'),
                                decodingModeString="sequence")
        self.assertEqual(params.decodingMode, DecodingMode.sequenceFolded)
        self.assertEqual(params.compress, True)

    def test_sequence(self):
        params = DecodingParams(self.inFileRoot, self.decodingQuantFile,
                                decodingModeString="sequence",
                                useAncestral=True)
        self.assertEqual(params.decodingMode, DecodingMode.sequence)


if __name__ == "__main__":
    unittest.main()
