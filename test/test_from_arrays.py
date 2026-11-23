"""Tests for from_arrays, from_counts, load_genotypes, and decode_pairs_parallel."""

import numpy as np
import pathlib
import unittest

from asmc.asmc import ASMC, Data, DecodingParams

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
dq_file = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')
file_root = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')


def load_reference():
    """Load data from files and extract numpy arrays for from_arrays tests."""
    asmc = ASMC(file_root, dq_file)
    phys = np.array(asmc.get_physical_positions(), dtype=np.int32)
    gen = np.array(asmc.get_genetic_positions(), dtype=np.float32)
    n_sites = asmc.get_num_sites()
    n_haps = asmc.get_haploid_sample_size()

    # Read genotypes from the ASMC Data object via decode_pairs
    # We decode a few pairs from the file-based ASMC to get reference results
    params = DecodingParams(file_root, dq_file)
    data = Data(params)

    # Extract genotypes from Individual objects
    geno = np.zeros((n_sites, n_haps), dtype=np.int8)
    for ind_idx, ind in enumerate(data.individuals):
        for s in range(n_sites):
            geno[s, 2 * ind_idx] = int(ind.genotype1[s])
            geno[s, 2 * ind_idx + 1] = int(ind.genotype2[s])

    return asmc, phys, gen, geno, n_sites, n_haps


class TestFromArrays(unittest.TestCase):
    """Test ASMC.from_arrays construction and decoding."""

    @classmethod
    def setUpClass(cls):
        cls.ref_asmc, cls.phys, cls.gen, cls.geno, cls.n_sites, cls.n_haps = load_reference()

        # Decode reference pairs from file-based ASMC
        cls.pairs_a = [0, 2, 4]
        cls.pairs_b = [1, 3, 5]
        cls.ref_asmc.set_store_per_pair_posterior_mean(True)
        cls.ref_asmc.set_store_per_pair_map(True)
        cls.ref_asmc.decode_pairs(cls.pairs_a, cls.pairs_b)
        cls.ref_results = cls.ref_asmc.get_copy_of_results()

    def test_from_arrays_construction(self):
        """from_arrays should create a valid ASMC object with correct dimensions."""
        asmc = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        self.assertEqual(asmc.get_num_sites(), self.n_sites)
        self.assertEqual(asmc.get_haploid_sample_size(), self.n_haps)
        self.assertEqual(asmc.get_diploid_sample_size(), self.n_haps // 2)

    def test_from_arrays_matches_file_based(self):
        """from_arrays decode results should match file-based construction."""
        asmc = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        asmc.set_store_per_pair_posterior_mean(True)
        asmc.set_store_per_pair_map(True)
        asmc.decode_pairs(self.pairs_a, self.pairs_b)
        results = asmc.get_copy_of_results()

        np.testing.assert_allclose(
            np.array(results.per_pair_posterior_means),
            np.array(self.ref_results.per_pair_posterior_means),
            rtol=1e-4,
            err_msg="from_arrays posterior means differ from file-based"
        )

    def test_from_arrays_validation_odd_haps(self):
        """from_arrays should reject odd number of haplotypes."""
        bad_geno = np.zeros((10, 3), dtype=np.int8)
        phys = np.arange(10, dtype=np.int32)
        gen = np.linspace(0, 0.01, 10, dtype=np.float32)
        with self.assertRaises(Exception):
            ASMC.from_arrays(bad_geno, phys, gen, dq_file)

    def test_from_arrays_validation_position_length(self):
        """from_arrays should reject mismatched position array lengths."""
        phys_short = self.phys[:10]
        with self.assertRaises(Exception):
            ASMC.from_arrays(self.geno, phys_short, self.gen, dq_file)


class TestFromCounts(unittest.TestCase):
    """Test ASMC.from_counts + load_genotypes streaming."""

    @classmethod
    def setUpClass(cls):
        _, cls.phys, cls.gen, cls.geno, cls.n_sites, cls.n_haps = load_reference()

        # Compute allele counts from genotypes
        cls.dac = cls.geno.sum(axis=1).astype(np.int32)
        cls.tsc = np.full(cls.n_sites, cls.n_haps, dtype=np.int32)

    def test_from_counts_construction(self):
        """from_counts should create a valid ASMC with zero sample size initially."""
        asmc = ASMC.from_counts(self.dac, self.tsc, self.phys, self.gen, dq_file)
        self.assertEqual(asmc.get_num_sites(), self.n_sites)
        self.assertEqual(asmc.get_haploid_sample_size(), 0)

    def test_load_and_decode(self):
        """from_counts + load_genotypes should produce same results as from_arrays."""
        # Reference from from_arrays
        ref = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        ref.set_store_per_pair_posterior_mean(True)
        pairs_a, pairs_b = [0, 2], [1, 3]
        ref.decode_pairs(pairs_a, pairs_b)
        ref_means = np.array(ref.get_copy_of_results().per_pair_posterior_means)

        # Streaming via from_counts
        asmc = ASMC.from_counts(self.dac, self.tsc, self.phys, self.gen, dq_file)
        asmc.set_store_per_pair_posterior_mean(True)
        # Load subset: first 4 haplotypes (2 diploid samples)
        subset = self.geno[:, :4].copy()
        asmc.load_genotypes(subset)
        asmc.decode_pairs(pairs_a, pairs_b)
        stream_means = np.array(asmc.get_copy_of_results().per_pair_posterior_means)

        np.testing.assert_allclose(ref_means, stream_means, rtol=1e-4,
                                   err_msg="streaming decode differs from from_arrays")

    def test_load_clear_reload(self):
        """load, decode, clear, reload, decode again should work and give consistent results."""
        asmc = ASMC.from_counts(self.dac, self.tsc, self.phys, self.gen, dq_file)
        asmc.set_store_per_pair_posterior_mean(True)

        # First batch
        asmc.load_genotypes(self.geno[:, :4].copy())
        asmc.decode_pairs([0, 2], [1, 3])
        means1 = np.array(asmc.get_copy_of_results().per_pair_posterior_means).copy()

        # Clear and reload same data
        asmc.clear_genotypes()
        asmc.load_genotypes(self.geno[:, :4].copy())
        asmc.decode_pairs([0, 2], [1, 3])
        means2 = np.array(asmc.get_copy_of_results().per_pair_posterior_means)

        np.testing.assert_allclose(means1, means2, rtol=1e-5,
                                   err_msg="reload after clear gives different results")

    def test_load_genotypes_validation_wrong_sites(self):
        """load_genotypes should reject wrong number of sites."""
        asmc = ASMC.from_counts(self.dac, self.tsc, self.phys, self.gen, dq_file)
        bad = np.zeros((10, 4), dtype=np.int8)
        with self.assertRaises(Exception):
            asmc.load_genotypes(bad)

    def test_load_genotypes_validation_odd_haps(self):
        """load_genotypes should reject odd haplotype count."""
        asmc = ASMC.from_counts(self.dac, self.tsc, self.phys, self.gen, dq_file)
        bad = np.zeros((self.n_sites, 3), dtype=np.int8)
        with self.assertRaises(Exception):
            asmc.load_genotypes(bad)

    def test_from_counts_validation_mismatched_lengths(self):
        """from_counts should reject mismatched array lengths."""
        with self.assertRaises(Exception):
            ASMC.from_counts(self.dac[:10], self.tsc, self.phys, self.gen, dq_file)


class TestDecodePairsParallel(unittest.TestCase):
    """Test decode_pairs_parallel matches sequential decode_pairs."""

    @classmethod
    def setUpClass(cls):
        _, cls.phys, cls.gen, cls.geno, cls.n_sites, cls.n_haps = load_reference()

    def test_parallel_matches_sequential(self):
        """Parallel decoding should produce identical results to sequential."""
        pairs_a = [0, 2, 4, 6, 8]
        pairs_b = [1, 3, 5, 7, 9]

        # Sequential
        seq = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        seq.set_store_per_pair_posterior_mean(True)
        seq.set_store_per_pair_map(True)
        seq.decode_pairs(pairs_a, pairs_b)
        seq_results = seq.get_copy_of_results()

        # Parallel
        par = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        par.set_store_per_pair_posterior_mean(True)
        par.set_store_per_pair_map(True)
        par.decode_pairs_parallel(pairs_a, pairs_b, n_threads=2)
        par_results = par.get_copy_of_results()

        np.testing.assert_allclose(
            np.array(par_results.per_pair_posterior_means),
            np.array(seq_results.per_pair_posterior_means),
            rtol=1e-4,
            err_msg="parallel posterior means differ from sequential"
        )

    def test_parallel_single_thread_matches(self):
        """Parallel with n_threads=1 should fall back to sequential."""
        pairs_a, pairs_b = [0, 2], [1, 3]

        seq = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        seq.set_store_per_pair_posterior_mean(True)
        seq.decode_pairs(pairs_a, pairs_b)
        seq_means = np.array(seq.get_copy_of_results().per_pair_posterior_means)

        par = ASMC.from_arrays(self.geno, self.phys, self.gen, dq_file)
        par.set_store_per_pair_posterior_mean(True)
        par.decode_pairs_parallel(pairs_a, pairs_b, n_threads=1)
        par_means = np.array(par.get_copy_of_results().per_pair_posterior_means)

        np.testing.assert_allclose(seq_means, par_means, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
