"""Tests for asmc.asmc.io: from_pgen, from_vcf and from_haps."""

import gzip
import importlib.util
import pathlib
import tempfile
import unittest

import numpy as np

from asmc.asmc import ASMC, io

data_dir = pathlib.Path(__file__).resolve().parent.parent / 'ASMC_data'
dq_file = str(data_dir / 'decoding_quantities' / '30-100-2000_CEU.decodingQuantities.gz')
file_root = str(data_dir / 'examples' / 'asmc' / 'exampleFile.n300.array')
vcf_file = file_root + '.vcf.gz'
map_file = file_root + '.map.gz'

have_cyvcf2 = importlib.util.find_spec('cyvcf2') is not None
have_pgenlib = importlib.util.find_spec('pgenlib') is not None

pairs_a = [0, 2, 4]
pairs_b = [1, 3, 5]


def decode(asmc):
    asmc.decode_pairs(pairs_a, pairs_b)
    return np.array(asmc.get_copy_of_results().per_pair_posterior_means)


def assert_same_as_file_based(test, asmc):
    """asmc must hold exactly the sites, positions and decoding of the example."""
    ref = ASMC(file_root, dq_file)
    test.assertEqual(asmc.get_num_sites(), ref.get_num_sites())
    test.assertEqual(asmc.get_haploid_sample_size(), ref.get_haploid_sample_size())
    np.testing.assert_array_equal(np.array(asmc.get_physical_positions()), np.array(ref.get_physical_positions()))
    np.testing.assert_array_equal(np.array(asmc.get_genetic_positions()), np.array(ref.get_genetic_positions()))
    np.testing.assert_array_equal(decode(asmc), decode(ref))


def read_hap_file():
    """Positions and (sites, haplotypes) genotypes of the example's .hap.gz."""
    positions = []
    rows = []
    with gzip.open(file_root + '.hap.gz', 'rt') as f:
        for line in f:
            parts = line.split()
            positions.append(int(parts[2]))
            rows.append(np.array(parts[5:], dtype=np.int32))
    return np.array(positions), np.stack(rows)


def write_text(path, lines):
    with open(path, 'w') as f:
        f.write(''.join(line + '\n' for line in lines))


def write_small_map(path):
    write_text(path, ['1 a 0.0 1', '1 b 1.0 1000000'])


def write_small_vcf(path, sites, chrom='1'):
    """sites: (pos, ref, alt, calls), with calls such as '0|1' for each of four samples."""
    lines = ['##fileformat=VCFv4.2', '##contig=<ID=1>', '##contig=<ID=2>',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4']
    for site in sites:
        site_chrom = site[4] if len(site) > 4 else chrom
        pos, ref, alt, calls = site[:4]
        lines.append('\t'.join([site_chrom, str(pos), '.', ref, alt, '.', '.', '.', 'GT'] + calls))
    write_text(path, lines)


SMALL_CALLS = ['0|1', '1|0', '0|0', '1|1']


class TestGeneticMap(unittest.TestCase):

    def test_header_line_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / 'g.map'
            write_text(path, ['chr id cM bp', '1 a 0.5 10', '1 b 1.5 20'])
            phys, cm = io._read_genetic_map(path)
        np.testing.assert_array_equal(phys, [10, 20])
        np.testing.assert_array_equal(cm, [0.5, 1.5])

    def test_rejects_short_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / 'g.map'
            write_text(path, ['1 a 0.5 10', '1 b 20', '1 c 1.5 30'])
            with self.assertRaisesRegex(ValueError, 'line 2'):
                io._read_genetic_map(path)


@unittest.skipUnless(have_cyvcf2, 'needs cyvcf2')
class TestFromVcf(unittest.TestCase):

    def test_matches_file_based(self):
        assert_same_as_file_based(self, io.from_vcf(vcf_file, map_file, dq_file))

    def small(self, sites):
        """from_vcf on a small VCF, without the CSFS (it needs at least 300 haplotypes)."""
        with tempfile.TemporaryDirectory() as tmp:
            vcf = pathlib.Path(tmp) / 's.vcf'
            small_map = pathlib.Path(tmp) / 's.map'
            write_small_vcf(vcf, sites)
            write_small_map(small_map)
            return io.from_vcf(vcf, small_map, dq_file, skip_csfs_distance=float('inf'))

    def test_skips_indels_and_warns(self):
        sites = [(100, 'A', 'G', SMALL_CALLS), (200, 'A', 'AT', SMALL_CALLS),
                 (300, 'C', 'T', SMALL_CALLS), (400, 'G', '<DEL>', SMALL_CALLS)]
        with self.assertWarnsRegex(UserWarning, 'skipped 2 sites'):
            asmc = self.small(sites)
        np.testing.assert_array_equal(np.array(asmc.get_physical_positions()), [100, 300])

    def test_rejects_unphased_heterozygous_call(self):
        sites = [(100, 'A', 'G', SMALL_CALLS), (200, 'C', 'T', ['0|1', '1/0', '0|0', '1|1'])]
        with self.assertRaisesRegex(ValueError, 'unphased heterozygous call at 1:200'):
            self.small(sites)

    def test_accepts_unphased_homozygous_call(self):
        sites = [(100, 'A', 'G', SMALL_CALLS), (200, 'C', 'T', ['0|1', '1|0', '0/0', '1/1'])]
        self.assertEqual(self.small(sites).get_num_sites(), 2)

    def test_rejects_several_chromosomes(self):
        sites = [(100, 'A', 'G', SMALL_CALLS), (200, 'C', 'T', SMALL_CALLS, '2')]
        with self.assertRaisesRegex(ValueError, 'several chromosomes'):
            self.small(sites)

    def test_rejects_decreasing_positions(self):
        sites = [(300, 'A', 'G', SMALL_CALLS), (200, 'C', 'T', SMALL_CALLS)]
        with self.assertRaisesRegex(ValueError, 'increasing physical position'):
            self.small(sites)


@unittest.skipUnless(have_pgenlib, 'needs Pgenlib')
class TestFromPgen(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import pgenlib
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = pathlib.Path(cls.tmp.name)
        cls.positions, alleles = read_hap_file()
        n_samples = alleles.shape[1] // 2

        cls.pgen = cls.dir / 'phased.pgen'
        with pgenlib.PgenWriter(bytes(cls.pgen), n_samples, len(cls.positions), False,
                                hardcall_phase_present=True) as writer:
            writer.append_alleles_batch(np.ascontiguousarray(alleles), all_phased=True)
        cls.pvar_lines = ['#CHROM\tPOS\tID\tREF\tALT'] + [f'1\t{p}\tsnp{i}\tA\tG' for i, p in enumerate(cls.positions)]
        write_text(cls.dir / 'phased.pvar', cls.pvar_lines)

        cls.unphased_pgen = cls.dir / 'unphased.pgen'
        with pgenlib.PgenWriter(bytes(cls.unphased_pgen), n_samples, len(cls.positions), False) as writer:
            writer.append_alleles_batch(np.ascontiguousarray(alleles), all_phased=False)
        write_text(cls.dir / 'unphased.pvar', cls.pvar_lines)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def copy_pgen(self, name):
        """A copy of the phased .pgen under another name, for a different .pvar."""
        path = self.dir / f'{name}.pgen'
        path.write_bytes(self.pgen.read_bytes())
        return path

    def test_matches_file_based(self):
        assert_same_as_file_based(self, io.from_pgen(self.pgen, map_file, dq_file))

    def test_reads_headerless_pvar(self):
        pgen = self.copy_pgen('headerless')
        write_text(pgen.with_suffix('.pvar'), [f'1\tsnp{i}\t0\t{p}\tG\tA' for i, p in enumerate(self.positions)])
        assert_same_as_file_based(self, io.from_pgen(pgen, map_file, dq_file))

    def test_rejects_unparseable_pvar_row(self):
        pgen = self.copy_pgen('bad_row')
        lines = list(self.pvar_lines)
        lines[5] = '1\tnot_a_position\tsnp4\tA\tG'
        write_text(pgen.with_suffix('.pvar'), lines)
        with self.assertRaisesRegex(ValueError, 'line 6'):
            io.from_pgen(pgen, map_file, dq_file)

    def test_rejects_pvar_with_wrong_variant_count(self):
        pgen = self.copy_pgen('short')
        write_text(pgen.with_suffix('.pvar'), self.pvar_lines[:-1])
        with self.assertRaisesRegex(ValueError, 'variants but'):
            io.from_pgen(pgen, map_file, dq_file)

    def test_rejects_unphased_heterozygous_call(self):
        with self.assertRaisesRegex(ValueError, 'unphased heterozygous call'):
            io.from_pgen(self.unphased_pgen, map_file, dq_file)

    def test_variant_indices(self):
        keep = np.arange(0, len(self.positions), 2)
        asmc = io.from_pgen(self.pgen, map_file, dq_file, variant_indices=keep)
        np.testing.assert_array_equal(np.array(asmc.get_physical_positions()), self.positions[keep])

    def test_rejects_decreasing_variant_indices(self):
        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            io.from_pgen(self.pgen, map_file, dq_file, variant_indices=[3, 1, 2])


class TestFromHaps(unittest.TestCase):

    def test_matches_file_based(self):
        assert_same_as_file_based(self, io.from_haps(file_root, dq_file))

    def test_passes_options(self):
        asmc = io.from_haps(file_root, dq_file, use_ancestral=True, skip_csfs_distance=float('inf'))
        params = asmc.get_decoding_params()
        self.assertTrue(params.useAncestral)
        self.assertFalse(params.foldData)
        self.assertEqual(params.skipCSFSdistance, float('inf'))
        self.assertFalse(params.usingCSFS)


if __name__ == '__main__':
    unittest.main()
