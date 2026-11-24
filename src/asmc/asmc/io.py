"""Construct ASMC from common phased-haplotype formats.

``ASMC.from_arrays`` takes genotypes and positions as arrays. These helpers
read them from the formats other tools in the workflow already produce:

* ``from_pgen``: a PLINK 2 ``.pgen`` / ``.pvar`` / ``.psam`` trio (needs
  ``Pgenlib``).
* ``from_vcf``: a phased VCF or BCF (needs ``cyvcf2``).
* ``from_haps``: an Eagle / SHAPEIT ``.hap[.gz]`` / ``.samples`` /
  ``.map[.gz]`` set, read by the file-based ``ASMC`` constructor.

Install the optional readers with ``pip install asmc-asmc[io]``.

The data must be phased: an unphased heterozygous call is an error.

``from_pgen`` and ``from_vcf`` take genetic positions from a PLINK-style map
(chromosome, SNP ID, cM, bp), interpolated linearly at each site, as the
file-based reader does. Sites beyond either end of the map take the cM of
that end. ASMC works in Morgans, so the cM values are divided by 100.

Every ``ASMC`` constructor, including these, turns on per-pair posterior
mean storage.
"""

from __future__ import annotations

import gzip
import importlib
import os
import pathlib
import warnings
from typing import Optional, Tuple

import numpy as np

from .asmc_python_bindings import ASMC, DecodingParams


__all__ = ["from_pgen", "from_vcf", "from_haps"]


PathLike = os.PathLike[str] | str


def _import_optional(module: str, package: str, caller: str):
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"{caller} needs the '{package}' package: "
            f"`pip install {package}` or `pip install asmc-asmc[io]`") from e


def _read_genetic_map(map_path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read a PLINK-style genetic map (chromosome, SNP ID, cM, bp).

    A first line that does not parse is taken as a header. Returns
    (bp, cM), sorted by bp.
    """
    phys = []
    cm = []
    opener = gzip.open if str(map_path).endswith(".gz") else open
    with opener(map_path, "rt") as f:
        for line_number, line in enumerate(f, start=1):
            parts = line.split()
            if not parts:
                continue
            try:
                if len(parts) < 4:
                    raise ValueError
                cm_value, bp = float(parts[2]), int(parts[3])
            except ValueError:
                if line_number == 1:
                    continue
                raise ValueError(
                    f"{map_path}, line {line_number}: expected chromosome, "
                    f"SNP ID, cM and bp, got {line.strip()!r}") from None
            phys.append(bp)
            cm.append(cm_value)
    if len(phys) < 2:
        raise ValueError(f"{map_path}: a genetic map needs at least two rows")
    phys = np.asarray(phys, dtype=np.int64)
    cm = np.asarray(cm, dtype=np.float64)
    order = np.argsort(phys, kind="stable")
    phys, cm = phys[order], cm[order]
    if np.any(np.diff(cm) < 0):
        raise ValueError(f"{map_path}: cM decreases as bp increases")
    return phys, cm


def _from_genotypes(geno: np.ndarray, phys: np.ndarray, map_file: PathLike,
                    dq_file: PathLike, use_ancestral: bool,
                    skip_csfs_distance: float, decoding_mode: str) -> ASMC:
    """Check positions, add genetic positions and call ASMC.from_arrays.

    geno is (sites, haplotypes) in {0, 1}; phys is in bp.
    """
    if phys.size == 0:
        raise ValueError("no sites to decode")
    if np.any(np.diff(phys) < 0):
        site = int(np.argmax(np.diff(phys) < 0)) + 1
        raise ValueError(
            "sites must be in increasing physical position, but "
            f"site {site} is at {phys[site]} bp after {phys[site - 1]} bp")
    if phys[0] < 0 or phys[-1] > np.iinfo(np.int32).max:
        raise ValueError("physical positions must fit in a 32-bit integer")

    map_phys, map_cm = _read_genetic_map(map_file)
    cm = np.interp(phys.astype(np.float64), map_phys.astype(np.float64),
                   map_cm)
    # Round to float32 before dividing, as the file reader does, so that both
    # paths give identical genetic positions.
    morgans = cm.astype(np.float32) / np.float32(100.0)

    return ASMC.from_arrays(
        genotypes=geno, physical_positions=phys.astype(np.int32),
        genetic_positions=morgans, dq_file=str(dq_file),
        use_ancestral=use_ancestral, skip_csfs_distance=skip_csfs_distance,
        decoding_mode=decoding_mode,
    )


def _read_pvar_positions(pvar: pathlib.Path) -> np.ndarray:
    """Read the POS column of a .pvar file.

    With no '#CHROM' header line, PLINK 2 reads the columns in .bim order
    (CHROM, ID, CM, POS, ...), so POS is the fourth column.
    """
    positions = []
    pos_column = None
    with open(pvar) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].split()
                if "POS" not in header:
                    raise ValueError(f"{pvar}: header line has no POS column")
                pos_column = header.index("POS")
                continue
            parts = line.split()
            if not parts:
                continue
            if pos_column is None:
                pos_column = 3
            try:
                positions.append(int(parts[pos_column]))
            except (IndexError, ValueError):
                raise ValueError(
                    f"{pvar}, line {line_number}: cannot read a position "
                    f"from {line.strip()!r}") from None
    return np.asarray(positions, dtype=np.int64)


def from_pgen(pgen: PathLike, map_file: PathLike, dq_file: PathLike,
              *, decoding_mode: str = "array",
              skip_csfs_distance: float = 0.0,
              use_ancestral: bool = False,
              variant_indices: Optional[np.ndarray] = None) -> ASMC:
    """Construct ASMC from a PLINK 2 (.pgen / .pvar / .psam) trio.

    pgen: path to the .pgen file; the .pvar is found beside it.
    map_file: PLINK-style genetic map (chromosome, SNP ID, cM, bp).
    dq_file: ASMC decoding quantities.
    variant_indices: optional increasing .pgen variant indices to keep
        (for example, to drop singletons).

    Every kept variant must be biallelic, with no missing calls and no
    unphased heterozygous calls.
    """
    pgenlib = _import_optional("pgenlib", "Pgenlib", "from_pgen")

    pgen = pathlib.Path(pgen)
    pvar = pgen.with_suffix(".pvar")
    if not pvar.exists():
        raise FileNotFoundError(f"missing .pvar beside {pgen}")
    pvar_phys = _read_pvar_positions(pvar)

    if variant_indices is None:
        variant_indices = np.arange(len(pvar_phys), dtype=np.int64)
    else:
        variant_indices = np.asarray(variant_indices, dtype=np.int64)
        if variant_indices.ndim != 1:
            raise ValueError("variant_indices must be one-dimensional")
        if variant_indices.size and (variant_indices[0] < 0 or
                                     variant_indices[-1] >= len(pvar_phys)):
            raise ValueError("variant_indices out of range for the .pvar")
        if np.any(np.diff(variant_indices) <= 0):
            raise ValueError("variant_indices must be strictly increasing")
    phys = pvar_phys[variant_indices]

    with pgenlib.PgenReader(bytes(pgen)) as reader:
        if reader.get_variant_ct() != len(pvar_phys):
            raise ValueError(
                f"{pgen} has {reader.get_variant_ct()} variants but {pvar} "
                f"has {len(pvar_phys)}")
        n_samples = reader.get_raw_sample_ct()
        alleles = np.empty((len(variant_indices), 2 * n_samples),
                           dtype=np.int32)
        phased = np.empty((len(variant_indices), n_samples), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_list(
            variant_indices.astype(np.uint32), alleles, phased)

    if alleles.size and alleles.min() < 0:
        raise ValueError(f"{pgen} has missing calls; ASMC needs hard calls")
    if alleles.size and alleles.max() > 1:
        raise ValueError(f"{pgen} has multiallelic variants; ASMC needs "
                         "biallelic variants")
    unphased = (alleles[:, 0::2] != alleles[:, 1::2]) & (phased == 0)
    if unphased.any():
        row = int(np.argmax(unphased.any(axis=1)))
        raise ValueError(
            f"{pgen} has an unphased heterozygous call at variant "
            f"{variant_indices[row]} ({phys[row]} bp); ASMC needs phased "
            "haplotypes")

    return _from_genotypes(alleles.astype(np.int8), phys, map_file, dq_file,
                           use_ancestral, skip_csfs_distance, decoding_mode)


def from_vcf(vcf: PathLike, map_file: PathLike, dq_file: PathLike,
             *, decoding_mode: str = "array",
             skip_csfs_distance: float = 0.0,
             use_ancestral: bool = False,
             region: Optional[str] = None) -> ASMC:
    """Construct ASMC from a phased VCF or BCF, read with cyvcf2.

    vcf: path to the VCF or BCF; region needs it to be indexed.
    map_file: PLINK-style genetic map (chromosome, SNP ID, cM, bp).
    dq_file: ASMC decoding quantities.
    region: optional region to read, such as "chr2" or "chr2:1-5000000".

    Sites that are not biallelic SNPs (indels, symbolic or multiple ALT
    alleles) and sites with a missing call are skipped, with a warning.
    The sites read must all be on one chromosome. Genotypes must be
    diploid, and heterozygous calls must be phased.
    """
    cyvcf2 = _import_optional("cyvcf2", "cyvcf2", "from_vcf")

    reader = cyvcf2.VCF(str(vcf))
    records = reader(region) if region else reader
    n_haps = 2 * len(reader.samples)
    chromosomes = set()
    positions = []
    rows = []
    skipped = 0
    for v in records:
        if (len(v.REF) != 1 or len(v.ALT) != 1 or len(v.ALT[0]) != 1
                or v.ALT[0] in ("*", ".")):
            skipped += 1
            continue
        gt = v.genotype.array()  # (samples, 3): two alleles, then phased
        if gt.shape[1] != 3:
            raise ValueError(f"{vcf}: {v.CHROM}:{v.POS} is not diploid")
        alleles = gt[:, :2]
        if (alleles < 0).any():
            skipped += 1
            continue
        if np.any((alleles[:, 0] != alleles[:, 1]) & (gt[:, 2] == 0)):
            raise ValueError(
                f"{vcf}: unphased heterozygous call at {v.CHROM}:{v.POS}; "
                "ASMC needs phased haplotypes")
        chromosomes.add(v.CHROM)
        positions.append(v.POS)
        rows.append(alleles.reshape(-1).astype(np.int8))
    reader.close()

    if len(chromosomes) > 1:
        raise ValueError(
            f"{vcf}: sites are on several chromosomes "
            f"({', '.join(sorted(chromosomes))}); choose one with region")
    if skipped:
        warnings.warn(f"from_vcf skipped {skipped} sites that are not "
                      "biallelic SNPs or have a missing call")

    geno = (np.stack(rows) if rows
            else np.empty((0, n_haps), dtype=np.int8))
    phys = np.asarray(positions, dtype=np.int64)
    return _from_genotypes(geno, phys, map_file, dq_file, use_ancestral,
                           skip_csfs_distance, decoding_mode)


def from_haps(haps_root: PathLike, dq_file: PathLike,
              *, decoding_mode: str = "array",
              skip_csfs_distance: float = 0.0,
              use_ancestral: bool = False) -> ASMC:
    """Construct ASMC from an Eagle / SHAPEIT .hap[.gz] / .samples / .map set.

    Uses the file-based reader with the same parameters as
    ``ASMC(haps_root, dq_file)``, adding the options given here.
    """
    haps_root = str(haps_root)
    params = DecodingParams(
        in_file_root=haps_root, dq_file=str(dq_file),
        out_file_root=haps_root, decoding_mode_string=decoding_mode,
        using_CSFS=False,  # turned on unless skip_csfs_distance is infinite
        use_ancestral=use_ancestral, skip_CSFS_distance=skip_csfs_distance,
        do_posterior_sums=True, do_major_minor_posterior_sums=True,
        do_per_pair_MAP=True)
    return ASMC(params)
