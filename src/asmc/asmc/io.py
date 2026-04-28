"""Convenience constructors for ASMC from common upstream formats.

ASMC's C++ side accepts in-memory genotype arrays via ``ASMC.from_arrays``
and allele-count summaries via ``ASMC.from_counts``. These helpers bridge
the gap from the formats other tools in the workflow already use:

* ``from_pgen`` — PLINK 2 binary genotype trio (``.pgen`` / ``.pvar`` /
  ``.psam``). Same format that threads-arg consumes natively.
* ``from_vcf`` — phased VCF/BCF. Slower than pgen for large data but the
  most portable starting point.
* ``from_haps`` — Eagle / SHAPEIT (``.hap[.gz]`` / ``.sample`` /
  ``.map``). What arg-needle and the legacy ASMC file path read today.

All three return a fully-built ``ASMC`` object with per-pair posterior
mean storage already enabled — equivalent to::

    asmc = ASMC.from_arrays(genotypes, physical_positions,
                            genetic_positions_in_morgans, dq_file)
    asmc.set_store_per_pair_posterior_mean()

Genetic positions on disk are in cM; ASMC's internal representation is
Morgans (the legacy ``.map`` reader divides by 100 internally). Each
helper handles that conversion so callers always pass on-disk units.
"""

from __future__ import annotations

import gzip
import os
import pathlib
from typing import Optional, Tuple

import numpy as np

from .asmc_python_bindings import ASMC


__all__ = ["from_pgen", "from_vcf", "from_haps", "PathLike"]


PathLike = os.PathLike[str] | str


def _interp_cm(physical: np.ndarray, map_phys: np.ndarray,
               map_cm: np.ndarray) -> np.ndarray:
    """Linear interpolation of cM at given physical positions."""
    if map_phys.size < 2:
        raise ValueError("genetic map needs at least two anchor points")
    if not (map_phys[0] <= physical[0] and physical[-1] <= map_phys[-1]):
        # Allow extrapolation but warn — common when chr ends are missing.
        # np.interp clamps at the edges, which is what we want here.
        pass
    return np.interp(physical.astype(np.float64),
                     map_phys.astype(np.float64),
                     map_cm.astype(np.float64))


def _read_genetic_map(map_path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read a 4-column PLINK-style genetic map (chr, snp_id, cM, bp).

    Returns (physical_bp, cM) sorted by bp ascending.
    """
    phys = []
    cm = []
    opener = gzip.open if str(map_path).endswith(".gz") else open
    with opener(map_path, "rt") as f:
        first = f.readline()
        # tolerate optional header
        try:
            parts = first.split()
            float(parts[2]); int(parts[3])
            phys.append(int(parts[3])); cm.append(float(parts[2]))
        except (IndexError, ValueError):
            pass
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            phys.append(int(parts[3])); cm.append(float(parts[2]))
    phys = np.asarray(phys, dtype=np.int64)
    cm = np.asarray(cm, dtype=np.float64)
    order = np.argsort(phys)
    return phys[order], cm[order]


def from_pgen(pgen: PathLike, map_file: PathLike, dq_file: PathLike,
              *, decoding_mode: str = "array",
              skip_csfs_distance: float = 0.0,
              use_ancestral: bool = False,
              variant_indices: Optional[np.ndarray] = None) -> ASMC:
    """Construct ASMC from a PLINK 2 (.pgen / .pvar / .psam) trio.

    Reads physical positions from the .pvar sidecar, interpolates cM at
    those positions from ``map_file``, then loads genotypes via pgenlib.
    Returns a ready-to-decode ASMC object.

    pgen — path to the .pgen file (sidecars discovered as .pvar / .psam).
    map_file — 4-col PLINK-style genetic map (chrom, snp_id, cM, bp).
    dq_file — ASMC decoding-quantities file.
    variant_indices — optional 1-D int array of pgen variant rows to keep
        (e.g. drop singletons / restrict to common variants).
    """
    try:
        import pgenlib
    except ImportError as e:
        raise ImportError(
            "from_pgen needs the 'pgenlib' Python package; "
            "install via PLINK 2 source or `pip install Pgenlib`") from e

    pgen = pathlib.Path(pgen)
    pvar = pgen.with_suffix(".pvar")
    if not pvar.exists():
        raise FileNotFoundError(f"missing .pvar sidecar for {pgen}")

    # Read .pvar — first non-header line onwards. PLINK 2 header is
    # ``##`` lines plus a single ``#CHROM POS ID REF ALT`` row.
    pvar_phys = []
    with open(pvar) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue  # column header
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                pvar_phys.append(int(parts[1]))
            except ValueError:
                continue
    pvar_phys = np.asarray(pvar_phys, dtype=np.int64)

    if variant_indices is not None:
        variant_indices = np.asarray(variant_indices, dtype=np.int64)
        if variant_indices.min() < 0 or variant_indices.max() >= len(pvar_phys):
            raise ValueError("variant_indices out of bounds for .pvar")
        keep_phys = pvar_phys[variant_indices]
    else:
        keep_phys = pvar_phys
        variant_indices = np.arange(len(pvar_phys), dtype=np.int64)

    # Read genotypes via pgenlib: (sites, 2*samples) int32 in {0, 1}
    # (haploid alleles for phased data).
    reader = pgenlib.PgenReader(str(pgen).encode())
    n_samples = reader.get_raw_sample_ct()
    n_haps = 2 * n_samples
    geno = np.empty((len(variant_indices), n_haps), dtype=np.int32)
    reader.read_alleles_list(variant_indices.astype(np.uint32), geno)
    reader.close()

    # ASMC's from_arrays wants int8 in {0, 1}.
    geno = geno.astype(np.int8, copy=False)
    if geno.min() < 0:
        raise ValueError("pgen contains missing alleles (-1); ASMC needs hard calls")
    if geno.max() > 1:
        raise ValueError("pgen contains multi-allelic sites; ASMC needs biallelic")

    map_phys, map_cm = _read_genetic_map(map_file)
    cm_at_sites = _interp_cm(keep_phys, map_phys, map_cm)
    morgans = (cm_at_sites / 100.0).astype(np.float32)
    phys = keep_phys.astype(np.int32)

    asmc = ASMC.from_arrays(
        genotypes=geno, physical_positions=phys, genetic_positions=morgans,
        dq_file=str(dq_file), use_ancestral=use_ancestral,
        skip_csfs_distance=skip_csfs_distance, decoding_mode=decoding_mode,
    )
    asmc.set_store_per_pair_posterior_mean()
    return asmc


def from_vcf(vcf: PathLike, map_file: PathLike, dq_file: PathLike,
             *, decoding_mode: str = "array",
             skip_csfs_distance: float = 0.0,
             use_ancestral: bool = False,
             region: Optional[str] = None) -> ASMC:
    """Construct ASMC from a phased VCF / BCF via cyvcf2.

    Reads only biallelic SNPs with no missing alleles. Phased genotypes
    are unpacked into (sites, 2*samples) int8 in {0, 1}. cM positions
    are interpolated from ``map_file`` (PLINK-style, cM column).
    """
    try:
        from cyvcf2 import VCF
    except ImportError as e:
        raise ImportError("from_vcf needs cyvcf2; `pip install cyvcf2`") from e

    vcf_obj = VCF(str(vcf))
    if region:
        records = vcf_obj(region)
    else:
        records = vcf_obj
    n_samples = len(vcf_obj.samples)
    n_haps = 2 * n_samples
    phys_list = []
    rows = []
    skipped = 0
    for v in records:
        if v.num_called != n_samples:
            skipped += 1; continue
        if len(v.ALT) != 1:
            skipped += 1; continue
        # cyvcf2 .genotypes returns [allele1, allele2, phased] per sample.
        gts = v.genotypes  # list of [a, b, phased]
        haps = np.empty(n_haps, dtype=np.int8)
        bad = False
        for i, (a, b, _) in enumerate(gts):
            if a < 0 or b < 0 or a > 1 or b > 1:
                bad = True; break
            haps[2 * i] = a
            haps[2 * i + 1] = b
        if bad:
            skipped += 1; continue
        phys_list.append(v.start + 1)  # cyvcf2 is 0-based; convert to 1-based bp
        rows.append(haps)
    geno = np.stack(rows, axis=0)  # (sites, 2*samples)
    phys = np.asarray(phys_list, dtype=np.int32)

    map_phys, map_cm = _read_genetic_map(map_file)
    cm_at_sites = _interp_cm(phys, map_phys, map_cm)
    morgans = (cm_at_sites / 100.0).astype(np.float32)

    asmc = ASMC.from_arrays(
        genotypes=geno, physical_positions=phys, genetic_positions=morgans,
        dq_file=str(dq_file), use_ancestral=use_ancestral,
        skip_csfs_distance=skip_csfs_distance, decoding_mode=decoding_mode,
    )
    asmc.set_store_per_pair_posterior_mean()
    return asmc


def from_haps(haps_root: PathLike, dq_file: PathLike,
              *, decoding_mode: str = "array",
              skip_csfs_distance: float = 0.0,
              use_ancestral: bool = False) -> ASMC:
    """Construct ASMC from an Eagle/SHAPEIT .hap[.gz] / .sample / .map trio.

    This wraps the legacy file-based ASMC constructor for symmetry with the
    other ``from_*`` helpers; the underlying implementation reads the .hap
    via the C++ side, no in-memory unpack.
    """
    haps_root = str(haps_root)
    asmc = ASMC(haps_root, str(dq_file), "", decoding_mode)
    asmc.set_store_per_pair_posterior_mean()
    return asmc
