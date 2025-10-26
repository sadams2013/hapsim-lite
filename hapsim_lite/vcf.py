"""
Contains functions for generating random sample IDs, vectored translation of
haplotypes into Nploid calls, and writing VCF output to stdout.

This approach has fewer guardrails than using pysam to generate and write the VCF records,
but it is orders-of-magnitude faster.
"""

import uuid

from typing import Optional
from itertools import product
import numpy as np

from .population import PopulationData
from . import VCF_HEADER


def _make_contig_row(contig: str, length: Optional[str] = None) -> str:
    """For a given contig, generate the VCF header row with optional length"""
    return f"##contig=<ID={contig}{',length=' + length if length is not None else ''}>"


def generate_sample_ids(n_samples: int) -> list[str]:
    """Generates a list of random sample IDs"""
    return [uuid.uuid4().hex[0:10] for i in range(0, n_samples)]


def generate_call_matrix(
    haplotype_matrix: np.array, ploidy: int, num_samples: int, unphased: bool
) -> np.array:
    """Reshape simulated haplotypes according to ploidy and sample count,
    returning GT strings for each variant."""
    # reshape to (samples, ploidy, variants)
    gt_shaped = haplotype_matrix.reshape(num_samples, ploidy, haplotype_matrix.shape[1])
    # compute genotype indices
    powers = 2 ** np.arange(ploidy)[::-1]  # [2^(p-1), ..., 1]
    # get the dot product the powers and gt reshaped, which will map to the string GT later
    gt_calls = np.tensordot(
        gt_shaped, powers, axes=([1], [0])
    )  # shape (samples, variants)
    # choose separator and lookup
    gt_delim = "/" if unphased else "|"
    lookup = np.array([gt_delim.join(i) for i in product(["0", "1"], repeat=ploidy)])
    return lookup[gt_calls]


def write_vcf(
    gt_matrix: np.array, population_data: PopulationData, sample_ids: list[str]
) -> None:
    """Writes a VCF to stdout so it can be sorted/bgzipped by bcftools
    (or written to uncompressed VCF)"""
    # slice haplotype matrix into sets of 2 rows
    contig_rows = []
    for contig in population_data.contigs:
        contig_rows.append(_make_contig_row(contig))
    samples = "\t".join(sample_ids)
    header = VCF_HEADER.format(contigs="\n".join(contig_rows), samples=samples)
    print(header)
    for i, (chrom, pos, ref, alt) in enumerate(population_data.variant_info):
        gt = "\t".join(gt_matrix[:, i])
        print(
            f"{chrom}\t{pos}\t{chrom}_{pos}_{ref}_{alt}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}"
        )
