"""
Contains PopulationData class which is used for parsing and management of
MAF and LD data
"""

from typing import Self
from collections import defaultdict

import numpy as np
import numpy.typing as npt

from scipy.sparse import csr_matrix


class PopulationData:  # pylint: disable=too-few-public-methods
    """Parse and store population LD and MAF data

    Most efficiently created from Plink2 outputs for
    --r-unphased and --freq using a standard variant
    nomenclature of chrom_pos_ref_alt

    Attributes:
        variant_index (dict): connects variant IDs with their array positions
        variant_info (list): index-consistent variant info tuples (helpful for writing output)
        mafs (np.array): index-consistent 1d array of minor allele frequencies with length N
        positions (np.array): index-consistent 1d array of variant positions with length N
        prob_matrix (csr_matrix): index-consistent sparse array of P(i,j) values derived from MAF
            and signed R with shape (N, N)
    """

    variant_index: dict[str, int]
    variant_info: list[tuple]
    mafs: npt.NDArray[np.float64]
    positions: npt.NDArray[np.int32]
    contigs: set[str]
    prob_matrix: csr_matrix[np.float32]

    def __init__(self, **kwargs) -> None:
        """default init for PopulationData, generally useful for testing only"""
        self.variant_index = None
        self.variant_info = None
        self.mafs = None
        self.positions = None
        self.prob_matrix = None
        self.contigs = None
        for key, value in kwargs.items():
            setattr(self, key, value)

    @classmethod
    def from_plink2_afreq_vcor(
        cls, plink_afreq_path: str, plink_vcor_path: str
    ) -> Self:
        """ideal constructor for PopulationData using outputs from plink2"""
        self = cls()
        self._parse_plink2_afreq(plink_afreq_path)
        self._parse_plink2_vcor(plink_vcor_path)
        return self

    def _parse_plink2_afreq(self, plink_afreq_path: str) -> None:
        """Parses plink2 afreq output"""
        self.variant_index = defaultdict(int)
        self.contigs = set()
        maf_array = []
        pos_array = []
        self.variant_info = []
        with open(plink_afreq_path, "r", encoding="utf-8") as freq_file:
            # skip the first (header) row
            last_pos = 0
            freq_file.readline()
            i = 0
            for row in freq_file:
                tokens = row.strip().split("\t")
                var_id = tokens[1]
                try:
                    chrom, pos, ref, alt = var_id.split("_")
                except ValueError:
                    continue
                self.contigs.add(chrom)
                self.variant_info.append((chrom, int(pos), ref, alt))
                if int(pos) < last_pos:
                    # happens if it is not sorted
                    raise ValueError("Plink Input files are not sorted by position!")
                last_pos = int(pos)
                maf = float(tokens[4])
                self.variant_index[var_id] = i
                maf_array.append(maf)
                pos_array.append(pos)
                i += 1
        self.mafs = np.array(maf_array, dtype=np.float32).clip(1e-9, 1 - 1e-9)
        self.positions = np.array(pos_array, dtype=np.int32)

    def _derive_joint_proba(self, a: int, b: int, r_signed: float) -> float:
        """Uses MAFs for 2 variants (a and b) and their R value to derive P(a,b)"""
        maf_a = self.mafs[a]
        maf_b = self.mafs[b]
        s = np.sqrt(maf_a * (1.0 - maf_a) * maf_b * (1.0 - maf_b))
        return maf_a * maf_b + r_signed * s

    def _parse_plink2_vcor(self, plink_vcor_path: str) -> None:
        """Parses plink vcor file to csr_matrix"""
        # initialize the frequency matrix
        # ld_csr = csr_matrix((num_vars, num_vars), dtype = np.float32)
        ld_data = []
        ld_rows = []
        ld_columns = []
        with open(plink_vcor_path, "r", encoding="utf-8") as ld_file:
            # skip header
            ld_file.readline()
            for row in ld_file:
                tokens = row.strip().split("\t")
                var_id_a_address = self.variant_index[tokens[2]]
                var_id_b_address = self.variant_index[tokens[6]]
                pa_b = self._derive_joint_proba(
                    var_id_a_address, var_id_b_address, float(tokens[8])
                )
                # set both sides of the matrix with the batch
                # so we don't have to manually mirror
                ld_data += [pa_b, pa_b]
                ld_rows += [var_id_a_address, var_id_b_address]
                ld_columns += [var_id_b_address, var_id_a_address]
        num_vars = len(self.mafs)
        self.prob_matrix = csr_matrix(
            (ld_data, (ld_rows, ld_columns)),
            shape=(num_vars, num_vars),
            dtype=np.float32,
        )
