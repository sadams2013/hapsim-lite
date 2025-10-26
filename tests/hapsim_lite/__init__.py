import pytest

import numpy as np
from scipy.sparse import csr_matrix

from hapsim_lite.population import PopulationData

@pytest.fixture()
def pop_data():
    """Creates a test instance of the haplotype generator for use in tests"""
    # we have to define all of the paramters needed for the population data class
    variants = [
        "1_1_A_G",
        "1_10_C_T",
        "1_50_A_T",
        "1_70_T_A",
        "1_100_T_A"
        ]
    variant_index = {k: i for i, k in enumerate(variants)}
    variant_info = []
    positions = []
    for v in variants:
        chrom, pos, ref, alt = v.split("_")
        variant_info.append((chrom, int(pos), ref, alt))
        positions.append(int(pos))
    mafs = np.array([0.1, 0.05, 0.01, 0.3, 0.04])
    # make our csr matrix and set variants 1 and 3 to 0.9
    # and 2 and 4 to 0.1
    prob_matrix = csr_matrix(
        ([0.9, 0.1, 0.9, 0.1], ([1, 2, 3, 4], [3, 4, 1, 2])),
            dtype = np.float32)
    yield PopulationData(
        variant_index = variant_index,
        variant_info = variant_info,
        positions = np.array(positions),
        mafs = np.array(mafs), 
        prob_matrix = prob_matrix
    )