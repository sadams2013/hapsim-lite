import pytest 

import numpy as np

from hapsim_lite.vcf import generate_sample_ids, generate_call_matrix

def test_generate_sample_ids_length():
    """Sample ID generation creates the number expected"""
    sample_ids = generate_sample_ids(10)
    assert len(sample_ids) == 10

def test_generate_sample_ids_unique():
    """Sample IDs are unique"""
    sample_ids = set(generate_sample_ids(10000))
    assert len(sample_ids) == 10000

def test_generate_call_matrix_haploid():
    """String call matrix generates expected ploidy and number of calls"""
    haplotype_matrix = np.zeros((10, 10000), dtype = np.int8)
    gt_matrix = generate_call_matrix(haplotype_matrix, 1, 10, False)
    assert gt_matrix.shape == (10, 10000)

def test_generate_call_matrix_diploid():
    """String call matrix generates expected ploidy and number of calls"""
    haplotype_matrix = np.zeros((10, 10000), dtype = np.int8)
    gt_matrix = generate_call_matrix(haplotype_matrix, 2, 5, False)
    assert gt_matrix.shape == (5, 10000)

def test_generate_call_matrix_triploid():
    """String call matrix generates expected ploidy and number of calls"""
    haplotype_matrix = np.zeros((15, 10000), dtype = np.int8)
    gt_matrix = generate_call_matrix(haplotype_matrix, 3, 5, False)
    assert gt_matrix.shape == (5, 10000)

