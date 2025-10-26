import pytest

from hapsim_lite.generate import HaplotypeGenerator
from hapsim_lite.population import PopulationData

from . import pop_data

def test_init(pop_data):
    """Initial haplotye simulation produces expected size"""
    hap = HaplotypeGenerator(pop_data, 2)
    assert hap.hap_matrix.shape == (2, 5)

def test_forward_pass(pop_data):
    """Forward pass doesn't alter the shape of the hap_matrix or crash"""
    hap = HaplotypeGenerator(pop_data, 2)
    hap.forward_pass()
    assert hap.hap_matrix.shape == (2, 5)

def test_reverse_pass(pop_data):
    """Reverse pass doesn't alter the shape of the hap_matrix or crash"""
    hap = HaplotypeGenerator(pop_data, 2)
    hap.reverse_pass()
    assert hap.hap_matrix.shape == (2, 5)

def test__get_context(pop_data):
    """Ensure that context windows are properly calculated and do not under/overflow"""
    hap = HaplotypeGenerator(pop_data, 2)
    context = hap._get_context(4, 3, 'left')
    assert context == [1, 2, 3]
    context = hap._get_context(4, 4, 'left')
    assert context == [0, 1, 2, 3]
    context = hap._get_context(4, 4, 'both')
    assert context == [0, 1, 2, 3]
    context = hap._get_context(1, 4, 'both')
    assert context == [0, 2, 3, 4]
    context = hap._get_context(1, 3, 'right')
    assert context == [2, 3, 4]

def test_get_prob_in_window(pop_data):
    """Ensure probability calculations are of right size and reasonable"""
    hap = HaplotypeGenerator(pop_data, 2)
    context = hap._get_context(4, 3, 'left')
    prob = hap.get_prob_in_window(4, context)
    assert prob.shape[0] == 2
    assert prob[0] > 0.0