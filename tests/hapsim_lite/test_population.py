import pytest

from . import pop_data

def test_init(pop_data):
    """Population Data init without error and with expected variant index length"""
    assert len(pop_data.variant_index) == 5