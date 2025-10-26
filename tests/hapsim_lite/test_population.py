import pytest

from . import pop_data

def test_init(pop_data):
    assert len(pop_data.variant_index) == 5