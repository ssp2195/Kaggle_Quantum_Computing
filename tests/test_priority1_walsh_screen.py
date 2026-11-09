"""Walsh direct phase screening and exact FWHT transform inversion."""
import numpy as np
from src.priority1_walsh_screen import walsh_hadamard, screen


def test_walsh_involution():
    rng = np.random.default_rng(20261109)
    x = rng.standard_normal(256)
    np.testing.assert_allclose(walsh_hadamard(walsh_hadamard(x)) / 256,
                               x, atol=1e-13)


def test_full_normalized_mantissa_screen():
    rows = screen(bits=20, terms=(1000, 1700, 1800, 2000))
    by_count = {row["terms"]: row for row in rows}
    assert all(row["normalized_values"] == (1 << 20) for row in rows)
    assert not by_count[1000]["numerical_pass"]
    assert not by_count[1700]["numerical_pass"]
    assert by_count[1800]["numerical_pass"]
    assert by_count[2000]["numerical_pass"]
