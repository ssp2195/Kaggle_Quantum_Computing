"""Pure-Python tests; do not require kickmix."""
from math import sqrt
from itertools import combinations
import numpy as np
import pytest
from src.direct_phase import TOL, dyadic_partition, verify_partition
from src.kickmix_phase import monomial_angles


@pytest.mark.parametrize("degree", [1, 2, 3])
def test_small_domain_exhaustive(degree):
    blocks = dyadic_partition(degree, max_x=255)
    report = verify_partition(blocks, max_x=255)
    assert report["passes_numerical_tolerance"]
    assert report["max_tolerance_fraction"] <= .75
    assert blocks[0].start == 1 and blocks[-1].end == 255


def test_boolean_expansion_and_pruning_bound():
    for b in dyadic_partition(3, max_x=4095):
        terms, dropped = monomial_angles(b, .05)
        assert b.max_ratio + dropped*sqrt(b.end)/TOL < 1.01
        if b.size <= 64:
            for x in range(b.start, b.end+1):
                low = x-b.start
                phase = sum(a for mask, a in terms if low & mask == mask)
                assert abs(phase - 1/sqrt(x)) <= TOL/sqrt(x)+2e-13


def test_parity_identity():
    for n in range(1, 5):
        for value in range(1 << n):
            acc = 0
            for r in range(1, n+1):
                for s in combinations(range(n), r):
                    acc += (-1)**(r-1) * (sum((value >> i) & 1 for i in s)%2)/(1 << (n-1))
            assert abs(acc - (value == (1 << n)-1)) < 1e-12
