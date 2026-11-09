"""The novel phase and symmetric-square architectures use native kickmix."""
import pytest
import numpy as np
from src.benchmark1_piecewise_linear import LinearConfig, screen

def test_linear_full_domain_screen():
    report = screen(LinearConfig(10, 28, 7))
    assert report["checked_inputs"] == 3_133_452
    assert report["failed_inputs"] == 0, report
    assert report["worst_error_ratio"] < .85

def test_linear_rejects_inaccurate_cut():
    report = screen(LinearConfig(10, 28, 10))
    assert report["failed_inputs"] > 0

def test_symmetric_square_algebra_exact():
    from random import Random
    rng = Random(4503)
    for n in range(2, 15):
        for drop in range(0, 2*n):
            for _ in range(8):
                x = rng.randrange(1 << n)
                ordered = sum(((x >> i) & 1) * ((x >> j) & 1) * (1 << (i+j-drop))
                              for i in range(n) for j in range(n)
                              if i+j >= drop)
                symmetric = (
                    sum(((x >> i) & 1) * (1 << (2*i-drop))
                        for i in range(n) if 2*i >= drop)
                    + sum(((x >> i) & 1) * ((x >> j) & 1) *
                          (1 << (i+j+1-drop))
                          for i in range(n) for j in range(i+1,n)
                          if i+j >= drop)
                )
                assert ordered == symmetric
