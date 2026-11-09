"""Test integer coefficient tuning against every Benchmark-I input."""
from src.priority1_integer_tuned import tuned_matrix, DROPS
from src.priority1_model import Approximation
from src.priority1_exact import full_cut_screen


def test_exhaustive_integer_optimized_phase():
    cfg = Approximation("chebyshev-cut-tuned", 2, 7, 25)
    matrix = tuned_matrix()
    assert matrix.shape == (512, 3)
    assert matrix[128, 0] > 0
    r = full_cut_screen(cfg, DROPS)
    assert r["checked_inputs"] == 3_133_452
    assert r["failed_inputs"] == 0, r
    assert r["worst_error_over_tolerance"] < .96, r
