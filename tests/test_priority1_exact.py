import numpy as np
import pytest
from src.priority1_exact import truncated_product,full_cut_screen
from src.priority1_model import Approximation

def test_exact_partial_product_identity():
    for a in range(32):
        for b in range(32):
            for drop in range(7):
                expected=sum(
                    ((a>>i)&1)*((b>>j)&1)*(1<<(i+j-drop))
                    for i in range(5) for j in range(5) if i+j>=drop)
                assert int(truncated_product(np.array([a]),np.array([b]),drop)[0])==expected

@pytest.mark.parametrize("drops",[(10,8,15),(10,12,11),(10,11,12)])
def test_all_3133452_integer_inputs(drops):
    report=full_cut_screen(Approximation("chebyshev",2,7,26),drops)
    assert report["passes"] and report["worst_error_over_tolerance"]<.85,report
