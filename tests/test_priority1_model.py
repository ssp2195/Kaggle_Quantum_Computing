import pytest
import numpy as np
from src.priority1_model import Approximation, coefficient_matrix, pack_coefficients, full_domain_screen

@pytest.mark.parametrize("family,degree,p,q",[
    ("chebyshev",2,6,26),
    ("chebyshev",2,7,26),
    ("chebyshev",3,4,26),
    ("chebyshev",3,5,26),
    ("minimax",2,7,26),
])
def test_full_three_million_domain(family,degree,p,q):
    c=Approximation(family,degree,p,q)
    report=full_domain_screen(c)
    assert report["passes"],report
    table,widths=pack_coefficients(c)
    assert len(table)==(4<<p)
    assert all(0<=v<2**sum(widths) for v in table)

@pytest.mark.parametrize("family,degree,p,q",[
    ("minimax",2,6,26),
    ("minimax",3,4,26),
])
def test_quantized_minimax_must_be_rejected_when_it_fails(family,degree,p,q):
    report=full_domain_screen(Approximation(family,degree,p,q))
    assert report["failures"]>0 and not report["passes"]

def test_polynomial_signs_and_widths():
    for degree,p in [(2,6),(3,4)]:
        cfg=Approximation("chebyshev",degree,p,26)
        mat=coefficient_matrix(cfg)
        assert np.all(mat>=0)
        assert len(pack_coefficients(cfg)[1])==degree+1
