import pytest
from src.priority1_model import Approximation,full_domain_screen
from src.priority1_truncated import omitted_partial_products, choose_drop, plan

def test_omission_enumeration():
    for na,nb in [(3,4),(5,5)]:
        for a in range(1<<na):
            for b in range(1<<nb):
                for drop in range(na+nb):
                    omitted=sum(((a>>i)&1)*((b>>j)&1)*(1<<(i+j))
                        for i in range(na) for j in range(nb) if i+j<drop)
                    assert omitted<=omitted_partial_products(na,nb,drop)
                    assert (a*b-omitted)%(1<<drop)==0

@pytest.mark.parametrize("family,p,q",[
    ("chebyshev",6,26),("chebyshev-tuned",6,25),("chebyshev",7,26)
])
def test_analytical_error_budget(family,p,q):
    cfg=Approximation(family,2,p,q)
    approx=full_domain_screen(cfg)
    params=plan(cfg)
    assert approx["passes"],approx
    # Triangle inequality: mathematical approximation + omitted products
    # (all expressed relative to the target 2**-24/sqrt(x)).
    assert approx["worst_error_ratio"]+params["addition_error_bound"] < 1,params
