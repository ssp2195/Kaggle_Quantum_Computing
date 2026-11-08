"""Full-domain mathematical screening, independent of native kickmix."""
import numpy as np
import pytest
from src.hybrid_model import Config, coefficient_table, full_domain_screen, phase_model

@pytest.mark.parametrize("p,q", [(8,26),(8,28),(9,28),(9,30)])
def test_entire_benchmark_domain(p,q):
    r=full_domain_screen(Config(p,q))
    assert r["passed"] and r["failed_inputs"]==0

@pytest.mark.parametrize("p,q",[(7,30),(8,24)])
def test_reject_insufficient_configs(p,q):
    assert full_domain_screen(Config(p,q))["failed_inputs"]>0

def test_coefficient_lookup():
    for cfg in (Config(8,26),Config(9,28)):
        rows,widths=coefficient_table(cfg)
        assert len(rows)==1<<(cfg.segment_bits+2)
        assert rows[0]==0
        assert all(0 <= row < (1<<sum(widths)) for row in rows)

def test_exact_powers_of_four():
    xs=np.array([1,4,16,64,256,1024,4096,16384,65536,262144,1048576])
    np.testing.assert_allclose(phase_model(Config(8,28),xs),
                               2.**(-np.arange(11)),rtol=0,atol=1e-14)
