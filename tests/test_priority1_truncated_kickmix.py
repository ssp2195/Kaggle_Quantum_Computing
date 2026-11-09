import pytest
km=pytest.importorskip("kickmix")
from src.priority1_model import Approximation
from src.priority1_truncated import build_truncated,plan
from src.hybrid_verify import boundary_inputs,verify_samples

@pytest.mark.parametrize("family,p,q",[
    ("chebyshev",6,26),("chebyshev-tuned",6,25),
])
def test_native_truncated_phase(family,p,q,tmp_path):
    cfg=Approximation(family,2,p,q)
    circuit=build_truncated(cfg)
    dest=tmp_path/"candidate.kmb"
    circuit.to_file(dest,format="kmb")
    read=km.Circuit.from_file(dest,format="kmb")
    assert read==circuit
    outcome=verify_samples(read,boundary_inputs(random_count=1024))
    assert outcome["failures"]==0,outcome
