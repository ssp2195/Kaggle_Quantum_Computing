import pytest
km=pytest.importorskip("kickmix")
from src.priority1_model import Approximation
from src.priority1_kickmix import build_priority1
from src.hybrid_verify import verify_samples,boundary_inputs

@pytest.mark.parametrize("degree,p",[(2,6),(2,7),(3,4)])
def test_native_full_width_circuits(degree,p,tmp_path):
    c=build_priority1(Approximation("chebyshev",degree,p,26))
    assert len(c.register_data["x"])==22
    dest=tmp_path/f"candidate_{degree}_{p}.kmb"
    c.to_file(dest,format="kmb")
    restored=km.Circuit.from_file(dest,format="kmb")
    assert restored==c
    result=verify_samples(restored,boundary_inputs(random_count=512))
    assert result["failures"]==0,result
