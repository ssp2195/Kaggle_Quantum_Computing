import pytest
km=pytest.importorskip("kickmix")
from src.priority1_model import Approximation
from src.priority1_streamed import build_streamed
from src.hybrid_verify import boundary_inputs,verify_samples

@pytest.mark.parametrize("drops",[(10,12,11),(10,11,12)])
def test_coefficient_streaming(drops,tmp_path):
    cfg=Approximation("chebyshev",2,7,26)
    circuit=build_streamed(cfg,drops)
    assert list(circuit.register_data)==["x"]
    assert len(circuit.register_data["x"])==22
    path=tmp_path/"streamed.kmb"
    circuit.to_file(path,format="kmb")
    reread=km.Circuit.from_file(path,format="kmb")
    assert reread==circuit
    report=verify_samples(reread,boundary_inputs(random_count=2048))
    assert report["failures"]==0,report
