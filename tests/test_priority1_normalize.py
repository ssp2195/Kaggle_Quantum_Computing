import pytest
km=pytest.importorskip("kickmix")
from src.priority1_normalize import normalize_barrel,denormalize_barrel
from src.hybrid_verify import boundary_inputs,verify_samples

def test_normalizer_roundtrip(tmp_path):
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    preds,scale=normalize_barrel(b,x)
    denormalize_barrel(b,x,preds,scale)
    circuit=b.finish_circuit()
    result=verify_samples(circuit,boundary_inputs(random_count=4096))
    assert result["failures"]==0,result
