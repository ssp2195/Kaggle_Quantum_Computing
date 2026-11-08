"""Full 22-bit kickmix circuit and measurement trajectory integration tests."""
import pytest
km=pytest.importorskip("kickmix")
from src.hybrid_model import Config
from src.hybrid_kickmix import build_hybrid, resource_counts
from src.hybrid_verify import boundary_inputs, verify_samples

@pytest.mark.parametrize("p,q",[(8,26),(8,28),(9,28)])
def test_full_width_candidate(p,q,tmp_path):
    circuit=build_hybrid(Config(p,q))
    assert set(circuit.register_data)=={"x"}
    assert len(circuit.register_data["x"])==22
    path=tmp_path/f"hybrid_{p}_{q}.kmb"
    circuit.to_file(path,format="kmb")
    reloaded=km.Circuit.from_file(path,format="kmb")
    assert reloaded==circuit
    sample=boundary_inputs(random_count=160)
    report=verify_samples(reloaded,sample)
    assert report["failures"]==0,report
    assert report["verified_shots"]==len(sample)*64
    assert resource_counts(circuit)["qubits"]>=22
