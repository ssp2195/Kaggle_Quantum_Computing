"""Real C++ kickmix generation, round trip and coherent phase trajectories."""
import pytest
pytest.importorskip("kickmix")
from src.priority1_model import Approximation
from src.priority1_truncated import build_truncated
from src.hybrid_kickmix import resource_counts
from src.hybrid_verify import verify_samples, boundary_inputs


def test_tuned_cut_native_kickmix(tmp_path):
    cfg = Approximation("chebyshev-cut-tuned", 2, 7, 25)
    circuit = build_truncated(cfg, explicit_drops=(9, 10, 9))
    result = resource_counts(circuit)
    assert result["qubits"] < 165, result
    assert result["toffoli_like"] < 4490, result
    output = tmp_path / "quantized_minimax.kmb"
    circuit.to_file(output, format="kmb")
    import kickmix as km
    recovered = km.Circuit.from_file(output, format="kmb")
    report = verify_samples(recovered, boundary_inputs(random_count=8192))
    assert report["failures"] == 0, report
