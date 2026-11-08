"""Native kickmix roundtrip and small-domain trajectory test."""
import pytest
km = pytest.importorskip("kickmix")
from src.direct_phase import dyadic_partition
from src.kickmix_phase import build_circuit
from src.verify import verify_kmb


def test_build_serialize_simulate(tmp_path):
    circuit = build_circuit(dyadic_partition(2, max_x=31), input_bits=5)
    path = tmp_path / "toy.kmb"
    circuit.to_file(path, format="kmb")
    assert km.Circuit.from_file(path, format="kmb") == circuit
    result = verify_kmb(str(path), max_x=31, shots_per_x=64,
                        competition_contract=False)
    assert result["failures"] == 0
    assert result["worst_error_ratio"] < 1
