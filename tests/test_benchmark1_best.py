"""Pinned complete Benchmark I circuit must rebuild exactly and remain coherent."""
from pathlib import Path
import hashlib
import pytest


def test_best_rebuild_and_representative_coherence(tmp_path):
    km = pytest.importorskip("kickmix")
    from src.benchmark1_best import build, VERIFIED_SHA256
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import verify_samples, boundary_inputs
    circuit = build()
    stats = resource_counts(circuit)
    assert (stats["qubits"], stats["toffoli_like"], stats["fine_rotations"]) == (112, 3218, 88)
    out = tmp_path / "best_invsqrt.kmb"
    circuit.to_file(out, format="kmb")
    assert hashlib.sha256(out.read_bytes()).hexdigest() == VERIFIED_SHA256
    recovered = km.Circuit.from_file(out, format="kmb")
    report = verify_samples(recovered, boundary_inputs(random_count=8192))
    assert report["failures"] == 0, report
    assert report["verified_shots"] == report["inputs_checked"] * 64
