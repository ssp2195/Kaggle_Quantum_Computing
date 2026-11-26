"""Regression for fully verified 5.65x one-lookup circuit.

The first-generation 4.83x reference is intentionally retained.
"""
import hashlib
import pytest

def test_verified_v2_exact_rebuild_and_sampled_coherence(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_single_lookup_verified import build,VERIFIED_SHA256,VERIFIED_COST
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import verify_samples,boundary_inputs
    circuit=build()
    stats=resource_counts(circuit)
    assert (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])==(112,2658,88)
    assert VERIFIED_COST==367309.6152203282
    filename=tmp_path/"verified_v2.kmb"
    circuit.to_file(filename,format="kmb")
    assert hashlib.sha256(filename.read_bytes()).hexdigest()==VERIFIED_SHA256
    report=verify_samples(km.Circuit.from_file(filename,format="kmb"),
                          boundary_inputs(random_count=4096))
    assert report["failures"]==0,report
