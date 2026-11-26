"""Byte-exact full-verifier reproduction and native phase coherency regression."""
import hashlib
import pytest

def test_best_joint_v3_exact_rebuild_and_coherence(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_joint_v3 import build_verified,VERIFIED_SHA256,VERIFIED_COST
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    c=build_verified()
    stats=resource_counts(c)
    assert (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])==(89,2764,87)
    assert VERIFIED_COST==302641.04424652964
    out=tmp_path/"best_joint_v3.kmb"
    c.to_file(out,format="kmb")
    assert hashlib.sha256(out.read_bytes()).hexdigest()==VERIFIED_SHA256
    report=verify_samples(km.Circuit.from_file(out,format="kmb"),
                          boundary_inputs(random_count=1024))
    assert report["failures"]==0,report
