"""Independent native regression for full-verifier v4 winner."""
import hashlib
import pytest

def test_verified_v4_exact_hash_resource_and_coherent_sample(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v4 import build_verified,VERIFIED_QC,VERIFIED_SHA256
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    circ=build_verified()
    stats=resource_counts(circ)
    assert (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])==(86,2652,86)
    assert VERIFIED_QC==283167.97477128694
    f=tmp_path/"v4.kmb"
    circ.to_file(f,format="kmb")
    assert hashlib.sha256(f.read_bytes()).hexdigest()==VERIFIED_SHA256
    report=verify_samples(km.Circuit.from_file(f,format="kmb"),boundary_inputs(random_count=1024))
    assert report["failures"]==0,report

def test_integer_phase_approximation_covers_all_normalized_inputs():
    from src.benchmark1_cost_cooptim import full_screen
    r=full_screen((10,13,12))
    assert r["passes"] and r["failures"]==0 and r["all_normalized_inputs"]==3145728,r
