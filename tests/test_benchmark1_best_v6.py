"""Regression tests for the independently full-verified v6 oracle."""
import hashlib
import pytest

def test_full_domain_binary_normalization():
    from src.benchmark1_binary_block_normalizer import numeric_normalization_check
    r=numeric_normalization_check()
    assert r["passes"] and r["checked"]==3133452 and r["max_shift"]==10

def test_unchanged_integer_phase_approximation():
    from src.benchmark1_cost_cooptim import full_screen
    r=full_screen((10,13,12))
    assert r["passes"] and r["failures"]==0 and r["all_normalized_inputs"]==3145728

def test_verified_v6_hash_resources_and_coherent_sample(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v6 import (
        build_verified,VERIFIED_SHA256,VERIFIED_QC,EXPECTED_RESOURCES)
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    assert VERIFIED_QC<272300.78036771005
    circuit=build_verified()
    r=resource_counts(circuit)
    assert (r["qubits"],r["toffoli_like"],r["fine_rotations"])==EXPECTED_RESOURCES
    out=tmp_path/"v6.kmb"
    circuit.to_file(out,format="kmb")
    assert hashlib.sha256(out.read_bytes()).hexdigest()==VERIFIED_SHA256
    report=verify_samples(km.Circuit.from_file(out,format="kmb"),
        boundary_inputs(random_count=256))
    assert report["failures"]==0,report

def test_v5_sha_pinned_fallback():
    from src.benchmark1_best_verified_v5 import VERIFIED_QC,VERIFIED_SHA256
    assert VERIFIED_QC==272300.78036771005
    assert len(VERIFIED_SHA256)==64
