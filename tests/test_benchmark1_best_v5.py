"""Native and exhaustive regression tests for the fully verified original v5."""
import hashlib
import pytest

def test_padding_minimality_for_all_competition_inputs():
    from src.benchmark1_zero_residual_pad import PADS,full_pad_certificate
    r=full_pad_certificate()
    assert r["passes"] and r["checked"]==3133452,r
    assert (r["min_linear_pad"],r["min_quadratic_pad"])==(5,5)
    assert PADS==(10,5,5)
    assert all(item["required_linear_pad"]<=5 and
               item["required_quadratic_pad"]<=5
               for item in r["max_bit_heights_by_shift"])

def test_pinned_native_kmb_and_coherent_phase(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v5 import (
       build_verified,VERIFIED_SHA256,VERIFIED_QC,EXPECTED_RESOURCES)
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import verify_samples,boundary_inputs
    circuit=build_verified()
    s=resource_counts(circuit)
    assert (s["qubits"],s["toffoli_like"],s["fine_rotations"])==EXPECTED_RESOURCES
    assert VERIFIED_QC<283167.97477128694
    out=tmp_path/"verified_v5.kmb"
    circuit.to_file(out,format="kmb")
    assert hashlib.sha256(out.read_bytes()).hexdigest()==VERIFIED_SHA256
    loaded=km.Circuit.from_file(out,format="kmb")
    report=verify_samples(loaded,boundary_inputs(random_count=256))
    assert report["failures"]==0,report

def test_v4_incumbent_remains_available():
    from src.benchmark1_best_verified_v4 import VERIFIED_SHA256,VERIFIED_QC
    assert VERIFIED_QC==283167.97477128694
    assert len(VERIFIED_SHA256)==64
