"""Exact math, Boolean-don't-care, SHA, and coherent tests for v7."""
import hashlib
import pytest

def test_unreachable_address_boolean_extension():
    import numpy as np
    from src.benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
    from src.benchmark1_cost_cooptim import coefficients,full_screen
    original=coefficients((10,13,12))
    optimized,packed=optimized_coefficients("affine")
    assert np.array_equal(original[128:],optimized[128:])
    assert proxy(packed)["top_interaction_monomials"]==0
    assert proxy(packed)["active_nontrivial_prefixes"]==374
    assert full_screen((10,13,12))["passes"]

def test_v7_native_fingerprint_and_exact_coherence(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v7 import (
      build_verified,VERIFIED_SHA256,VERIFIED_QC,EXPECTED_RESOURCES)
    from src.hybrid_verify import verify_samples,boundary_inputs
    from src.hybrid_kickmix import resource_counts
    c=build_verified()
    stats=resource_counts(c)
    assert (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])==EXPECTED_RESOURCES
    assert VERIFIED_QC==256476.78036771005
    out=tmp_path/"original-v7.kmb"
    c.to_file(out,format="kmb")
    assert hashlib.sha256(out.read_bytes()).hexdigest()==VERIFIED_SHA256
    sample=verify_samples(km.Circuit.from_file(out,format="kmb"),
        boundary_inputs(random_count=256))
    assert sample["failures"]==0,sample

def test_prior_pinned_winners_preserved():
    from src.benchmark1_best_verified_v6 import VERIFIED_QC as qc6,VERIFIED_SHA256 as sha6
    from src.benchmark1_best_verified_v5 import VERIFIED_QC as qc5
    assert qc6==268172.78036771005 and qc5==272300.78036771005
    assert len(sha6)==64
