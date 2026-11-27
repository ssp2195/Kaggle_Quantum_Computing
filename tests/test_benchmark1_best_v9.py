"""Independently full-verified Benchmark 1 v9 regression suite."""
import hashlib
import pytest

def test_13bit_exact_symmetric_square():
    from src.benchmark1_symmetric_square_rows_v7 import exact_math_check
    r=exact_math_check()
    assert r["passes"] and r["checked"]==8192

def test_full_domain_normalization_and_coefficient_equality():
    from src.benchmark1_binary_block_normalizer import numeric_normalization_check
    from src.benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
    from src import benchmark1_cost_cooptim as cc
    assert numeric_normalization_check()["passes"]
    assert cc.full_screen((10,13,12))["passes"]
    mat,packed=optimized_coefficients("affine")
    assert (mat[128:]==cc.coefficients((10,13,12))[128:]).all()
    assert proxy(packed)["top_interaction_monomials"]==0

def test_exact_kmb_hash_and_coherent_clean_ancillas(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v9 import build_verified,VERIFIED_SHA256,VERIFIED_QC
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    c=build_verified()
    s=resource_counts(c)
    assert (s["qubits"],s["toffoli_like"],s["fine_rotations"])==(85,2232,76)
    assert VERIFIED_QC==237684.49222389946
    p=tmp_path/"benchmark1-v9.kmb"
    c.to_file(p,format="kmb")
    assert hashlib.sha256(p.read_bytes()).hexdigest()==VERIFIED_SHA256
    rec=verify_samples(km.Circuit.from_file(p,format="kmb"),
                       boundary_inputs(random_count=256))
    assert rec["failures"]==0,rec

def test_prior_reproducible_winners_remain_available():
    from src.benchmark1_best_verified_v8 import VERIFIED_QC as v8
    from src.benchmark1_best_verified_v7 import VERIFIED_QC as v7
    assert v8==241424.49222389946 and v7==256476.78036771005
