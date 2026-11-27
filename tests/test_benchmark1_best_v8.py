"""Proof, SHA, and coherence regressions for original 8.60x Benchmark I v8."""
import hashlib
import pytest

def test_symmetric_square_all_13_bit_words():
    from src.benchmark1_symmetric_square_rows_v7 import exact_math_check
    v=exact_math_check()
    assert v["passes"] and v["checked"]==8192 and v["rows"]==12

def test_intercept_overflow_is_only_whole_turn():
    from src.benchmark1_phase_pad9_symmetric import overflow_certificate
    v=overflow_certificate()
    assert v["passes"]
    assert v["max_extra_relative_error_ratio"]<.001

def test_v8_sha_and_exact_coherent_cleanup(tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_best_verified_v8 import build_verified,VERIFIED_SHA256,VERIFIED_QC
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    c=build_verified()
    st=resource_counts(c)
    assert (st["qubits"],st["toffoli_like"],st["fine_rotations"])==(85,2276,76)
    assert VERIFIED_QC<256476.78036771005
    p=tmp_path/"v8.kmb"
    c.to_file(p,format="kmb")
    assert hashlib.sha256(p.read_bytes()).hexdigest()==VERIFIED_SHA256
    rep=verify_samples(km.Circuit.from_file(p,format="kmb"),
                       boundary_inputs(random_count=256))
    assert rep["failures"]==0,rep

def test_v7_v6_preserved():
    from src.benchmark1_best_verified_v7 import VERIFIED_QC as q7
    from src.benchmark1_best_verified_v6 import VERIFIED_QC as q6
    assert q7==256476.78036771005 and q6==268172.78036771005
