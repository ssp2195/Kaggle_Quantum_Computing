"""Test exact proxy combinatorics and native one-term price."""

from src.benchmark1_measure_streaming_cost import assess,raw_phase_term_counts


def test_direct_product_phase_term_counts():
    d=raw_phase_term_counts()
    assert d["intercept_bits"]==27
    assert d["linear_bit_pair_terms"]==179
    assert d["quadratic_unique_terms"]==202
    assert d["quadratic_crossbit_degree3_terms"]==161
    assert d["quadratic_squarebit_degree2_terms"]==41
    assert d["total_nonconstant_raw_phase_terms"]==381


def test_cost_proxy_disclaims_exact_oracle():
    d=assess()
    assert "not a complete valid phase oracle" in d["status"]
    assert d["estimated_rotation_cost_only"]["50"]["rotations_only_cost_if_one_fine_angle_per_raw_term"]>d["target_qc"]


def test_native_measurement_based_phase_term():
    import pytest
    pytest.importorskip("kickmix")
    from src.benchmark1_measure_streaming_cost import native_single_term
    inline=native_single_term("inline")
    mbuc=native_single_term("mbuc")
    assert inline["toffoli_like"]==mbuc["toffoli_like"]
    assert inline["fine_rotations"]==mbuc["fine_rotations"]
    print("STREAMED_MEASURED_TERM",{"inline":inline,"mbuc":mbuc})
