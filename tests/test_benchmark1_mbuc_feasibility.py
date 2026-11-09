"""Actual native kickmix MBUC causal correction tests."""
import pytest


def test_staged_cubic_cleanup():
    pytest.importorskip("kickmix")
    from src.benchmark1_mbuc_feasibility import compare
    result=compare()
    assert result["mbuc"]["shots_verified"]==512
    assert result["coherent"]["shots_verified"]==512
    assert result["mbuc"]["toffoli_like"]==2
    assert result["coherent"]["toffoli_like"]==4
    assert result["mbuc"]["qubits"]==result["coherent"]["qubits"]
