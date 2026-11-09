"""Exhaustive classical Horner and actual native kickmix comparisons."""
import pytest
from src.benchmark1_horner import HornerConfig, screen


@pytest.mark.parametrize("guard,passes",[(2,False),(3,True),(4,True)])
def test_all_3133452_inputs_horner(guard,passes):
    result=screen(HornerConfig(guard=guard))
    assert result["checked_inputs"]==3_133_452
    assert result["passes"]==passes, result
    if passes: assert result["worst_error_ratio"]<.90


@pytest.mark.parametrize("guard,cleanup",[(3,"mbuc"),(4,"mbuc"),(4,"coherent")])
def test_native_horner_phase(guard,cleanup,tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_horner import build
    from src.hybrid_verify import boundary_inputs,verify_samples
    from src.hybrid_kickmix import resource_counts
    cfg=HornerConfig(guard=guard,cleanup=cleanup)
    circuit=build(cfg)
    path=tmp_path / f"horner_g{guard}_{cleanup}.kmb"
    circuit.to_file(path,format="kmb")
    loaded=km.Circuit.from_file(path,format="kmb")
    report=verify_samples(loaded,boundary_inputs(random_count=4096))
    stats=resource_counts(circuit)
    print("HORNER_NATIVE",dict(config=vars(cfg),resources=stats,verification=report))
    assert report["failures"]==0,report
