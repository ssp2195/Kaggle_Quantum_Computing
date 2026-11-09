"""End-to-end numeric and native kickmix 8-bit XAG/Horner tests."""
import pytest
from src.benchmark1_horner_xag import XAGHornerConfig,screen


@pytest.mark.parametrize("g,passes",[(4,True),(5,True),(6,True)])
def test_full_domain_3133452_inputs(g,passes):
    r=screen(XAGHornerConfig(g,"xag"))
    assert r["checked_inputs"]==3_133_452
    assert r["passes"]==passes, r
    if g>=5: assert r["worst_error_ratio"]<.86


@pytest.mark.parametrize("lookup",["xag","qrom"])
def test_native_full_22bit_kickmix(lookup,tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_horner_xag import build
    from src.hybrid_verify import verify_samples,boundary_inputs
    from src.hybrid_kickmix import resource_counts
    c=build(XAGHornerConfig(5,lookup,"mbuc"))
    path=tmp_path/(lookup+".kmb")
    c.to_file(path,format="kmb")
    recovered=km.Circuit.from_file(path,format="kmb")
    report=verify_samples(recovered,boundary_inputs(random_count=8192))
    stats=resource_counts(c)
    print("NATIVE_XAG_HORNER",dict(lookup=lookup,stats=stats,check=report))
    assert report["failures"]==0,report
