"""Full-scale mathematical and native signed-digit phase verification."""
import pytest
from src.benchmark1_centered_horner import CenteredConfig,screen


@pytest.mark.parametrize("p,g,expected",[(6,3,True),(6,4,True),(7,3,True),(7,4,True)])
def test_all_3133452_inputs(p,g,expected):
    r=screen(CenteredConfig(p,g))
    assert r["checked_inputs"]==3_133_452
    assert r["passes"]==expected,r
    assert r["worst_error_ratio"]<.9,r


@pytest.mark.parametrize("p,g,lookup",[(6,4,"qrom"),(6,4,"xag"),(7,3,"qrom")])
def test_native_signed_horner(tmp_path,p,g,lookup):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_centered_horner import build
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    cfg=CenteredConfig(p,g,lookup)
    c=build(cfg)
    path=tmp_path/"centered.kmb"
    c.to_file(path,format="kmb")
    copy=km.Circuit.from_file(path,format="kmb")
    report=verify_samples(copy,boundary_inputs(random_count=4096))
    stats=resource_counts(copy)
    print("CENTERED_NATIVE",dict(config=vars(cfg),stats=stats,report=report))
    assert report["failures"]==0,report
