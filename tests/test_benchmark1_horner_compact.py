"""Native kickmix in-place Horner with early measurement-erasure normalization."""
import pytest


@pytest.mark.parametrize("guard", [3,4])
def test_compact_horner_with_sampled_coherence(guard,tmp_path):
    km=pytest.importorskip("kickmix")
    from src.benchmark1_horner import HornerConfig,build,screen
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import boundary_inputs,verify_samples
    config=HornerConfig(guard,"mbuc","compact")
    math=screen(config)
    assert math["passes"] and math["worst_error_ratio"]<.9
    c=build(config)
    path=tmp_path/"horner_compact.kmb"
    c.to_file(path,format="kmb")
    recovered=km.Circuit.from_file(path,format="kmb")
    report=verify_samples(recovered,boundary_inputs(random_count=8192))
    stats=resource_counts(c)
    print("COMPACT_HORNER",dict(config=vars(config),resources=stats,
                                mathematical=math,native=report))
    assert report["failures"]==0,report
