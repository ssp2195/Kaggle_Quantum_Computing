import pytest
km=pytest.importorskip("kickmix")
from src.priority1_model import Approximation
from src.priority1_grouped import build_grouped
from src.hybrid_verify import verify_samples,boundary_inputs

@pytest.mark.parametrize("grouping",["ab_c","a_bc"])
@pytest.mark.parametrize("normalization",["serial","barrel"])
def test_grouped_coefficient_circuit(grouping,normalization,tmp_path):
    c=build_grouped(Approximation("chebyshev",2,7,26),
                    (10,12,11),grouping=grouping,normalization=normalization)
    path=tmp_path/"grouped.kmb"
    c.to_file(path,format="kmb")
    restored=km.Circuit.from_file(path,format="kmb")
    assert restored==c
    report=verify_samples(restored,boundary_inputs(random_count=1024))
    assert report["failures"]==0,report
