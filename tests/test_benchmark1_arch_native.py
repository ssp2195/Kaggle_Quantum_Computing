"""Native 22-bit kickmix architectural comparison, coherent 64-shot tests."""
import json
import pytest
km = pytest.importorskip("kickmix")
from src.benchmark1_piecewise_linear import LinearConfig, build, screen
from src.priority1_grouped import build_grouped
from src.priority1_model import Approximation
from src.hybrid_verify import boundary_inputs, verify_samples
from src.hybrid_kickmix import resource_counts


@pytest.mark.parametrize("architecture", ["linear-joint", "linear-streamed", "square-symmetric"])
def test_kickmix_architecture(architecture, tmp_path):
    if architecture.startswith("linear"):
        circ = build(LinearConfig(10, 28, 7),
                     grouping=architecture.split("-")[1])
    else:
        circ = build_grouped(Approximation("chebyshev", 2, 7, 26),
                             drops=(10,12,11), grouping="ab_c",
                             normalization="barrel", square_method="symmetric")
    path = tmp_path / (architecture + ".kmb")
    circ.to_file(path, format="kmb")
    recovered = km.Circuit.from_file(path, format="kmb")
    assert recovered == circ
    result = verify_samples(recovered, boundary_inputs(random_count=1500))
    cost = resource_counts(circ)
    print(json.dumps(dict(architecture=architecture, **cost, **result)))
    assert result["failures"] == 0, result
