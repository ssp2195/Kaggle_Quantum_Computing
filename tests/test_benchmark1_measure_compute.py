"""Check coherent, measurement-in-computation, and cleanup costs."""

from fractions import Fraction
import math
import pytest

from src.benchmark1_measure_compute import ANGLE, phase_gadget


def test_phase_gadgets_all_native_trajectories():
    km = pytest.importorskip("kickmix")
    from src.benchmark1_measure_compute import assess
    rows = assess()
    assert len(rows) == 6
    assert all(r["failures"] == 0 for r in rows), rows
    assert all(r["worst_phase_turns_error"] < 1e-10 for r in rows)
    for degree in (2, 3):
        d = {r["method"]: r for r in rows if r["degree"] == degree}
        # The inline measurement changes an instrument, not its required
        # nonlinear gate count. The explicit CZ correction is Clifford.
        assert d["inline"]["toffoli_like"] == d["mbuc"]["toffoli_like"]
        assert d["inline"]["toffoli_like"] < d["coherent"]["toffoli_like"]


def test_physical_data_measurement_is_not_coherence_preserving():
    """Z-measuring |00>+|11> destroys off-diagonal coherence."""
    import numpy as np
    bell = np.array([1,0,0,1], dtype=complex)/math.sqrt(2)
    rho = np.outer(bell, bell.conj())
    # Measuring AND(a,b) in Z projects on |00>,|01>,|10> and |11>.
    p0 = np.diag([1,1,1,0])
    p1 = np.diag([0,0,0,1])
    dephased = p0 @ rho @ p0 + p1 @ rho @ p1
    assert abs(rho[0, 3]) > 0.49
    assert abs(dephased[0, 3]) < 1e-12


def test_phase_measurement_must_have_correctable_byproduct():
    """Odd/cubic byproducts aren't quadratic on original input bits."""
    import itertools
    n = 3
    # 1, x_i, x_i*x_j form the full diagonal Clifford phase basis over F2
    columns = [
        [1 for _ in range(1<<n)],
        *[[(x>>i)&1 for x in range(1<<n)] for i in range(n)],
        *[[((x>>i)&1)*((x>>j)&1) for x in range(1<<n)]
          for i in range(n) for j in range(i+1,n)]
    ]
    target = [((x&1)>0 and ((x>>1)&1)>0 and ((x>>2)&1)>0)
              for x in range(1<<n)]
    target = list(map(int, target))
    assert not any(
        all((sum(((mask>>k)&1)*cols[x] for k,cols in enumerate(columns))%2)==target[x]
            for x in range(1<<n))
        for mask in range(1<<len(columns))
    )
