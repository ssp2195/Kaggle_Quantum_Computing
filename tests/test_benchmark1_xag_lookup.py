"""Boolean XAG lookup truth table and native coherent measurement checks."""
import pytest
from src.benchmark1_xag_lookup import anf_coefficients


def test_binary_truth_table_and_anf_inverse():
    for n in range(1, 10):
        table = [((v * 109 + (v >> 2) * 97) ^ 127) & 65535
                 for v in range(1 << n)]
        coeff = anf_coefficients(table, n)
        for x in range(1 << n):
            reconstructed = 0
            m = x
            while True:
                reconstructed ^= coeff[m]
                if m == 0:
                    break
                m = (m - 1) & x
            assert reconstructed == table[x]


def test_invalid_table_rejected():
    with pytest.raises(ValueError):
        anf_coefficients([1, 2, 3], 2)


def test_kickmix_anf_reversible_roundtrip():
    km = pytest.importorskip("kickmix")
    from src.benchmark1_xag_lookup import comparison_circuits
    from src.priority1_model import Approximation, pack_coefficients
    table, widths = pack_coefficients(
        Approximation("chebyshev-tuned", 2, 6, 25))
    assert len(table) == 256
    generic, xag = comparison_circuits(table, 8, sum(widths))
    for c in (generic, xag):
        assert set(c.register_data) == {"address", "output"}
        sim = km.Simulator(batch_size=64)
        sim.use_same_registers_as(c)
        for x in (0, 1, 2, 3, 15, 31, 63, 64, 127, 128, 192, 255):
            sim.clear_for_shot()
            for shot in range(64):
                sim.write_within_shot("address", shot, x)
            sim.do(c)
            for shot in range(64):
                assert sim.read_within_shot("address", shot, out=int) == x
                assert sim.read_within_shot("output", shot, out=int) == 0
                assert abs(float(sim.read_shot_phase(shot))) < 1e-10
    base = generic.operation_statistics()
    candidate = xag.operation_statistics()
    print("NATIVE_QROM_COMPARISON",
          dict(reference_qubits=generic.num_qubits,
               reference_toffolis=base.num_toffoli_like,
               xag_qubits=xag.num_qubits,
               xag_toffolis=candidate.num_toffoli_like))
