"""Isolated, exact kickmix resource profiling of the current best architecture.

Each stage is compiled as its own real native kickmix circuit. Gate/rotation
counts can be added to identify dominant operations; the isolated stage's Q
is NOT the peak Q of the full circuit (workspace lifetimes overlap).
"""
import argparse
import json
from pathlib import Path


def profile():
    import kickmix as km
    from .priority1_model import Approximation, coefficient_matrix
    from .priority1_truncated import _mul_cut
    from .priority1_normalize import normalize_barrel, denormalize_barrel
    from .benchmark1_symmetric_square import square_symmetric
    from .hybrid_kickmix import _phase
    cfg = Approximation("chebyshev", 2, 7, 26)
    mat = coefficient_matrix(cfg)
    widths = [max(1, int(mat[:, i].max()).bit_length()) for i in range(3)]
    lo, q = cfg.low_bits, cfg.precision
    dl, ds, dq = 10, 12, 11
    rows = {}

    def measure(name, create):
        circuit = create()
        stats = circuit.operation_statistics()
        rows[name] = dict(toffoli_like=stats.num_toffoli_like,
                          fine_rotations=stats.num_z_pow(min_clifford_hierarchy_level=4),
                          t_like=stats.num_t_like,
                          instructions=stats.num_instructions,
                          isolated_qubit_peak=circuit.num_qubits)

    def normalization():
        b = km.CircuitBuilder()
        x = b.create_quantum_register(22, name="x")
        flags, shift = normalize_barrel(b, x)
        denormalize_barrel(b, x, flags, shift)
        return b.finish_circuit()
    measure("normalization_and_restore", normalization)

    def lookup_block(indices):
        b = km.CircuitBuilder()
        x = b.create_quantum_register(22, name="x")
        vals = []
        for row in mat:
            n, pos = 0, 0
            for i in indices:
                n |= int(row[i]) << pos
                pos += widths[i]
            vals.append(n)
        w = sum(widths[i] for i in indices)
        store = b.alloc_qubits(w)
        b.init_lookup(table=vals, address=x[lo:], target=store)
        b.del_lookup(table=vals, address=x[lo:], target=store)
        b.free(store)
        return b.finish_circuit()
    measure("AB_lookup_and_unlookup", lambda: lookup_block((0, 1)))
    measure("C_lookup_and_unlookup", lambda: lookup_block((2,)))

    def constant():
        b = km.CircuitBuilder()
        a = b.create_quantum_register(widths[0], name="a")
        shift = b.create_quantum_register(4, name="shift")
        _phase(b, a, shift, q, +1)
        return b.finish_circuit()
    measure("constant_phase", constant)

    def linear():
        b = km.CircuitBuilder()
        beta = b.create_quantum_register(widths[1], name="beta")
        u = b.create_quantum_register(lo, name="u")
        shift = b.create_quantum_register(4, name="shift")
        product = b.alloc_qubits(len(beta) + lo - dl)
        b.reset(product)
        _mul_cut(b, beta, u, product, dl)
        _phase(b, product, shift, q + lo - dl, -1)
        _mul_cut(b, beta, u, product, dl, inverse=True)
        b.free(product)
        return b.finish_circuit()
    measure("linear_phase_with_reversible_product", linear)

    def quadratic(method="rows"):
        b = km.CircuitBuilder()
        gamma = b.create_quantum_register(widths[2], name="gamma")
        u = b.create_quantum_register(lo, name="u")
        shift = b.create_quantum_register(4, name="shift")
        square = b.alloc_qubits(2*lo-ds)
        b.reset(square)
        if method == "rows":
            _mul_cut(b, u, u, square, ds)
        else:
            square_symmetric(b, u, square, ds)
        product = b.alloc_qubits(len(gamma) + len(square) - dq)
        b.reset(product)
        _mul_cut(b, gamma, square, product, dq)
        _phase(b, product, shift, q+2*lo-ds-dq, +1)
        _mul_cut(b, gamma, square, product, dq, inverse=True)
        b.free(product)
        if method == "rows":
            _mul_cut(b, u, u, square, ds, inverse=True)
        else:
            square_symmetric(b, u, square, ds, inverse=True)
        b.free(square)
        return b.finish_circuit()
    measure("quadratic_phase_rows", lambda: quadratic("rows"))
    measure("quadratic_phase_symmetric", lambda: quadratic("symmetric"))

    totals = dict(
        baseline_stage_toffolis=sum(v["toffoli_like"] for k, v in rows.items()
                                    if k != "quadratic_phase_symmetric"),
        symmetric_stage_toffolis=sum(v["toffoli_like"] for k, v in rows.items()
                                     if k != "quadratic_phase_rows"),
    )
    return dict(coefficient_widths=widths, segment_bits=cfg.segment_bits,
                stage_costs=rows, totals=totals)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    data = profile()
    print(json.dumps(data, indent=2))
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
