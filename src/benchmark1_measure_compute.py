"""Measurement *during* phase computation: native kickmix resource study.

This is not just ordinary uncomputation: the phase rotation and the
X-basis measurement together implement an unbiased phase-transfer
instrument. The measured temporary need not persist as a logical output.

For g(x)=a AND b, immediately after Z_POW(t,e) the X-measurement
outcome m introduces (-1)^(m*a*b). A classically conditioned CZ(a,b)
corrects it. In kickmix, ccz(a,b,m) becomes CZ_IF.

Compare:
 - coherent: compute AND, rotate, coherent CCX erase;
 - mbuc:     compute AND, rotate, native del_and cleanup;
 - inline:   compute AND, rotate, explicit X measurement + conditional CZ;
 - cubic:    compute staged ANDs, rotate, inline-measure terminal AND,
             correct while intermediate a*b remains live, then del_and.

All four are deterministic. The phase is 2*pi*exponent/2 times g(x),
and the input must retain coherence. Neither postselection nor
Z-measurement of a data-dependent temporary is allowed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json

ANGLE = Fraction(3, 16)  # Z_POW exponent, yields 3/32 turns on target |1>


def phase_gadget(degree: int, method: str):
    """Return a real native kickmix circuit with named live input x."""
    if degree not in (2, 3):
        raise ValueError("degree must be 2 or 3")
    if method not in ("coherent", "mbuc", "inline"):
        raise ValueError("method must be coherent, mbuc or inline")
    import kickmix as km
    b = km.CircuitBuilder()
    x = b.create_quantum_register(degree, name="x")
    temp = b.alloc_qubits(1)[0]
    b.init_and(x[0], x[1], temp)
    controls = (x[0], x[1])
    if degree == 3:
        stage = b.alloc_qubits(1)[0]
        b.init_and(temp, x[2], stage)
        controls = (temp, x[2])
    else:
        stage = temp

    # The measurement is part of phase *generation* here.
    b.z_pow(stage, ANGLE)
    if method == "coherent":
        if degree == 3:
            b.ccx(temp, x[2], stage)
        else:
            b.ccx(x[0], x[1], stage)
    elif method == "mbuc":
        b.del_and(controls[0], controls[1], stage)
    else:
        m = b.alloc_bits(1)[0]
        b.hmr(stage, m)
        b.ccz(controls[0], controls[1], m)
        b.free(m)

    if degree == 3:
        b.free(stage)
        if method == "coherent":
            b.ccx(x[0], x[1], temp)
        else:
            b.del_and(x[0], x[1], temp)
    b.free(temp)
    return b.finish_circuit()


def assess(shots: int = 64):
    """All 2^n computational inputs, multiple stochastic trajectories."""
    import kickmix as km
    from .hybrid_kickmix import resource_counts

    rows = []
    for degree in (2, 3):
        for method in ("coherent", "mbuc", "inline"):
            circuit = phase_gadget(degree, method)
            stats = resource_counts(circuit)
            sim = km.Simulator(batch_size=shots)
            sim.use_same_registers_as(circuit)
            worst = 0.0
            failures = 0
            for x in range(1 << degree):
                sim.clear_for_shot()
                for shot in range(shots):
                    sim.write_within_shot("x", shot, x)
                sim.do(circuit)
                expected = float(ANGLE) / 2 if x == (1 << degree) - 1 else 0.
                for shot in range(shots):
                    if sim.read_within_shot("x", shot, out=int) != x:
                        failures += 1
                    phase = (float(sim.read_shot_phase(shot)) / 2) % 1
                    error = abs((phase - expected + .5) % 1 - .5)
                    worst = max(worst, error)
                    if error > 1e-10:
                        failures += 1
            rows.append(dict(
                degree=degree, method=method,
                verified_trajectories=(1 << degree)*shots,
                failures=failures,
                worst_phase_turns_error=worst,
                **stats
            ))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output")
    args = p.parse_args()
    rows = assess()
    print(json.dumps(rows, indent=2, sort_keys=True))
    if args.output:
        from pathlib import Path
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(rows, indent=2)+"\n")
    if any(r["failures"] for r in rows):
        raise SystemExit("Native instrument left an input or phase error.")


if __name__ == "__main__":
    main()
