"""Native kickmix trajectory checks (not a substitute for Kaggle verifier)."""
from __future__ import annotations
from math import sqrt
import random
from .hybrid_model import MAX_X, TOL, NORMAL_BITS


def boundary_inputs(max_x=MAX_X, random_count=256, seed=20261108):
    points = {1, 2, 3, 4, 5, 6, 7, 8, max_x, max_x - 1}
    for p in range(23):
        for center in (1 << p, 3 * (1 << p)):
            for offset in (-2, -1, 0, 1, 2):
                if 1 <= center + offset <= max_x:
                    points.add(center + offset)
    for p in range(8, 12):
        stride = 1 << (20 - p)
        for center in range(stride, max_x + 1, stride * 123):
            for offset in (-1, 0, 1):
                if 1 <= center + offset <= max_x:
                    points.add(center + offset)
    rng = random.Random(seed)
    points.update(rng.randint(1, max_x) for _ in range(random_count))
    return sorted(points)


def verify_samples(circuit, xs, shots=64, print_progress=False):
    """64 trajectories per x, 4 values per SIMD batch (256 trajectories).

    Compare the named input register, every unused qubit and every fine phase.
    Not authoritative: official verifier uses its own fixed-point phase
    reference and private random salt, and analytically prices rotations.
    """
    import kickmix as km
    xs = list(xs) if not isinstance(xs, range) else xs
    if shots != 64:
        raise ValueError("Use 64 shots per input")
    if set(circuit.register_data) != {"x"} or len(circuit.register_data["x"]) != NORMAL_BITS:
        raise ValueError("Expected one named 22-qubit input x")
    if not len(xs):
        return dict(inputs_checked=0, verified_shots=0, failures=0)
    qubits = set(circuit.register_data["x"])
    scratch = km.array(sorted(set(km.q(i) for i in range(circuit.num_qubits)) - qubits))
    sim = km.Simulator(batch_size=256, debug_print=False)
    sim.use_same_registers_as(circuit)
    lane = (1 << 64) - 1
    failures, worst, first_bad = 0, 0.0, []
    for offset in range(0, len(xs), 4):
        group = xs[offset:offset + 4]
        sim.clear_for_shot()
        expected_masks = []
        for bit in range(NORMAL_BITS):
            mask = sum(lane << (i * 64) for i, x in enumerate(group)
                       if x & (1 << bit))
            expected_masks.append(mask)
            sim.write_across_shots(circuit.register_data["x"][bit], mask)
        sim.do(circuit)
        got_masks = sim.read_across_shots("x", out=int)
        free_masks = sim.read_across_shots(scratch, out=int)
        if any(int(a) != b for a, b in zip(got_masks, expected_masks)) or any(int(a) for a in free_masks):
            failures += len(group) * 64
            first_bad.extend(group[:max(0, 12 - len(first_bad))])
        for i, x in enumerate(group):
            target = 1.0 / sqrt(x)
            allowance = TOL * target
            for shot in range(i * 64, (i + 1) * 64):
                # kickmix read_shot_phase is in HALF-turns; divide by two.
                actual = (float(sim.read_shot_phase(shot)) / 2) % 1
                error = abs((actual - target + .5) % 1 - .5)
                ratio = error / allowance
                worst = max(worst, ratio)
                if ratio > 1:
                    failures += 1
                    if len(first_bad) < 12:
                        first_bad.append(x)
        if print_progress and (offset + 4) % 200_000 == 0:
            print(f"checked {min(offset+4,len(xs)):,}/{len(xs):,}; failures={failures}; worst={worst:.6g}", flush=True)
    stats = circuit.operation_statistics()
    return dict(inputs_checked=len(xs), shots_per_x=shots,
                verified_shots=len(xs) * shots, failures=failures,
                worst_error_ratio=worst, first_bad_inputs=first_bad[:12],
                qubits=circuit.num_qubits,
                native_toffoli_like=stats.num_toffoli_like, fine_z_pow=stats.num_z_pow(min_clifford_hierarchy_level=4),
                all_z_pow=stats.num_z_pow())
