"""Development verifier for all Benchmark I trajectory conditions.

The competition's official verifier remains authoritative for submission
acceptance, salted sampling, rotation synthesis and Q*C score.
"""
import math
from .direct_phase import INPUT_BITS, MAX_X, TOL


def verify_kmb(path, *, max_x=MAX_X, shots_per_x=64, competition_contract=True):
    import kickmix as km
    circuit = km.Circuit.from_file(path, format="kmb")
    regs = circuit.register_data
    if set(regs) != {"x"}:
        raise ValueError(f"Expected only named register x, found {set(regs)}")
    if competition_contract and (len(regs["x"]) != INPUT_BITS
                                 or max_x != MAX_X or shots_per_x != 64):
        raise ValueError("Full contract: 22-bit x, 3,133,452 inputs, 64 shots/x")
    if max_x >= (1 << len(regs["x"])):
        raise ValueError("Input register too narrow")
    scratch = km.array(sorted(set(km.q(i) for i in range(circuit.num_qubits))
                              - set(regs["x"])))
    sim = km.Simulator(batch_size=shots_per_x)
    sim.use_same_registers_as(circuit)
    failures, worst = 0, 0.0
    for x in range(1, max_x+1):
        sim.clear_for_shot()
        for shot in range(shots_per_x):
            sim.write_within_shot("x", shot, x)
        sim.do(circuit)
        # Integer-sqrt reference: floor(2**65 / sqrt(x)) / 2**65.
        expected = math.isqrt((1 << 130)//x) / float(1 << 65)
        tolerance = TOL*expected
        # kickmix reports the complete phase in half-turns (Z_POW exponent).
        for shot in range(shots_per_x):
            read_x = sim.read_within_shot("x", shot, out=int)
            clean = sim.read_within_shot(scratch, shot, out=int) == 0
            phase = (float(sim.read_shot_phase(shot)) / 2.0) % 1.0
            error = abs((phase-expected+0.5) % 1-0.5)
            worst = max(worst, error/tolerance)
            if read_x != x or not clean or error > tolerance:
                failures += 1
                if failures <= 3:
                    print(f"FAIL x={x} shot={shot} x_out={read_x} "
                          f"clean={clean} error_ratio={error/tolerance:.4g}")
    stats = circuit.operation_statistics()
    return dict(inputs_checked=max_x, shots_per_x=shots_per_x, failures=failures,
                worst_error_ratio=worst, qubits=circuit.num_qubits,
                toffoli_like=stats.num_toffoli_like, t_like=stats.num_t_like,
                z_pow=stats.num_z_pow())
