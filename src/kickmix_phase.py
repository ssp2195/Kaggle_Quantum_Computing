"""Kickmix phase implementation using prefix matching and parity gadgets.

Experimental only. All scratch is uncomputed, but no official verification
or cost claim is made until the native Kaggle verifier passes.
"""
from itertools import combinations
from math import sqrt
from .direct_phase import Block, INPUT_BITS, TOL


def monomial_angles(block: Block, prune_fraction: float = 0.05):
    """Expand a local polynomial into ANDs of low input bits.

    Each (mask, turns) means turns * AND(bits in mask), conditioned
    on the block's fixed high-bit prefix. mask=0 is constant.
    """
    import numpy as np
    k = block.size.bit_length()-1
    def p(offset):
        return float(block.approximation(np.array([block.start+offset]))[0])
    angles = []
    for r in range(min(k, block.degree)+1):
        for pos in combinations(range(k), r):
            coeff = 0.0
            for t in range(r+1):
                for selected in combinations(pos, t):
                    coeff += (-1)**(r-t) * p(sum(1 << i for i in selected))
            angles.append((sum(1 << i for i in pos), coeff))
    # The sum of absolute discarded coefficients bounds added phase error.
    budget = prune_fraction*TOL/sqrt(block.end)
    dropped = 0.0
    kept = list(angles)
    for mask, coeff in sorted(angles, key=lambda a: abs(a[1])):
        if dropped + abs(coeff) <= budget:
            kept.remove((mask, coeff))
            dropped += abs(coeff)
    return kept, dropped


def estimated_rotations(blocks, prune_fraction=0.05):
    """Non-fused gadget Z_POW count, not the official synthesis cost."""
    total = 0
    for b in blocks:
        terms, _ = monomial_angles(b, prune_fraction)
        total += sum((1 << (mask.bit_count()+1))-1 for mask, _ in terms)
    return total


def _and_phase_via_parities(builder, wires, phase_turns):
    """Apply exp(2 pi i * phase_turns * AND(wires)).

    Boolean AND is a signed linear combination of nonempty XOR parities.
    The kickmix Z_POW exponent is in HALF-turns, hence factor 2.
    """
    m = len(wires)
    for r in range(1, m+1):
        for ids in combinations(range(m), r):
            subset = [wires[i] for i in ids]
            target = subset[-1]
            for q in subset[:-1]:
                builder.cx(q, target)
            exponent = 2*phase_turns*((-1)**(r-1))/(1 << (m-1))
            builder.z_pow(target, format(exponent, ".17g"))
            for q in reversed(subset[:-1]):
                builder.cx(q, target)


def build_circuit(blocks: list[Block], input_bits: int = INPUT_BITS,
                  prune_fraction: float = 0.05):
    """Compile a direct-phase candidate into a native kickmix Circuit."""
    import kickmix as km
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(input_bits, name="x")
    flag = builder.alloc_qubits(1)[0]
    for block in blocks:
        k = block.size.bit_length()-1
        prefix = block.start >> k
        high, low = x[k:], x[:k]
        negative = [high[i] for i in range(len(high)) if not (prefix >> i) & 1]
        for q in negative:
            builder.x(q)
        builder.multi_controlled_x(high, flag)
        for q in negative:
            builder.x(q)
        for mask, turns in monomial_angles(block, prune_fraction)[0]:
            wires = [flag] + [low[j] for j in range(k) if mask & (1 << j)]
            _and_phase_via_parities(builder, wires, turns)
        for q in negative:
            builder.x(q)
        builder.multi_controlled_x(high, flag)
        for q in negative:
            builder.x(q)
    builder.reset(flag)
    builder.free(flag)
    return builder.finish_circuit()
