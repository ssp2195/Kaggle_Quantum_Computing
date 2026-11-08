"""Native kickmix normalized quadratic *phase* without Newton refinement."""
from __future__ import annotations
from fractions import Fraction
from .hybrid_model import Config, coefficient_table, NORMAL_BITS


def _normalize(builder, x):
    flags = []
    for _ in range(10):
        builder.x(x[20]); builder.x(x[21])
        flag = builder.alloc_qubits(1)[0]
        builder.init_and(x[20], x[21], flag)
        builder.x(x[20]); builder.x(x[21])
        builder.left_rotate(x, shift=2, control=flag)
        flags.append(flag)
    shift = builder.alloc_qubits(4)
    builder.reset(shift)
    for flag in flags:
        builder.iadd(1, target=shift, control=flag)
    return flags, shift


def _denormalize(builder, x, flags, shift):
    for flag in reversed(flags):
        builder.isub(1, target=shift, control=flag)
    builder.free(shift)
    for flag in reversed(flags):
        builder.right_rotate(x, shift=2, control=flag)
        builder.x(x[20]); builder.x(x[21])
        builder.del_and(x[20], x[21], flag)
        builder.x(x[20]); builder.x(x[21])
        builder.free(flag)


def _phase(builder, word, shift, frac_bits, sign):
    """Apply sign*(word/2**frac_bits)*2**(shift-10) turns.

    The 10 high padding qubits stay zero before and after variable rotation.
    Z_POW uses a half-turn exponent; exact fractions avoid float rounding.
    """
    pad = builder.alloc_qubits(10)
    builder.reset(pad)
    extended = word + pad
    builder.left_rotate(extended, shift=shift)
    denom = frac_bits + 10
    for bit in range(len(extended)):
        exp = Fraction(sign * (1 << (bit + 1)), 1 << denom)
        if exp % 2:
            builder.z_pow(extended[bit], exp)
    builder.right_rotate(extended, shift=shift)
    builder.free(pad)


def _mul(builder, a, b, dest, inverse=False):
    """Reversible unsigned product into dest using carry-safe row additions."""
    na, nb = len(a), len(b)
    if len(dest) != na + nb:
        raise ValueError("Product width must equal len(a)+len(b)")
    for j in (reversed(range(nb)) if inverse else range(nb)):
        ctl = b[j]
        row = builder.alloc_qubits(na + 1)
        builder.reset(row)
        for i in range(na):
            if a[i] == ctl:
                builder.cx(ctl, row[i])
            else:
                builder.init_and(a[i], ctl, row[i])
        if inverse:
            builder.isub(row, target=dest[j:j + na + 1])
        else:
            builder.iadd(row, target=dest[j:j + na + 1])
        for i in reversed(range(na)):
            if a[i] == ctl:
                builder.cx(ctl, row[i])
            else:
                builder.del_and(a[i], ctl, row[i])
        builder.free(row)


def build_hybrid(config: Config):
    import kickmix as km
    packed, widths = coefficient_table(config)
    wa, wb, wc = widths
    lo = config.low_bits
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(NORMAL_BITS, name="x")
    with builder.mark("normalize"):
        flags, shift = _normalize(builder, x)
    with builder.mark("coefficients"):
        coeff = builder.alloc_qubits(wa + wb + wc)
        builder.init_lookup(table=packed, address=x[lo:], target=coeff)
        a, b, c = coeff[:wa], coeff[wa:wa + wb], coeff[wa + wb:]

    with builder.mark("intercept_phase"):
        _phase(builder, a, shift, config.precision, +1)

    with builder.mark("linear_phase"):
        prod = builder.alloc_qubits(wb + lo)
        builder.reset(prod)
        _mul(builder, b, x[:lo], prod)
        _phase(builder, prod, shift, config.precision + lo, -1)
        _mul(builder, b, x[:lo], prod, inverse=True)
        builder.free(prod)

    with builder.mark("quadratic_phase"):
        square = builder.alloc_qubits(2 * lo)
        builder.reset(square)
        _mul(builder, x[:lo], x[:lo], square)
        prod = builder.alloc_qubits(wc + 2 * lo)
        builder.reset(prod)
        _mul(builder, c, square, prod)
        _phase(builder, prod, shift, config.precision + 2 * lo, +1)
        _mul(builder, c, square, prod, inverse=True)
        builder.free(prod)
        _mul(builder, x[:lo], x[:lo], square, inverse=True)
        builder.free(square)

    with builder.mark("cleanup"):
        builder.del_lookup(table=packed, address=x[lo:], target=coeff)
        builder.free(coeff)
        _denormalize(builder, x, flags, shift)
    return builder.finish_circuit()


def resource_counts(circuit):
    stats = circuit.operation_statistics()
    return dict(qubits=circuit.num_qubits, instructions=stats.num_instructions,
                toffoli_like=stats.num_toffoli_like, t_like=stats.num_t_like,
                fine_rotations=stats.num_z_pow(),
                lower_bound_qc=circuit.num_qubits *
                (stats.num_toffoli_like + stats.num_t_like / 4))
