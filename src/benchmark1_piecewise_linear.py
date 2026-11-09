"""Piecewise-linear *direct phase* architecture for Benchmark I.

Eliminates the quadratic multiplier and square register entirely, at the
cost of a finer mantissa interval QROM. All arithmetic remains reversible.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .hybrid_model import MAX_X, TOL, NORMAL_BITS, MANTISSA_BITS
from .hybrid_kickmix import _phase, resource_counts
from .priority1_normalize import normalize_barrel, denormalize_barrel
from .priority1_truncated import _mul_cut
from .priority1_exact import truncated_product


@dataclass(frozen=True)
class LinearConfig:
    segment_bits: int = 10
    precision: int = 28
    linear_drop: int = 7

    @property
    def low_bits(self):
        return MANTISSA_BITS - self.segment_bits

    def __post_init__(self):
        if not 8 <= self.segment_bits <= 12:
            raise ValueError("segment bits must be 8..12")
        if not 24 <= self.precision <= 32:
            raise ValueError("fractional precision must be 24..32")
        if not 0 <= self.linear_drop <= self.low_bits:
            raise ValueError("invalid linear partial-product cutoff")


@lru_cache(maxsize=16)
def coefficients(cfg: LinearConfig):
    """[A_j, B_j] where phase is 2**(-k)*(A_j-B_j*u)/2**q."""
    p, q = cfg.segment_bits, cfg.precision
    mat = np.zeros((4 << p, 2), dtype=np.int64)
    nodes = .5 + .5 * np.cos(np.pi * (np.arange(2) + .5) / 2)
    vander = np.polynomial.polynomial.polyvander(nodes, 1)
    for j in range(1 << p, 4 << p):
        z0 = j / (1 << p)
        f = 1.0 / np.sqrt(z0 + nodes / (1 << p))
        polynomial = np.linalg.solve(vander, f)
        mat[j] = np.rint(polynomial * (1 << q)).astype(np.int64) * [1, -1]
    if np.any(mat < 0):
        raise ArithmeticError("Unexpected signed interpolation coefficient")
    mat.setflags(write=False)
    return mat


def screen(cfg: LinearConfig, max_x=MAX_X):
    """Exhaustive finite-domain phase screen with exact partial-product omissions."""
    mat = coefficients(cfg)
    lo, q, drop = cfg.low_bits, cfg.precision, cfg.linear_drop
    worst = 0.0
    worst_x = 1
    failed = 0
    for start in range(1, max_x + 1, 200_000):
        x = np.arange(start, min(start + 200_000, max_x + 1), dtype=np.int64)
        exponent = np.floor(np.log2(x)).astype(np.int64) // 2
        normalized = x << (20 - 2 * exponent)
        idx, offset = normalized >> lo, normalized & ((1 << lo) - 1)
        a, b = mat[idx, 0], mat[idx, 1]
        linear = truncated_product(b, offset, drop)
        phase = np.ldexp(a / (1 << q) -
                         linear / (1 << (q + lo - drop)), -exponent)
        target = 1.0 / np.sqrt(x)
        ratio = np.abs(phase / target - 1) / TOL
        i = int(np.argmax(ratio))
        if ratio[i] > worst:
            worst, worst_x = float(ratio[i]), int(x[i])
        failed += int(np.count_nonzero(ratio > 1))
    return dict(family="piecewise-linear-Chebyshev", segment_bits=cfg.segment_bits,
                precision=cfg.precision, linear_drop=drop,
                intervals=3 * (1 << cfg.segment_bits),
                lookup_rows=4 << cfg.segment_bits,
                worst_error_ratio=worst, worst_x=worst_x,
                failed_inputs=failed, checked_inputs=max_x, passes=failed == 0)


def build(cfg: LinearConfig, grouping="joint"):
    """Construct native kickmix circuit; joint or sequential coefficient QROM."""
    if grouping not in ("joint", "streamed"):
        raise ValueError("grouping must be joint or streamed")
    result = screen(cfg)
    if not result["passes"] or result["worst_error_ratio"] >= .9:
        raise ValueError(f"Mathematical phase candidate insufficiently accurate: {result}")
    import kickmix as km
    mat = coefficients(cfg)
    widths = [max(1, int(mat[:, i].max()).bit_length()) for i in range(2)]
    lo = cfg.low_bits
    builder = km.CircuitBuilder()
    x = builder.create_quantum_register(NORMAL_BITS, name="x")
    flags, shift = normalize_barrel(builder, x)

    def load(indices):
        word_width = sum(widths[i] for i in indices)
        table = []
        for row in mat:
            val, pos = 0, 0
            for i in indices:
                val |= int(row[i]) << pos
                pos += widths[i]
            table.append(val)
        buf = builder.alloc_qubits(word_width)
        builder.init_lookup(table=table, address=x[lo:], target=buf)
        values = []
        pos = 0
        for i in indices:
            values.append(buf[pos:pos + widths[i]])
            pos += widths[i]
        return buf, table, values

    def release(buf, table):
        # Reapply XOR QROM to erase coefficients coherently.
        # Measurement-based del_lookup needs additional unary-cleanup workspace
        # for the large 12-bit address space (often unavailable at peak Q).
        builder.init_lookup(table=table, address=x[lo:], target=buf)
        builder.free(buf)

    def linear_phase(b):
        product = builder.alloc_qubits(len(b) + lo - cfg.linear_drop)
        builder.reset(product)
        _mul_cut(builder, b, x[:lo], product, cfg.linear_drop)
        _phase(builder, product, shift, cfg.precision + lo - cfg.linear_drop, -1)
        _mul_cut(builder, b, x[:lo], product, cfg.linear_drop, inverse=True)
        builder.free(product)

    if grouping == "joint":
        with builder.mark("joint-coefficients"):
            buf, table, (a, b) = load((0, 1))
            _phase(builder, a, shift, cfg.precision, +1)
            linear_phase(b)
            release(buf, table)
    else:
        with builder.mark("constant"):
            buf, table, (a,) = load((0,))
            _phase(builder, a, shift, cfg.precision, +1)
            release(buf, table)
        with builder.mark("linear"):
            buf, table, (b,) = load((1,))
            linear_phase(b)
            release(buf, table)
    denormalize_barrel(builder, x, flags, shift)
    return builder.finish_circuit()
