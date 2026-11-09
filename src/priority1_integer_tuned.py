"""Quantized coefficients optimized for the *actual truncated* reversible circuit.

This is different from rounding a real-valued minimax polynomial: each trial
coefficient is evaluated with the same omitted partial products as _mul_cut.
The deterministic coordinate search is restricted to degree=2, p=7, q=25,
drops=(9,11,11). It is a candidate generator, NOT a proof of optimality.
"""
from functools import lru_cache
import numpy as np

DROPS = (9, 11, 11)
SEGMENT_BITS = 7
PRECISION = 25
LOW_BITS = 13


def truncated_product(a, b, drop):
    """Exactly match omitted partial products, including absent low carries."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    omitted = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.int64)
    for bit in range(drop):
        omitted += ((b >> bit) & 1) * (a & ((1 << (drop - bit)) - 1)) * (1 << bit)
    return (a * b - omitted) >> drop


@lru_cache(maxsize=1)
def tuned_matrix():
    """Return 4*2**p rows of nonnegative (A,B,C) integer coefficients."""
    from numpy.polynomial import Chebyshev, Polynomial
    p, q, low = SEGMENT_BITS, PRECISION, LOW_BITS
    dl, ds, dq = DROPS
    offsets = np.arange(1 << low, dtype=np.int64)
    u = offsets.astype(float) / (1 << low)
    square = truncated_product(offsets, offsets, ds)
    result = np.zeros((4 << p, 3), dtype=np.int64)
    for index in range(1 << p, 4 << p):
        z0 = index / (1 << p)
        z = z0 + u / (1 << p)
        target = 1 / np.sqrt(z)
        poly = Chebyshev.interpolate(lambda v: 1 / np.sqrt(z0 + v / (1 << p)),
                                     2, domain=[0, 1]).convert(kind=Polynomial)
        initial = np.rint(poly.coef * (1 << q)).astype(np.int64)
        c = np.array([initial[0], -initial[1], initial[2]], dtype=np.int64)

        def error(v):
            lin = truncated_product(v[1], offsets, dl)
            quad = truncated_product(v[2], square, dq)
            approximation = (v[0] / (1 << q)
                             - lin / (1 << (q + low - dl))
                             + quad / (1 << (q + 2 * low - ds - dq)))
            return float(np.max(np.abs(approximation / target - 1)) * (1 << 24))

        for _ in range(8):
            changed = False
            for power in (0, 1, 2, 1, 0):
                candidates = []
                for delta in range(-4, 5):
                    trial = c.copy()
                    trial[power] += delta
                    candidates.append(error(trial))
                best = int(np.argmin(candidates)) - 4
                if best:
                    c[power] += best
                    changed = True
            if not changed:
                break
        if np.any(c < 0):
            raise ValueError(f"Negative coefficient at segment {index}")
        result[index] = c
    result.setflags(write=False)
    return result
