"""Piecewise phase approximations for Benchmark I; no circuit-level shortcuts.

Each normalized input z=(j+u)/2**segments uses f(z)=1/sqrt(z).
Coefficients are integer-quantized monomial coefficients in u in [0,1).
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from math import cos, pi, sqrt
import numpy as np
from .hybrid_model import MAX_X, NORMAL_BITS, MANTISSA_BITS, TOL


@dataclass(frozen=True)
class Approximation:
    family: str = "chebyshev"
    degree: int = 2
    segment_bits: int = 6
    precision: int = 26

    def __post_init__(self):
        if self.family not in ("taylor", "chebyshev", "minimax", "chebyshev-tuned"):
            raise ValueError("family must be taylor, chebyshev, minimax, or chebyshev-tuned")
        if self.degree not in (2, 3) or self.segment_bits not in range(4, 11):
            raise ValueError("degree 2 or 3, segment_bits 4..10 required")
        if self.family == "chebyshev-tuned" and (self.degree,self.segment_bits,self.precision) != (2,6,25):
            raise ValueError("The quantization-aware adjustment is certified only for 2/6/25")
        if not 24 <= self.precision <= 34:
            raise ValueError("precision 24..34 required")

    @property
    def low_bits(self):
        return MANTISSA_BITS - self.segment_bits


def _coefficients(family: str, degree: int, z0: float, h: float):
    if family == "chebyshev-tuned":
        family = "chebyshev"
    if family == "taylor":
        # binomial series of (1 + h*u/z0)^(-1/2).
        a = [1.0]
        for n in range(1, degree+1):
            a.append(a[-1] * (-0.5-(n-1)) / n * h/z0)
        return np.array(a)*z0**(-.5)

    nodes = 0.5 + 0.5*np.cos(
        np.pi*(np.arange(degree+1)+0.5)/(degree+1))
    if family == "chebyshev":
        return np.linalg.solve(
            np.polynomial.polynomial.polyvander(nodes, degree),
            1 / np.sqrt(z0+h*nodes))
    # Weighted minimax linear program on a dense normalized interpolation mesh.
    # The subsequent integer quantization can worsen the minimax solution;
    # full-domain validation is mandatory and failure is never suppressed.
    from scipy.optimize import linprog
    u = np.r_[np.linspace(0, 1, 65),
              0.5+0.5*np.cos(np.pi*np.arange(45)/44)]
    y = 1/np.sqrt(z0+h*u)
    w = np.sqrt(z0+h*u)
    v = np.polynomial.polynomial.polyvander(u, degree)*w[:,None]
    aub = np.concatenate([
        np.column_stack([v, -np.ones(len(u))]),
        np.column_stack([-v, -np.ones(len(u))])])
    bub = np.concatenate([w*y, -w*y])
    bounds = [(None,None)]*(degree+1)+[(0,None)]
    res = linprog(np.r_[np.zeros(degree+1),1.],
                  A_ub=aub,b_ub=bub,bounds=bounds,method="highs")
    if not res.success:
        raise ArithmeticError(res.message)
    return res.x[:-1]


@lru_cache(maxsize=32)
def coefficient_matrix(config: Approximation) -> np.ndarray:
    """Nonnegative packed magnitude coefficients; signs alternate + - + -."""
    p, q, d = config.segment_bits, config.precision, config.degree
    mat = np.zeros((4 << p, d+1), dtype=np.int64)
    h = 2.0**(-p)
    for j in range(1 << p, 4 << p):
        signed = np.rint(_coefficients(config.family, d, j*h, h) * 2**q)
        # Alternating signs are a property of these smooth approximations.
        # Enforce this, not abs(...), to avoid silently changing the polynomial.
        for i, a in enumerate(signed):
            if ((-1)**i) * a < 0:
                raise ArithmeticError(
                    f"Unexpected sign at interval {j} coefficient {i}: {a}")
            mat[j,i] = int(((-1)**i)*a)
    if config.family == "chebyshev-tuned":
        # Integer coefficient adjustments found by deterministic coordinate
        # search on all 2**14 normalized offsets, then exhaustively checked
        # on all 3,133,452 benchmark inputs. Entries are (interval, power):
        # signed correction to the nonnegative coefficient magnitude.
        improvements = {
            (78,1): +1, (172,1): +1, (203,0): +1,
            (224,1): -1, (231,1): -1, (246,1): +1,
            (255,1): -1,
        }
        for (interval,power),delta in improvements.items():
            mat[interval,power] += delta
    mat.setflags(write=False)
    return mat


def pack_coefficients(config: Approximation):
    mat = coefficient_matrix(config)
    widths = tuple(max(1,int(mat[:,i].max()).bit_length())
                   for i in range(config.degree+1))
    table = []
    for row in mat:
        val, shift = 0, 0
        for v, bits in zip(row,widths):
            val |= int(v) << shift
            shift += bits
        table.append(val)
    return table, widths


def full_domain_screen(config: Approximation, max_x: int = MAX_X):
    """Exhaustive per-integer check of exact quantized polynomial model."""
    mat = coefficient_matrix(config)
    lo, p, q, d = (config.low_bits,config.segment_bits,
                    config.precision,config.degree)
    worst, worst_x, failures = 0., 1, 0
    for begin in range(1,max_x+1,150_000):
        x = np.arange(begin,min(begin+150_000,max_x+1),dtype=np.int64)
        k = np.floor(np.log2(x)).astype(np.int64)//2
        zint = x << (20 - 2*k)
        j, low = zint >> lo, zint & ((1 << lo)-1)
        u = low.astype(float)/(1 << lo)
        c = mat[j]
        phase = np.zeros_like(u)
        for i in reversed(range(d+1)):
            phase = phase*u + ((-1)**i)*c[:,i]
        phase = np.ldexp(phase/(1 << q),-k)
        exact = 1 / np.sqrt(x.astype(float))
        ratio = np.abs(phase-exact) / (TOL*exact)
        m = int(np.argmax(ratio))
        if ratio[m] > worst:
            worst, worst_x = float(ratio[m]),int(x[m])
        failures += int(np.count_nonzero(ratio>1))
    _,widths = pack_coefficients(config)
    return dict(family=config.family,degree=d,segment_bits=p,precision=q,
                intervals=3*(1<<p),coefficient_widths=widths,
                worst_error_ratio=worst,worst_x=worst_x,
                failures=failures,passes=failures==0)
