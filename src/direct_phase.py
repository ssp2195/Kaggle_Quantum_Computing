"""Numerical dyadic direct-phase approximations for Benchmark I."""
from dataclasses import dataclass
from math import sqrt
import numpy as np
from numpy.polynomial import Chebyshev

MAX_X = 3_133_452
INPUT_BITS = 22
TOL = 2.0 ** -24


@dataclass
class Block:
    start: int
    size: int
    degree: int
    coeffs: tuple[float, ...]
    max_ratio: float

    @property
    def end(self):
        return self.start + self.size - 1

    def approximation(self, x):
        if self.size == 1:
            return np.full_like(np.asarray(x, dtype=float), self.coeffs[0])
        return Chebyshev(self.coeffs, domain=[self.start, self.end])(x)


def fit_block(start: int, size: int, degree: int) -> Block:
    end = start + size - 1
    if size == 1:
        return Block(start, size, degree, (1 / sqrt(start),), 0.0)
    poly = Chebyshev.interpolate(lambda x: 1.0 / np.sqrt(x), degree,
                                  domain=[start, end])
    xx = np.arange(start, end + 1, dtype=np.float64)
    ratio = np.max(np.abs(poly(xx) - 1.0 / np.sqrt(xx)) * np.sqrt(xx) / TOL)
    return Block(start, size, degree, tuple(float(a) for a in poly.coef),
                 float(ratio))


def dyadic_partition(degree: int, max_x: int = MAX_X, safety: float = 0.75,
                     max_k: int = 18) -> list[Block]:
    """Greedily take the longest aligned interval passing EVERY integer."""
    if degree not in (1, 2, 3) or not 0 < safety < 1:
        raise ValueError("degree must be 1..3 and safety must be in (0,1)")
    blocks = []
    start = 1
    while start <= max_x:
        allowed_k = min(max_k, (start & -start).bit_length()-1,
                        (max_x-start+1).bit_length()-1)
        for k in range(allowed_k, -1, -1):
            block = fit_block(start, 1 << k, degree)
            if block.max_ratio <= safety:
                break
        blocks.append(block)
        start += block.size
    return blocks


def verify_partition(blocks: list[Block], max_x: int = MAX_X) -> dict:
    expected = 1
    worst = 0.0
    worst_x = 1
    for block in blocks:
        if block.start != expected:
            raise ValueError(f"gap/overlap at {expected}")
        xx = np.arange(block.start, block.end+1, dtype=np.float64)
        ratio = np.abs(block.approximation(xx) - 1 / np.sqrt(xx)) * np.sqrt(xx) / TOL
        idx = int(np.argmax(ratio))
        if ratio[idx] > worst:
            worst, worst_x = float(ratio[idx]), int(xx[idx])
        expected = block.end+1
    if expected != max_x+1:
        raise ValueError("partition does not cover domain")
    return {"blocks": len(blocks), "max_x": max_x, "worst_x": worst_x,
            "max_tolerance_fraction": worst, "passes_numerical_tolerance": worst < 1}
