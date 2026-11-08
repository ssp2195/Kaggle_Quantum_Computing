"""Quantized segmented quadratic phase model for Benchmark I."""
from __future__ import annotations
from dataclasses import dataclass
from math import sqrt
import numpy as np

MAX_X = 3_133_452
TOL = 2**-24
NORMAL_BITS = 22
MANTISSA_BITS = 20

@dataclass(frozen=True)
class Config:
    segment_bits: int = 8
    precision: int = 28

    def __post_init__(self):
        if not 4 <= self.segment_bits <= 12 or not 22 <= self.precision <= 36:
            raise ValueError("segment_bits in 4..12; precision in 22..36")

    @property
    def low_bits(self):
        return MANTISSA_BITS - self.segment_bits


def coefficient_table(config: Config):
    """QROM rows contain A, positive B, positive C for A - Bu + Cu^2."""
    p, q = config.segment_bits, config.precision
    rows = []
    for j in range(1 << (p + 2)):
        if j < (1 << p):
            rows.append((0, 0, 0))
            continue
        z = j / (1 << p)
        a = round(1 / sqrt(z) * (1 << q))
        b = round(.5 / (1 << p) / z**1.5 * (1 << q))
        c = round(.375 / (1 << (2*p)) / z**2.5 * (1 << q))
        rows.append((a, b, c))
    widths = tuple(max(1, max(row[k] for row in rows).bit_length()) for k in range(3))
    packed = [a | (b << widths[0]) | (c << (widths[0] + widths[1]))
              for a, b, c in rows]
    return packed, widths


def phase_model(config: Config, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.int64)
    if np.any(x < 1) or np.any(x >= (1 << NORMAL_BITS)):
        raise ValueError("Expected positive 22-bit inputs")
    k = np.floor(np.log2(x)).astype(np.int64) // 2
    norm = x << (2 * (10 - k))
    lo = config.low_bits
    j = norm >> lo
    r = (norm & ((1 << lo) - 1)).astype(np.float64) / (1 << lo)
    z = j.astype(np.float64) / (1 << config.segment_bits)
    factor = float(1 << config.precision)
    a = np.rint((1 / np.sqrt(z)) * factor)
    b = np.rint((.5 / (1 << config.segment_bits) / z**1.5) * factor)
    c = np.rint((.375 / (1 << (2 * config.segment_bits)) / z**2.5) * factor)
    return np.ldexp((a - b * r + c * r * r) / factor, -k)


def full_domain_screen(config: Config, max_x: int = MAX_X, chunk: int = 200_000):
    worst = -1.0
    worst_x = 0
    failed = 0
    for left in range(1, max_x + 1, chunk):
        x = np.arange(left, min(max_x + 1, left + chunk), dtype=np.int64)
        target = 1 / np.sqrt(x.astype(float))
        ratio = np.abs(phase_model(config, x) - target) / target / TOL
        idx = int(np.argmax(ratio))
        if ratio[idx] > worst:
            worst, worst_x = float(ratio[idx]), int(x[idx])
        failed += int(np.count_nonzero(ratio > 1))
    _, widths = coefficient_table(config)
    return dict(segment_bits=config.segment_bits, precision=config.precision,
                interval_count=3 * (1 << config.segment_bits),
                table_rows=1 << (config.segment_bits + 2),
                coefficient_widths=widths, worst_ratio=worst, worst_x=worst_x,
                failed_inputs=failed, passed=failed == 0)
