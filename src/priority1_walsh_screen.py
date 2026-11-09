"""Full-grid Walsh parity-phase prescreen for normalized 1/sqrt(z).

Investigates direct diagonal-phase synthesis WITHOUT reversible polynomial
multipliers. Walsh basis: (-1)**(bitwise_dot(mask, normalized_input)).
The scalar constant phase is retained in numerical screening but may be
ignored as a global phase for a single fixed-exponent block.

This screen is an *architecture feasibility check*, not a Kaggle circuit
verifier or a proof of lower bounds for all possible diagonal compilers.
"""
from __future__ import annotations
import argparse
import json
import numpy as np


def walsh_hadamard(a):
    """Unnormalized FWHT. Input length must be a power of two."""
    v = np.array(a, dtype=np.float64, copy=True)
    n = v.size
    if n <= 0 or n & (n - 1):
        raise ValueError("Expected positive power-of-two length")
    step = 1
    while step < n:
        blocks = v.reshape(-1, 2 * step)
        left = blocks[:, :step].copy()
        right = blocks[:, step:].copy()
        blocks[:, :step] = left + right
        blocks[:, step:] = left - right
        step *= 2
    return v


def screen(bits: int = 20, terms=(500, 1000, 1500, 1700, 1800, 2000, 3000)):
    """Check every normalized value 1 + m*2^-bits, 0<=m<2^bits."""
    if not 4 <= bits <= 21:
        raise ValueError("bits must be between 4 and 21")
    n = 1 << bits
    mantissa = np.arange(n, dtype=np.float64) / n
    z = 1.0 + mantissa
    target = 1 / np.sqrt(z)
    coeff = walsh_hadamard(target) / n
    order = np.argsort(np.abs(coeff))[::-1]
    rows = []
    for count in terms:
        if not 0 < count <= n:
            raise ValueError("Term count outside valid range")
        subset = np.zeros(n, dtype=np.float64)
        subset[order[:count]] = coeff[order[:count]]
        approx = walsh_hadamard(subset)
        error = np.abs(approx / target - 1) * (1 << 24)
        index = int(np.argmax(error))
        rows.append(dict(terms=int(count), normalized_values=n,
                         worst_error_over_tolerance=float(error[index]),
                         worst_mantissa=index,
                         numerical_pass=bool(error[index] < 1)))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bits", type=int, default=20)
    p.add_argument("--terms", type=int, nargs="+",
                   default=[500, 1000, 1500, 1700, 1800, 2000, 3000])
    p.add_argument("--output")
    a = p.parse_args()
    data = screen(a.bits, a.terms)
    print(json.dumps(data, indent=2))
    if a.output:
        from pathlib import Path
        dest = Path(a.output)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
