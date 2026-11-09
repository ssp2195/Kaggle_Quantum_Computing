"""Full 20-bit normalized-mantissa screen of a sparse Boolean ANF phase.

This is a mathematical rejection screen of naive independent monomial phase
gadgets, not a proof against shared-predicate synthesis or a quantum circuit.
"""
import argparse
import json
import numpy as np


def mobius(values, inverse=False):
    """Boolean-lattice Mobius transform (inverse=True applies zeta sum)."""
    out = np.array(values, dtype=np.float64, copy=True)
    n = out.size
    if n < 1 or n & (n-1):
        raise ValueError("size must be a power of two")
    for bit in range(n.bit_length()-1):
        step = 1 << bit
        blocks = out.reshape(-1, 2*step)
        if inverse:
            blocks[:, step:] += blocks[:, :step]
        else:
            blocks[:, step:] -= blocks[:, :step]
    return out


def screen(bits=20, counts=(1000,1500,2000,3000)):
    if not 4 <= bits <= 21:
        raise ValueError("bits must lie in 4..21")
    n = 1 << bits
    target = 1 / np.sqrt(1 + np.arange(n, dtype=np.float64)/n)
    coefficient = mobius(target)
    ranked = np.argsort(np.abs(coefficient))[::-1]
    result = []
    for count in counts:
        if count < 1 or count > n:
            raise ValueError("bad count")
        retained = np.zeros(n, dtype=np.float64)
        retained[ranked[:count]] = coefficient[ranked[:count]]
        approximate = mobius(retained, inverse=True)
        ratio = np.abs(approximate/target-1) * (1 << 24)
        x = int(np.argmax(ratio))
        result.append(dict(terms=int(count), domain_points=n,
                           worst_error_ratio=float(ratio[x]), worst_index=x,
                           failed_points=int(np.count_nonzero(ratio>1)),
                           passes_numerical=bool(ratio[x]<1)))
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--bits", type=int, default=20)
    p.add_argument("--output")
    a = p.parse_args()
    result = screen(a.bits)
    print(json.dumps(result, indent=2))
    if a.output:
        from pathlib import Path
        dest = Path(a.output)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(result, indent=2)+"\n")
