"""Low-workspace, in-place Horner phase-oracle candidate (actual kickmix).

Normalized quadratic phase P(u) = A - B*u + C*u*u, u=low/2**13.
We load A<<g and B<<g in the existing coefficient table, then:
  R <- B<<g - CUT(C*low, low_bits-g)
  Y <- A<<g - CUT(R*low, low_bits)
  Z_POW(Y scaled by 2**(-q-g-k))
  undo Y and R, then erase the coefficient lookup, normalize x back.

CUT is exactly the sum of selected bit partial-products, not a rounded
ordinary product. Direct updates eliminate independent full product/square
registers and share row temporaries. MBUC del_and or coherent CCX cleanup
can be compared without changing the mathematical function.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .priority1_model import Approximation, coefficient_matrix
from .priority1_exact import truncated_product
from .priority1_normalize import normalize_barrel, denormalize_barrel
from .benchmark1_compact_normalize import normalize_compact, denormalize_compact
from .hybrid_model import MAX_X, TOL, NORMAL_BITS
from .hybrid_kickmix import _phase


@dataclass(frozen=True)
class HornerConfig:
    guard: int = 4
    cleanup: str = "mbuc"
    normalization: str = "barrel"

    def __post_init__(self):
        if self.guard not in (2, 3, 4, 5, 6):
            raise ValueError("guard bits must be between 2 and 6")
        if self.cleanup not in ("mbuc", "coherent"):
            raise ValueError("cleanup must be mbuc or coherent")
        if self.normalization not in ("barrel", "serial", "compact"):
            raise ValueError("normalization must be barrel, serial or compact")


COEFF_CONFIG = Approximation("chebyshev", 2, 7, 26)


def screen(config: HornerConfig, max_x=MAX_X):
    """All allowed x, exact omitted-partial-product integer arithmetic."""
    mat = coefficient_matrix(COEFF_CONFIG)
    lo, q, g = COEFF_CONFIG.low_bits, COEFF_CONFIG.precision, config.guard
    worst, worst_x, failed = 0., 1, 0
    for begin in range(1, max_x+1, 160_000):
        x = np.arange(begin, min(max_x+1, begin+160_000), dtype=np.int64)
        k = np.floor(np.log2(x)).astype(np.int64)//2
        z = x << (20-2*k)
        idx = z >> lo
        low = z & ((1 << lo)-1)
        a, beta, gamma = (mat[idx,i] for i in range(3))
        r = (beta << g) - truncated_product(gamma, low, lo-g)
        y = (a << g) - truncated_product(r, low, lo)
        turns = np.ldexp(y.astype(float)/(1 << (q+g)), -k)
        exact = 1/np.sqrt(x.astype(float))
        ratio = np.abs(turns / exact - 1) / TOL
        j = int(np.argmax(ratio))
        if ratio[j] > worst:
            worst, worst_x = float(ratio[j]), int(x[j])
        failed += int(np.count_nonzero(ratio > 1))
    return dict(architecture="horner-in-place", guard=config.guard,
                cleanup=config.cleanup, normalization=config.normalization,
                worst_error_ratio=worst, worst_x=worst_x,
                failed_inputs=failed, checked_inputs=max_x,
                passes=failed==0)


def _partial_product_update(builder, source, multiplier, dest, drop, *, subtract, cleanup):
    """Update dest in place by selected positive bit partial products.

    The destination must not alias source or multiplier. The max high
    column is included in dest; final extra row bit handles ripple carry.
    For undo, call again with opposite subtract and reverse=True.\n    Critically update dest[start:] instead of a row-sized slice: a borrow\n    from a nonzero coefficient must propagate through all higher bits.
    """
    n, m = len(source), len(multiplier)
    if not 0 <= drop < n+m:
        raise ValueError("bad cut")
    for j in range(m):
        first = max(0, drop-j)
        if first >= n:
            continue
        start = j+first-drop
        length = n-first+1
        if start + length > len(dest):
            raise ValueError("Destination too narrow for partial products")
        row = builder.alloc_qubits(length)
        builder.reset(row)
        for i in range(first,n):
            if source[i] == multiplier[j]:
                builder.cx(source[i], row[i-first])
            else:
                builder.init_and(source[i], multiplier[j], row[i-first])
        if subtract:
            builder.isub(row, target=dest[start:])
        else:
            builder.iadd(row, target=dest[start:])
        for i in reversed(range(first,n)):
            if source[i] == multiplier[j]:
                builder.cx(source[i], row[i-first])
            elif cleanup == "mbuc":
                builder.del_and(source[i], multiplier[j], row[i-first])
            else:
                builder.ccx(source[i], multiplier[j], row[i-first])
        builder.free(row)


def _inplace_product(builder, source, multiplier, dest, drop, *,
                     subtract=True, cleanup="mbuc", inverse=False):
    """Forward or inverse row updates in exact reverse order."""
    n,m = len(source),len(multiplier)
    js = range(m-1,-1,-1) if inverse else range(m)
    for j in js:
        first=max(0,drop-j)
        if first >= n: continue
        start=j+first-drop
        length=n-first+1
        if start+length>len(dest):
            raise ValueError("Destination too narrow")
        row=builder.alloc_qubits(length)
        builder.reset(row)
        for i in range(first,n):
            if source[i]==multiplier[j]:
                builder.cx(source[i],row[i-first])
            else:
                builder.init_and(source[i],multiplier[j],row[i-first])
        if subtract != inverse:
            builder.isub(row,target=dest[start:])
        else:
            builder.iadd(row,target=dest[start:])
        for i in reversed(range(first,n)):
            if source[i]==multiplier[j]:
                builder.cx(source[i],row[i-first])
            elif cleanup=="mbuc":
                builder.del_and(source[i],multiplier[j],row[i-first])
            else:
                builder.ccx(source[i],multiplier[j],row[i-first])
        builder.free(row)


@lru_cache(maxsize=8)
def packed_horner(guard):
    mat=coefficient_matrix(COEFF_CONFIG)
    widths=[max(1,int(mat[:,i].max()).bit_length()) for i in range(3)]
    wa,wb,wc=widths
    widths=(wa+guard, wb+guard, wc)
    table=[(int(a)<<guard) |
           ((int(b)<<guard)<<widths[0]) |
           (int(c)<<(widths[0]+widths[1]))
           for a,b,c in mat]
    return table,widths


def build(config: HornerConfig):
    result=screen(config)
    if not result["passes"] or result["worst_error_ratio"] >= .9:
        raise ValueError("Finite-domain Horner error screen failed: "+str(result))
    import kickmix as km
    from .hybrid_kickmix import _normalize, _denormalize
    g=config.guard
    table,widths=packed_horner(g)
    lo=COEFF_CONFIG.low_bits
    b=km.CircuitBuilder()
    x=b.create_quantum_register(NORMAL_BITS,name="x")
    norm,undo=((normalize_barrel,denormalize_barrel)
               if config.normalization=="barrel" else (normalize_compact,denormalize_compact)
               if config.normalization=="compact" else (_normalize,_denormalize))
    with b.mark("normalize"):
        flags,shift=norm(b,x)
    coeff=b.alloc_qubits(sum(widths))
    with b.mark("lookup"):
        b.init_lookup(table=table,address=x[lo:],target=coeff)
    wa,wb,wc=widths
    a=coeff[:wa]
    beta=coeff[wa:wa+wb]
    gamma=coeff[wa+wb:]
    with b.mark("horner_inner_R"):
        _inplace_product(b,gamma,x[:lo],beta,lo-g,
                         cleanup=config.cleanup)
    with b.mark("horner_outer_Y"):
        _inplace_product(b,beta,x[:lo],a,lo,cleanup=config.cleanup)
    with b.mark("phase"):
        _phase(b,a,shift,COEFF_CONFIG.precision+g,+1)
    with b.mark("undo_horner_outer_Y"):
        _inplace_product(b,beta,x[:lo],a,lo,
                         cleanup=config.cleanup,inverse=True)
    with b.mark("undo_horner_inner_R"):
        _inplace_product(b,gamma,x[:lo],beta,lo-g,
                         cleanup=config.cleanup,inverse=True)
    with b.mark("unlookup"):
        b.del_lookup(table=table,address=x[lo:],target=coeff)
        b.free(coeff)
    with b.mark("unnormalize"):
        undo(b,x,flags,shift)
    return b.finish_circuit()
