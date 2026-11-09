"""End-to-end 8-bit XAG / in-place Horner inverse-square-root phase oracle.

This is a full-width 22-bit candidate, not a toy lookup. It uses a
Chebyshev degree-2 approximation over 192 normalized intervals (p=6).
Coefficients A/B are loaded into q+guard fixed-point words. The Horner
updates are exact reversible updates using selected partial-product bits.

The ANF-XAG approach is a valid shared-prefix XOR/AND implementation;
it is not claimed to be an optimized multi-output 172-AND synthesis.
Two coefficient lookup modes are compared under otherwise identical
arithmetic: built-in kickmix QROM versus explicit reversible XAG.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .priority1_model import Approximation, coefficient_matrix
from .priority1_exact import truncated_product
from .benchmark1_horner import _inplace_product
from .benchmark1_xag_lookup import apply_anf_lookup
from .priority1_normalize import normalize_barrel, denormalize_barrel
from .hybrid_kickmix import _phase
from .hybrid_model import MAX_X, TOL, NORMAL_BITS

APPROX = Approximation("chebyshev", 2, 6, 26)


@dataclass(frozen=True)
class XAGHornerConfig:
    guard: int = 5
    lookup: str = "xag"
    cleanup: str = "mbuc"

    def __post_init__(self):
        if self.guard not in (4, 5, 6):
            raise ValueError("Guard width must be 4..6")
        if self.lookup not in ("xag", "qrom"):
            raise ValueError("Lookup must be xag or qrom")
        if self.cleanup not in ("mbuc", "coherent"):
            raise ValueError("Cleanup must be mbuc or coherent")


def screen(config: XAGHornerConfig, max_x=MAX_X):
    """Numerical proof-of-concept on the complete 3,133,452-input domain."""
    m=coefficient_matrix(APPROX)
    lo,q,g=APPROX.low_bits,APPROX.precision,config.guard
    worst,worst_x,failed=0.,1,0
    for start in range(1,max_x+1,150000):
        x=np.arange(start,min(start+150000,max_x+1),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        norm=x<<(20-2*k)
        j=norm>>lo
        u=norm&((1<<lo)-1)
        A,B,C=(m[j,i] for i in range(3))
        R=(B<<g)-truncated_product(C,u,lo-g)
        Y=(A<<g)-truncated_product(R,u,lo)
        approx=np.ldexp(Y.astype(float)/2**(q+g),-k)
        relative=np.abs(approx*np.sqrt(x)-1)/TOL
        index=int(np.argmax(relative))
        if relative[index]>worst:
            worst,worst_x=float(relative[index]),int(x[index])
        failed+=int(np.count_nonzero(relative>1))
    return dict(guard=g,lookup=config.lookup,cleanup=config.cleanup,
                checked_inputs=max_x,failed_inputs=failed,
                worst_error_ratio=worst,worst_x=worst_x,
                passes=failed==0)


@lru_cache(maxsize=8)
def packed_table(guard):
    mat=coefficient_matrix(APPROX)
    w=tuple(max(1,int(mat[:,i].max()).bit_length()) for i in range(3))
    wa,wb,wc=w[0]+guard,w[1]+guard,w[2]
    table=[(int(a)<<guard)|((int(b)<<guard)<<wa)|(int(c)<<(wa+wb))
           for a,b,c in mat]
    if len(table)!=256:
        raise AssertionError("Expected eight address bits")
    return table,(wa,wb,wc)


def build(config: XAGHornerConfig):
    report=screen(config)
    if not report["passes"] or report["worst_error_ratio"]>=.9:
        raise ValueError("XAG-Horner fails phase tolerance: "+str(report))
    import kickmix as km
    lo=APPROX.low_bits
    table,widths=packed_table(config.guard)
    b=km.CircuitBuilder()
    x=b.create_quantum_register(NORMAL_BITS,name="x")
    flags,shift=normalize_barrel(b,x)
    addr=x[lo:]
    storage=b.alloc_qubits(sum(widths))
    with b.mark("eight_bit_coefficient_XAG_or_QROM"):
        if config.lookup=="xag":
            apply_anf_lookup(b,addr,storage,table)
        else:
            b.init_lookup(table=table,address=addr,target=storage)
    wa,wb,wc=widths
    A,B,C=storage[:wa],storage[wa:wa+wb],storage[wa+wb:]
    with b.mark("inplace_horner_R"):
        _inplace_product(b,C,x[:lo],B,lo-config.guard,
                         cleanup=config.cleanup)
    with b.mark("inplace_horner_Y"):
        _inplace_product(b,B,x[:lo],A,lo,cleanup=config.cleanup)
    with b.mark("variable_phase"):
        _phase(b,A,shift,APPROX.precision+config.guard,+1)
    with b.mark("reverse_horner_Y"):
        _inplace_product(b,B,x[:lo],A,lo,cleanup=config.cleanup,
                         inverse=True)
    with b.mark("reverse_horner_R"):
        _inplace_product(b,C,x[:lo],B,lo-config.guard,
                         cleanup=config.cleanup,inverse=True)
    with b.mark("XAG_or_QROM_uncompute"):
        if config.lookup=="xag":
            apply_anf_lookup(b,addr,storage,table)
        else:
            b.del_lookup(table=table,address=addr,target=storage)
        b.free(storage)
    denormalize_barrel(b,x,flags,shift)
    return b.finish_circuit()
