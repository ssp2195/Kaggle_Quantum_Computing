"""Centered signed-digit quadratic Horner phase oracle on 22 input qubits.

Use v = u - 1/2, where u is the local 13/14-bit mantissa fraction.
For input low mantissa m, flipping its most significant bit gives the
two's-complement representation of d=m-2**(lo-1) on the same qubits.

Fit A-B*v+C*v² directly with Chebyshev interpolation on v in [-1/2,1/2].
Horner evaluation with guard g:
 R <- (B<<g) - CUT_signed(C*d, lo-g)
 Y <- (A<<g) - CUT_signed(R*d, lo)
 apply Y/2**(q+g+k) turns, then reverse every update.
The sign-bit row is done by a single quantum-controlled ADD/SUB instead
of a bank of CCX-created partial-product bits.

MBUC is used ONLY for temporary AND results with both live controls.
Arbitrary accumulator bits are NEVER measured away: their nonquadratic
measurement byproducts would generally require non-Clifford corrections.
"""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from .hybrid_model import MAX_X,TOL,NORMAL_BITS
from .hybrid_kickmix import _phase
from .priority1_normalize import normalize_barrel,denormalize_barrel
from .priority1_exact import truncated_product
from .benchmark1_xag_lookup import apply_anf_lookup


@dataclass(frozen=True)
class CenteredConfig:
    segment_bits: int = 6
    guard: int = 4
    lookup: str = "qrom"
    cleanup: str = "mbuc"
    precision: int = 26

    @property
    def low_bits(self):
        return 20-self.segment_bits

    def __post_init__(self):
        if self.segment_bits not in (6,7):
            raise ValueError("Choose p6 or p7")
        if self.guard not in (3,4,5):
            raise ValueError("Guard bits must be 3..5")
        if self.lookup not in ("qrom","xag"):
            raise ValueError("lookup must be qrom/xag")
        if self.cleanup not in ("mbuc","coherent"):
            raise ValueError("cleanup must be mbuc/coherent")
        if self.precision != 26:
            raise ValueError("Only q=26 screened")


@lru_cache(maxsize=8)
def centered_coefficients(segment_bits,precision=26):
    p,q=segment_bits,precision
    nodes=.5*np.cos(np.pi*(np.arange(3)+.5)/3)
    vander=np.polynomial.polynomial.polyvander(nodes,2)
    mat=np.zeros((4<<p,3),dtype=np.int64)
    for idx in range(1<<p,4<<p):
        z=idx/(1<<p)+(nodes+.5)/(1<<p)
        signed=np.linalg.solve(vander,1/np.sqrt(z))
        values=np.rint(signed*2**q).astype(np.int64)
        mat[idx]=values*np.array([1,-1,1],dtype=np.int64)
    if np.any(mat<0):
        raise ArithmeticError("Centered coefficient unexpectedly negative")
    mat.setflags(write=False)
    return mat


def _signed_cut(coefficient,operand,drop,low_bits):
    """Omitted partial-products of signed two's-complement d."""
    low=operand & ((1<<(low_bits-1))-1)
    sign=((operand>>(low_bits-1))&1)*(1<<(low_bits-1))
    return truncated_product(coefficient,low,drop) - truncated_product(coefficient,sign,drop)


def screen(cfg: CenteredConfig, max_x=MAX_X):
    mat=centered_coefficients(cfg.segment_bits,cfg.precision)
    lo,q,g=cfg.low_bits,cfg.precision,cfg.guard
    worst,worst_x,failed=0.,1,0
    for first in range(1,max_x+1,150000):
        x=np.arange(first,min(max_x+1,first+150000),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        zn=x<<(20-2*k)
        idx=zn>>lo
        low=zn&((1<<lo)-1)
        encoded=low ^ (1<<(lo-1))
        A,B,C=(mat[idx,i] for i in range(3))
        R=(B<<g)-_signed_cut(C,encoded,lo-g,lo)
        Y=(A<<g)-_signed_cut(R,encoded,lo,lo)
        approx=np.ldexp(Y.astype(float)/2**(q+g),-k)
        ratio=np.abs(approx*np.sqrt(x)-1)/TOL
        j=int(np.argmax(ratio))
        if ratio[j]>worst:
            worst,worst_x=float(ratio[j]),int(x[j])
        failed+=int(np.count_nonzero(ratio>1))
    return dict(architecture="centered-signed-digit-horner",
                segment_bits=cfg.segment_bits,guard=g,
                lookup=cfg.lookup,cleanup=cfg.cleanup,
                checked_inputs=max_x,failed_inputs=failed,
                worst_error_ratio=worst,worst_x=worst_x,passes=failed==0)


@lru_cache(maxsize=12)
def packed_centered(p,q,g):
    mat=centered_coefficients(p,q)
    widths=[max(1,int(mat[:,i].max()).bit_length()) for i in range(3)]
    wa,wb,wc=widths[0]+g,widths[1]+g,widths[2]
    values=[(int(a)<<g)|((int(b)<<g)<<wa)|(int(c)<<(wa+wb))
            for a,b,c in mat]
    return values,(wa,wb,wc)


def signed_product_update(builder,source,operand,dest,drop,*,cleanup,inverse=False):
    """In-place subtract signed truncated product, undo in reverse order.

    Non-sign rows: compute ANDs with native kickmix, add/sub exact to the
    entire destination suffix, measurement-correct AND cleanup. Sign row:
    quantum-controlled opposite-direction add, with no AND workspace.
    """
    n,m=len(source),len(operand)
    for j in (range(m-1,-1,-1) if inverse else range(m)):
        first=max(0,drop-j)
        if first>=n: continue
        start=j+first-drop
        if start+n-first>len(dest):
            raise ValueError("Destination too narrow for signed product")
        sign_bit=j==m-1
        if sign_bit:
            # Positive coefficient times a negative two's-complement bit.
            if not inverse:
                builder.iadd(source[first:],target=dest[start:],
                             control=operand[j])
            else:
                builder.isub(source[first:],target=dest[start:],
                             control=operand[j])
            continue
        row=builder.alloc_qubits(n-first+1)
        builder.reset(row)
        for i in range(first,n):
            if source[i]==operand[j]:
                builder.cx(source[i],row[i-first])
            else:
                builder.init_and(source[i],operand[j],row[i-first])
        if not inverse:
            builder.isub(row,target=dest[start:])
        else:
            builder.iadd(row,target=dest[start:])
        for i in reversed(range(first,n)):
            if source[i]==operand[j]:
                builder.cx(source[i],row[i-first])
            elif cleanup=="mbuc":
                builder.del_and(source[i],operand[j],row[i-first])
            else:
                builder.ccx(source[i],operand[j],row[i-first])
        builder.free(row)


def build(cfg: CenteredConfig):
    report=screen(cfg)
    if not report["passes"] or report["worst_error_ratio"]>=.9:
        raise ValueError("Full finite-domain phase certificate failed: "+str(report))
    import kickmix as km
    p,q,g,lo=cfg.segment_bits,cfg.precision,cfg.guard,cfg.low_bits
    tab,widths=packed_centered(p,q,g)
    b=km.CircuitBuilder()
    x=b.create_quantum_register(NORMAL_BITS,name="x")
    flags,shift=normalize_barrel(b,x)
    coeff=b.alloc_qubits(sum(widths))
    with b.mark("coefficients_XAG_or_QROM"):
        if cfg.lookup=="xag":
            apply_anf_lookup(b,x[lo:],coeff,tab)
        else:
            b.init_lookup(table=tab,address=x[lo:],target=coeff)
    wa,wb,wc=widths
    A,B,C=coeff[:wa],coeff[wa:wa+wb],coeff[wa+wb:]
    # Convert mantissa offset m to signed v=m/2**lo-1/2 in place.
    b.x(x[lo-1])
    with b.mark("centered_inner_horner"):
        signed_product_update(b,C,x[:lo],B,lo-g,cleanup=cfg.cleanup)
    with b.mark("centered_outer_horner"):
        signed_product_update(b,B,x[:lo],A,lo,cleanup=cfg.cleanup)
    with b.mark("phase"):
        _phase(b,A,shift,q+g,+1)
    with b.mark("reverse_centered_outer"):
        signed_product_update(b,B,x[:lo],A,lo,cleanup=cfg.cleanup,inverse=True)
    with b.mark("reverse_centered_inner"):
        signed_product_update(b,C,x[:lo],B,lo-g,cleanup=cfg.cleanup,inverse=True)
    b.x(x[lo-1])
    with b.mark("erase_coefficient_lookup"):
        if cfg.lookup=="xag":
            apply_anf_lookup(b,x[lo:],coeff,tab)
        else:
            b.del_lookup(table=tab,address=x[lo:],target=coeff)
        b.free(coeff)
    denormalize_barrel(b,x,flags,shift)
    return b.finish_circuit()
