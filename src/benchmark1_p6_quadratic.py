"""Independent p6 low-table quadratic direct phase with original binary normalization.

Uses 192 real mantissa intervals (256-entry padded QROM), instead of the
verified p7 architecture's 384 (512-entry) intervals. Coefficients are
generated afresh from a Chebyshev minimax approximation of the same
function, and exact partial-product cuts are verified on all possible
normalized words. The quadratic residual is computed directly, phased,
and measurement-corrected uncomputed. No competitor source is imported.
"""
from __future__ import annotations
from functools import lru_cache
from fractions import Fraction
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from .priority1_model import Approximation,coefficient_matrix
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled
from .priority1_truncated import _mul_cut

P=6;LO=20-P;Q=26;H=1<<LO
CUTS=(11,12,10)
TOL=2**-24
MAX_X=3133452
INCUMBENT=268172.78036771005
REFERENCE=2076880.0

@lru_cache(maxsize=1)
def coefficients():
    c=coefficient_matrix(Approximation("chebyshev",2,P,Q))
    c.setflags(write=False)
    return c

def full_normalized_screen():
    c=coefficients();r=np.arange(H,dtype=np.int64)
    square=truncated_product(r,r,CUTS[1])
    worst=0;bad=0;at_j=0;at_r=0
    for j in range(1<<P,4<<P):
        A,B,C=map(int,c[j])
        linear=truncated_product(B,r,CUTS[0])
        quadratic=truncated_product(C,square,CUTS[2])
        z=(j+r/H)/(1<<P)
        phase=A/(1<<Q)-linear/2**(Q+LO-CUTS[0])+quadratic/2**(Q+2*LO-CUTS[1]-CUTS[2])
        error=np.abs(phase*np.sqrt(z)-1)/TOL
        n=int(np.count_nonzero(error>=1));bad+=n
        ix=int(np.argmax(error))
        if float(error[ix])>worst:
            worst=float(error[ix]);at_j=j;at_r=int(r[ix])
    return dict(passes=bad==0,failures=bad,
        normalized_words=3*(1<<P)*H,
        worst_error_ratio=worst,worst_j=at_j,worst_r=at_r,
        cuts=CUTS)

def pad_certificate():
    c=coefficients()
    B,C=c[:,1],c[:,2]
    wb,wc=max(1,int(B.max()).bit_length()),max(1,int(C.max()).bit_length())
    wl=wb+LO-CUTS[0]
    wq=wc+(2*LO-CUTS[1])-CUTS[2]
    max_l=[0]*11;max_q=[0]*11
    for first in range(1,MAX_X+1,200000):
        x=np.arange(first,min(MAX_X+1,first+200000),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        s=10-k
        z=x<<(2*s)
        idx=z>>LO
        r=z&(H-1)
        linear=truncated_product(B[idx],r,CUTS[0])
        square=truncated_product(r,r,CUTS[1])
        quad=truncated_product(C[idx],square,CUTS[2])
        for sh in range(11):
            mask=s==sh
            if mask.any():
                max_l[sh]=max(max_l[sh],int(linear[mask].max()))
                max_q[sh]=max(max_q[sh],int(quad[mask].max()))
    padl=max(0,max(i+v.bit_length()-wl for i,v in enumerate(max_l)))
    padq=max(0,max(i+v.bit_length()-wq for i,v in enumerate(max_q)))
    return dict(legal_inputs=MAX_X,required_linear_pad=padl,
        required_quadratic_pad=padq,width_linear=wl,width_quadratic=wq,
        max_linear_by_shift=max_l,max_quadratic_by_shift=max_q)

def build():
    import kickmix as km
    math=full_normalized_screen()
    if not math["passes"]:raise ValueError("p6 failed stronger normalized math: "+str(math))
    pad=pad_certificate()
    if pad["required_linear_pad"]>10 or pad["required_quadratic_pad"]>10:
        raise ValueError("Pads beyond 10: "+str(pad))
    co=coefficients();A,B,C=co[:,0],co[:,1],co[:,2]
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    ta=[int(v) for v in A]
    tbc=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    table=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
            for a,b,c in zip(A,B,C)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(b,x)
    reg=b.alloc_qubits(wa+wb+wc)
    b.init_lookup(table=table,address=x[LO:],target=reg)
    a=reg[:wa];beta=reg[wa:wa+wb];gamma=reg[wa+wb:]
    with b.mark("coarse_intercept_phase"):
        phase_with_pad(b,a,shift,Q,+1,10)
        b.del_lookup(table=ta,address=x[LO:],target=a)
        b.free(a)
    with b.mark("coarse_linear_phase"):
        dl,ds,dq=CUTS
        prod=b.alloc_qubits(wb+LO-dl);b.reset(prod)
        mul_controlled(b,beta,x[:LO],prod,dl)
        phase_with_pad(b,prod,shift,Q+LO-dl,-1,pad["required_linear_pad"])
        mul_controlled(b,beta,x[:LO],prod,dl,inverse=True)
        b.free(prod)
    with b.mark("coarse_quadratic_phase"):
        sq=b.alloc_qubits(2*LO-ds);b.reset(sq)
        _mul_cut(b,x[:LO],x[:LO],sq,ds)
        prod=b.alloc_qubits(wc+len(sq)-dq);b.reset(prod)
        mul_controlled(b,sq,gamma,prod,dq)
        phase_with_pad(b,prod,shift,Q+2*LO-ds-dq,+1,pad["required_quadratic_pad"])
        mul_controlled(b,sq,gamma,prod,dq,inverse=True)
        b.free(prod)
        _mul_cut(b,x[:LO],x[:LO],sq,ds,inverse=True)
        b.free(sq)
    b.del_lookup(table=tbc,address=x[LO:],target=reg[wa:])
    b.free(reg[wa:])
    denormalize_binary_block(b,x,flags,shift)
    return b.finish_circuit()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native",action="store_true")
    ap.add_argument("--sample",type=int,default=128)
    ap.add_argument("--output",type=Path,default=Path("out/original_p6_quadratic.kmb"))
    a=ap.parse_args()
    report=dict(math=full_normalized_screen(),padding=pad_certificate())
    if a.native and report["math"]["passes"]:
        c=build()
        s=resource_counts(c)
        report["native"]=s
        report["can_beat_incumbent_at_toffoli_floor"]=s["lower_bound_qc"]<INCUMBENT
        report["can_reach_20x_at_toffoli_floor"]=s["lower_bound_qc"]<REFERENCE/20
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        report["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
        if a.sample and s["lower_bound_qc"]<INCUMBENT:
            from .hybrid_verify import verify_samples,boundary_inputs
            report["coherent"]=verify_samples(c,boundary_inputs(random_count=a.sample))
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
