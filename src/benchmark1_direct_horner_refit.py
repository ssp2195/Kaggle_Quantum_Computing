"""Independently refitted Horner residual phase oracle, no explicit square.

Normalized function:
    A_j/2**26 - ( (B4_j*2**(g+2) - T(C16_j*r,13-g-4))*r )_cut13
                  / 2**(26+g),
where T omits low partial-product columns (not a rounded multiply).
A, B4, C16 are refitted *jointly* against exactly that integer operation.
For g=4 the complete 3,145,728 normalized grid passes error <2**-24.

This is independent of other Kaggle participants' source/circuits.
The native circuit implements the exact operations via a single Horner
intermediate; it must pass the independent verifier before promotion.
"""
from __future__ import annotations
import argparse,hashlib,json
from functools import lru_cache
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .hybrid_model import MAX_X,TOL
from .hybrid_kickmix import resource_counts
from .benchmark1_compact_normalize import normalize_compact,denormalize_compact
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled

P=7;LO=13;Q=26;G=4;BZ=2;CZ=4;H=1<<LO
MESH=np.unique(np.r_[np.arange(0,H,101),
                     np.arange(H-300,H),np.linspace(0,H-1,101,dtype=int),
                     np.arange(0,256)]).astype(np.int64)
TARGET=2076880.0/20
INCUMBENT=272300.78036771005


def _ideal_seed():
    nodes=.5+.5*np.cos(np.pi*(np.arange(3)+.5)/3)
    vand=np.polynomial.polynomial.polyvander(nodes,2)
    j=np.arange(1<<P,4<<P)
    c=np.linalg.solve(vand,1/np.sqrt((j[:,None]+nodes[None,:])/(1<<P)).T).T
    base=np.zeros((4<<P,3),dtype=np.int64)
    base[1<<P:,0]=np.rint(c[:,0]*(1<<Q)).astype(np.int64)
    base[1<<P:,1]=np.rint(-c[:,1]*(1<<Q)).astype(np.int64)
    base[1<<P:,2]=np.rint(c[:,2]*(1<<Q)).astype(np.int64)
    return base


def _row(j,b4,c16,r):
    s=(b4<<(G+BZ))-truncated_product(c16,r,LO-G-CZ)
    if np.any(s<=0): raise ArithmeticError("Horner derivative became negative")
    prod=truncated_product(s,r,LO)
    z=(j+r/H)/(1<<P)
    desired=(1<<Q)/np.sqrt(z)+prod/(1<<G)
    weight=np.sqrt(z)/((1<<Q)*TOL)
    mid=(desired.min()+desired.max())/2
    guess=(int(np.floor(mid)),int(np.ceil(mid)))
    result=min((float(np.max(np.abs((a-desired)*weight))),(a,b4,c16)) for a in guess)
    return result


@lru_cache(maxsize=1)
def coefficients():
    seed=_ideal_seed()
    out=np.zeros_like(seed)
    allr=np.arange(H,dtype=np.int64)
    for j in range(1<<P,4<<P):
        _,b,c=(int(v) for v in seed[j])
        best=(float('inf'),None)
        for dc in range(-6,7):
            C=round(c/(1<<CZ))+dc
            if C<0:continue
            for db in range(-10,11):
                B=round(b/(1<<BZ))+db
                cand=_row(j,B,C,MESH)
                if cand[0]<best[0]:best=cand
        score=_row(j,int(best[1][1]),int(best[1][2]),allr)
        if score[0]>1:
            B0,C0=(int(v) for v in score[1][1:])
            for dc in range(-10,11):
                C=max(0,C0+dc)
                for db in range(-10,11):
                    B=B0+db
                    cand=_row(j,B,C,allr)
                    if cand[0]<score[0]:score=cand
        out[j]=score[1]
    out.setflags(write=False)
    return out


def full_screen():
    coeff=coefficients()
    r=np.arange(H,dtype=np.int64)
    worst=0.;worst_j=0;worst_r=0;failures=0
    for j in range(1<<P,4<<P):
        A,B,C=(int(v) for v in coeff[j])
        derivative=(B<<(G+BZ))-truncated_product(C,r,LO-G-CZ)
        out=truncated_product(derivative,r,LO)
        z=(j+r/H)/(1<<P)
        err=np.abs((A/(1<<Q)-out/2**(Q+G))*np.sqrt(z)-1)/TOL
        failures+=int(np.count_nonzero(err>=1))
        ix=int(np.argmax(err))
        if float(err[ix])>worst:
            worst=float(err[ix]);worst_j=int(j);worst_r=int(r[ix])
    return dict(passes=failures==0,failures=failures,
        normalized_inputs=(3<<P)*H,worst_error_ratio=worst,
        worst_j=worst_j,worst_r=worst_r,guard=G)


def _inplace_derivative(builder,c16,residual,derivative,drop,inverse=False):
    """Exact selected-column controlled adds into a nonzero register.

    Unlike computing a standalone product, carries must propagate through
    all higher derivative bits; keep target=derivative[start:] intact.
    """
    n,m=len(c16),len(residual)
    js=range(m-1,-1,-1) if inverse else range(m)
    for j in js:
        first=max(0,drop-j)
        if first>=n:continue
        start=j+first-drop
        if inverse:
            builder.iadd(c16[first:],target=derivative[start:],control=residual[j])
        else:
            builder.isub(c16[first:],target=derivative[start:],control=residual[j])


def build(padding=5):
    import kickmix as km
    cert=full_screen()
    if not cert["passes"]:raise ValueError("Full normalized grid failed: "+str(cert))
    co=coefficients()
    a,b,c=(co[:,0],co[:,1],co[:,2])
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (a,b,c)]
    ta=[int(v) for v in a]
    tbc=[int(bv)|(int(cv)<<wb) for bv,cv in zip(b,c)]
    packed=[int(av)|(int(bv)<<wa)|(int(cv)<<(wa+wb))
            for av,bv,cv in zip(a,b,c)]
    bld=km.CircuitBuilder()
    x=bld.create_quantum_register(22,name="x")
    flags,shift=normalize_compact(bld,x)
    word=bld.alloc_qubits(wa+wb+wc)
    bld.init_lookup(table=packed,address=x[LO:],target=word)
    a_reg=word[:wa]
    breg=word[wa:wa+wb]
    c_reg=word[wa+wb:]
    with bld.mark("constant_intercept_phase"):
        phase_with_pad(bld,a_reg,shift,Q,+1,10)
        bld.del_lookup(table=ta,address=x[LO:],target=a_reg)
        bld.free(a_reg)
    with bld.mark("coefficient_horner_inner"):
        zeros=bld.alloc_qubits(G+BZ)
        bld.reset(zeros)
        derivative=zeros+breg
        # B4 << (g+2) is exactly the register value since zeros are LSBs.
        _inplace_derivative(bld,c_reg,x[:LO],derivative,LO-G-CZ)
    with bld.mark("direct_horner_residual_phase"):
        # (derivative*r)>>LO, with exact dropped partial products
        prod=bld.alloc_qubits(len(derivative))
        bld.reset(prod)
        mul_controlled(bld,derivative,x[:LO],prod,LO)
        phase_with_pad(bld,prod,shift,Q+G,-1,padding)
        mul_controlled(bld,derivative,x[:LO],prod,LO,inverse=True)
        bld.free(prod)
    with bld.mark("horner_cleanup"):
        _inplace_derivative(bld,c_reg,x[:LO],derivative,LO-G-CZ,inverse=True)
        bld.free(zeros)
        bld.del_lookup(table=tbc,address=x[LO:],target=word[wa:])
        bld.free(word[wa:])
        denormalize_compact(bld,x,flags,shift)
    return bld.finish_circuit()


def full_padding_certificate():
    co=coefficients()
    width=max(1,int(co[:,1].max()).bit_length())+G+BZ
    maxw=[0]*11
    for first in range(1,MAX_X+1,200000):
        x=np.arange(first,min(MAX_X+1,first+200000),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        s=10-k
        norm=x<<(2*s)
        j=norm>>LO;r=norm&(H-1)
        B=co[j,1];C=co[j,2]
        derivative=(B<<(G+BZ))-truncated_product(C,r,LO-G-CZ)
        prod=truncated_product(derivative,r,LO)
        for shift in range(11):
            mask=s==shift
            if mask.any():
                m=int(prod[mask].max())
                if m>maxw[shift]:maxw[shift]=m
    needed=max((s+maxw[s].bit_length()-width for s in range(11)))
    return dict(passes=needed<=5,required_padding=max(0,needed),
                configured=5,values_by_shift=[int(m) for m in maxw],
                checked=MAX_X)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--math-only",action="store_true")
    p.add_argument("--output",type=Path,default=Path("out/horner_refitted.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    args=p.parse_args()
    res=dict(math=full_screen(),padding=full_padding_certificate())
    if not res["padding"]["passes"]:raise RuntimeError("Insufficient phase padding")
    if not args.math_only:
        circ=build()
        res["native"]=resource_counts(circ)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        circ.to_file(args.output,format="kmb")
        res["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
        if args.verify_sample:
            from .hybrid_verify import verify_samples,boundary_inputs
            res["sample"]=verify_samples(circ,boundary_inputs(random_count=args.verify_sample))
            if res["sample"]["failures"]:raise RuntimeError(str(res["sample"]))
    print(json.dumps(res,indent=2))

if __name__=="__main__":main()
