"""Original co-refit of intercept trailing zero bits with B/C under exact cuts.

The 26-fractional-bit intercept A normally allocates 27 qubits. Its
low 1-2 bits can potentially be made identically zero by jointly
refitting linear B and curvature C against *exact truncated-product*
arithmetic, not by independently rounding a continuous Chebyshev fit.

Each valid normalized interval is fitted on a stratified mesh and
independently screened on all 8192 residuals. Only passing source
coefficients can proceed to a quantum circuit. The ROM stores A>>K,
and the intercept Z_POW uses fractional scale Q-K.
"""
from __future__ import annotations
import argparse,json
from functools import lru_cache
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product

CUTS=(10,13,12)
LO=cc.LO;Q=cc.Q;TOL=cc.TOL;BZ=cc.BZ;CZ=cc.CZ
H=1<<LO
DENSE=np.arange(H,dtype=np.int64)
MESH=np.unique(np.r_[np.linspace(0,H-1,200,dtype=np.int64),
         np.arange(0,256),np.arange(H-256,H),
         np.arange(0,H,107)]).astype(np.int64)

def error(j,A,B,C,r,sq=None):
    if sq is None:sq=truncated_product(r,r,CUTS[1])
    linear=truncated_product(np.int64(B>>BZ),r,CUTS[0]-BZ)
    quad=truncated_product(np.int64(C>>CZ),sq,CUTS[2]-CZ)
    z=(j+r/H)/(1<<cc.P)
    phase=A/(1<<Q)-linear/2**(Q+LO-CUTS[0])+quad/2**(Q+2*LO-CUTS[1]-CUTS[2])
    return np.abs(phase*np.sqrt(z)-1)/TOL

def search_row(j,baseline,K,r,square,span_b,span_c):
    orig_A,orig_B,orig_C=map(int,baseline)
    best=(float("inf"),None)
    z=(j+r/H)/(1<<cc.P)
    weight=np.sqrt(z)/4
    ideal=(1<<Q)/np.sqrt(z)
    for dc in range(-span_c,span_c+1):
        C=orig_C+dc*(1<<CZ)
        if C<0:continue
        quad=truncated_product(np.int64(C>>CZ),square,CUTS[2]-CZ)
        for db in range(-span_b,span_b+1):
            B=orig_B+db*(1<<BZ)
            if B<0:continue
            lin=truncated_product(np.int64(B>>BZ),r,CUTS[0]-BZ)
            wanted=ideal+lin/(1<<(LO-CUTS[0]))-quad/(1<<(2*LO-CUTS[1]-CUTS[2]))
            candidate=int(np.floor(((wanted.max()+wanted.min())*.5)/(1<<K)+.5))*(1<<K)
            for A in (candidate-(1<<K),candidate,candidate+(1<<K)):
                if A<0:continue
                val=float(np.max(np.abs((A-wanted)*weight)))
                if val<best[0]:
                    best=(val,(A,B,C))
    return best

@lru_cache(maxsize=3)
def coefficients(K):
    if K not in (1,2,3):raise ValueError("A trim K must be 1..3")
    original=cc.coefficients(CUTS)
    result=np.array(original,copy=True)
    sq_mesh=truncated_product(MESH,MESH,CUTS[1])
    sq_dense=truncated_product(DENSE,DENSE,CUTS[1])
    failure_log=[]
    for j in range(1<<cc.P,4<<cc.P):
        score,abc=search_row(j,original[j],K,MESH,sq_mesh,span_b=5,span_c=4)
        if abc is None:raise RuntimeError("Refit no feasible candidate")
        actual=float(np.max(error(j,*abc,DENSE,sq_dense)))
        if actual>=1:
            # Coherent slope and curvature shifts can offset a quantized A,
            # but A's value at residual=0 is a strict, independent bound.
            score,abc=search_row(j,original[j],K,DENSE,sq_dense,span_b=14,span_c=8)
            actual=float(np.max(error(j,*abc,DENSE,sq_dense)))
        result[j]=abc
        if actual>=1: failure_log.append((j,actual))
    result.setflags(write=False)
    return result,tuple(failure_log)

def full_screen(K):
    c,failed=coefficients(K)
    square=truncated_product(DENSE,DENSE,CUTS[1])
    worst=0.;idx=0;failures=0
    for j in range(1<<cc.P,4<<cc.P):
        A,B,C=map(int,c[j])
        v=error(j,A,B,C,DENSE,square)
        failures+=int(np.count_nonzero(v>=1))
        k=int(np.argmax(v))
        if float(v[k])>worst:
            worst=float(v[k]);idx=j*H+k
    return dict(trim_low_A_bits=K,passes=failures==0,
        failed_normalized_words=failures,
        rows_with_failures=len(failed),failed_rows=list(failed)[:30],
        worst_error_over_tolerance=worst,
        worst_normalized_word=idx,
        normalized_words=(3<<cc.P)*H,
        A_physical_bits=max(1,int(c[:,0].max()).bit_length())-K,
        B_zeros=BZ,C_zeros=CZ,cuts=CUTS)

def main():
    args=argparse.ArgumentParser()
    args.add_argument('--trims',type=int,nargs='+',default=[1,2])
    a=args.parse_args()
    rows=[]
    for k in a.trims:
        v=full_screen(k);rows.append(v)
        print(json.dumps(v),flush=True)
    path=Path('runs/intercept_trailing_zero_refit.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
