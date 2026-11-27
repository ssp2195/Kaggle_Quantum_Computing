"""Independent B/C coefficient compression with joint exact-arithmetic refit.

Increase the number of guaranteed zero low bits in the signed-magnitude
slope B and curvature C. New integer candidates are jointly optimized
with intercept A to preserve the complete normalized-domain diagonal
phase under the *actual* omitted-column arithmetic (not rounded
multiplication).

The native B/C coefficient registers may then shrink; whether reduced
width offsets altered multiplier cuts is measured only after a full
3,145,728-word accuracy certificate succeeds.

No Kaggle participant code, table or circuit is used.
"""
from __future__ import annotations
import json,argparse
from functools import lru_cache
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product

CUTS=(10,13,12)
H=8192
DENSE=np.arange(H,dtype=np.int64)
MESH=np.unique(np.r_[np.arange(0,H,47),np.linspace(0,H-1,129,dtype=int),
    np.arange(0,256),np.arange(H-256,H)]).astype(np.int64)
PARAMS=((3,4),(2,5),(3,5),(4,4),(2,6),(4,5))
TOL=cc.TOL
Q=cc.Q;P=cc.P;LO=cc.LO

def row_error(j,bc,cz,A,B,C,r,sq):
    dl,ds,dq=CUTS
    y=(j+r/H)/(1<<P)
    lin=truncated_product(np.int64(B>>bc),r,dl-bc)
    qua=truncated_product(np.int64(C>>cz),sq,dq-cz)
    phase=A/(1<<Q)-lin/(1<<(Q+LO-dl))+qua/(1<<(Q+2*LO-ds-dq))
    return np.abs(phase*np.sqrt(y)-1)/TOL

def fit(j,prior,bc,cz,r,sq,span=7):
    _,baseB,baseC=map(int,prior)
    bcenter=round(baseB/(1<<bc))
    ccenter=round(baseC/(1<<cz))
    dl,ds,dq=CUTS
    y=(j+r/H)/(1<<P)
    wanted=(1<<Q)/np.sqrt(y)
    weight=np.sqrt(y)/((1<<Q)*TOL)
    best=(float('inf'),None)
    for dc in range(-span,span+1):
        C=(ccenter+dc)*(1<<cz)
        if C<0:continue
        qterm=truncated_product(np.int64(C>>cz),sq,dq-cz)/2**(2*LO-ds-dq)
        for db in range(-span,span+1):
            B=(bcenter+db)*(1<<bc)
            if B<0:continue
            lterm=truncated_product(np.int64(B>>bc),r,dl-bc)/2**(LO-dl)
            target=wanted+lterm-qterm
            a=int(np.floor((target.min()+target.max())/2))
            for A in (a-1,a,a+1,a+2):
                if A<0:continue
                loss=float(np.max(np.abs((A-target)*weight)))
                if loss<best[0]:best=(loss,(A,B,C))
    return best

@lru_cache(maxsize=7)
def coefficients(bc,cz):
    source=cc.coefficients(CUTS)
    dest=np.array(source,copy=True)
    mesh_sq=truncated_product(MESH,MESH,CUTS[1])
    dense_sq=truncated_product(DENSE,DENSE,CUTS[1])
    firstfails=[]
    for j in range(128,512):
        _,cand=fit(j,source[j],bc,cz,MESH,mesh_sq)
        if cand is None:raise ValueError("No fitted coeff j"+str(j))
        loss=float(np.max(row_error(j,bc,cz,*cand,DENSE,dense_sq)))
        if loss>=1:
            _,cand2=fit(j,source[j],bc,cz,DENSE,dense_sq,span=12)
            if cand2 is not None:
                loss2=float(np.max(row_error(j,bc,cz,*cand2,DENSE,dense_sq)))
                if loss2<loss:loss,cand=loss2,cand2
        dest[j]=cand
        if loss>=1:firstfails.append((j,loss))
    dest.setflags(write=False)
    return dest,firstfails

def full_screen(bc,cz):
    if not (0<=bc<CUTS[0] and 0<=cz<CUTS[2]):
        raise ValueError("Unrealizable cut")
    mat,fail=coefficients(bc,cz)
    sq=truncated_product(DENSE,DENSE,CUTS[1])
    nfail=0;worst=0.;where=None
    for j in range(128,512):
        A,B,C=(int(x) for x in mat[j])
        error=row_error(j,bc,cz,A,B,C,DENSE,sq)
        nfail+=int(np.count_nonzero(error>=1))
        k=int(np.argmax(error))
        if float(error[k])>worst:
            worst=float(error[k]);where=(j,k)
    A,B,C=mat[:,0],mat[:,1]>>bc,mat[:,2]>>cz
    widths=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    return dict(passes=nfail==0,failures=nfail,
      worst_error_over_tolerance=worst,worst_j_r=where,
      bc_zero_bits=bc,cz_zero_bits=cz,
      coefficient_bits=widths,
      total_coefficient_qubits=sum(widths),
      original_coefficient_bits=[27,11,12],
      rows_with_initial_fail=len(fail),
      failed_rows=fail[:20],normalized_inputs=384*H)

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--max",type=int,default=len(PARAMS))
    args=a.parse_args()
    result=[]
    for bc,cz in PARAMS[:args.max]:
        try:r=full_screen(bc,cz)
        except Exception as e:r=dict(bc_zero_bits=bc,cz_zero_bits=cz,error=type(e).__name__+":"+str(e)[:400])
        result.append(r)
        print(json.dumps(r),flush=True)
    p=Path("runs/quantized_slope_curvature_math.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2)+"\n")
if __name__=="__main__":main()
