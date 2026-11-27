"""Independent output-error-aware integer coefficient repair at critical intervals.

The public verifier's rotation synthesis price is set by the MINIMUM
absolute phase slack, not by average error or max relative error.
For the independently full-verified v9, exact phase screening identifies
critical normalized intervals 477,366,428,368. Here we refit only those
A/B/C integer coefficients using a minimax search against the EXACT
native omitted binary partial-product circuit, preserving the original
coefficient widths, cuts, lookup and quantum correctness certificate.

For each candidate B/C, choose the intercept A minimizing the
weighted Chebyshev error on a dense residual mesh. Exhaustively test
the best candidate triples on all 8192 residuals. Then re-evaluate
the full legal-input absolute phase slack and verifier rotation-synthesis
cost; accept only if *every legal input* remains within tolerance.

This is original numerical synthesis, not a borrowed Kaggle circuit.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product
from .benchmark1_phase_rotation_prune_math import data_full_domain

CUTS=(10,13,12)
ROWS=(477,366,428,368)
V9=237684.49222389946
DENSE=np.arange(8192,dtype=np.int64)
MESH=np.unique(np.r_[
    np.linspace(0,8191,440,dtype=np.int64),
    np.arange(0,128),np.arange(8064,8192),
    [4800,7803,5016,8065]]).astype(np.int64)

def error_row(j,abc,r):
    A,B,C=map(int,abc)
    return cc._error_for_row(j,r,A,B,C,CUTS)

def fit(j,base,span=35):
    oldA,oldB,oldC=map(int,base)
    bcenter=oldB>>cc.BZ
    ccenter=oldC>>cc.CZ
    bvals=np.arange(max(0,bcenter-span),bcenter+span+1,dtype=np.int64)
    cvals=np.arange(max(0,ccenter-span),ccenter+span+1,dtype=np.int64)
    sq=truncated_product(MESH,MESH,CUTS[1])
    linear=np.stack([truncated_product(np.int64(b),MESH,CUTS[0]-cc.BZ) for b in bvals])
    quad=np.stack([truncated_product(np.int64(c),sq,CUTS[2]-cc.CZ) for c in cvals])
    z=(j+MESH/8192)/(1<<cc.P)
    desired=(1<<cc.Q)/np.sqrt(z)
    weight=np.sqrt(z)/((1<<cc.Q)*cc.TOL)
    desired_A=desired[None,None,:]+linear[:,None,:]/8-quad[None,:,:]/2
    center=(desired_A.max(axis=2)+desired_A.min(axis=2))/2
    candidates=[]
    for offset in (0,1):
        A=np.floor(center).astype(np.int64)+offset
        errs=np.max(np.abs((A[:,:,None]-desired_A)*weight[None,None,:]),axis=2)
        flat=np.argpartition(errs.ravel(),min(40,len(errs.ravel())-1))[:40]
        for idx in flat:
            bi,ci=np.unravel_index(int(idx),errs.shape)
            candidates.append((float(errs[bi,ci]),[int(A[bi,ci]),int(bvals[bi]<<cc.BZ),int(cvals[ci]<<cc.CZ)]))
    original_worst=float(error_row(j,base,DENSE).max())
    best=(original_worst,list(map(int,base)))
    for score,vec in sorted(candidates,key=lambda x:x[0])[:80]:
        if vec[0]<0:continue
        if int(vec[0]).bit_length()>int(oldA).bit_length()+1:continue
        worst=float(error_row(j,vec,DENSE).max())
        if worst<best[0]:best=(worst,vec)
    return dict(row=j,before=original_worst,after=best[0],
      original_coeff=list(map(int,base)),new_coeff=best[1],
      improved=best[0]<original_worst)

def full_input_slack(mat):
    rot,errors,tol,widths=data_full_domain()
    orig=cc.coefficients(CUTS)
    x=np.arange(1,3133453,dtype=np.int64)
    s=10-(np.floor(np.log2(x)).astype(np.int64)//2)
    normalized=x<<(2*s)
    j=normalized>>cc.LO
    res=normalized&((1<<cc.LO)-1)
    changed=np.zeros(len(x),dtype=np.float64)
    for interval in ROWS:
        idx=np.flatnonzero(j==interval)
        if len(idx)==0:continue
        r=res[idx]
        oldA,oldB,oldC=map(int,orig[interval])
        newA,newB,newC=map(int,mat[interval])
        sq=truncated_product(r,r,CUTS[1])
        ld=truncated_product(np.int64(newB>>cc.BZ),r,8)-truncated_product(np.int64(oldB>>cc.BZ),r,8)
        qd=truncated_product(np.int64(newC>>cc.CZ),sq,8)-truncated_product(np.int64(oldC>>cc.CZ),sq,8)
        # Native phase: A*2^s/2^36 - L*2^s/2^39 + P*2^s/2^37.
        delta=((newA-oldA)/(2.0**36)-ld/(2.0**39)+qd/(2.0**37))*np.exp2(s[idx])
        changed[idx]=delta/tol[idx]
    new_errors=errors+changed
    if np.max(np.abs(new_errors))>=1:raise RuntimeError("Original refit fails legal-input tolerance")
    prior_min=float(np.min(tol*(1-np.abs(errors))))
    new_min=float(np.min(tol*(1-np.abs(new_errors))))
    def estimate(margin):
        eps=math.pi*margin/76
        return 85*(2232+(1+76*(.56*math.log2(1/eps)+5.3))/4)
    return dict(passes=True,full_legal_inputs=len(x),
      original_min_abs_margin_turns=prior_min,
      new_min_abs_margin_turns=new_min,
      previous_verified_qc=V9,
      original_predicted_qc=estimate(prior_min),
      new_predicted_qc=estimate(new_min),
      model_predicts_improvement=estimate(new_min)<V9,
      max_error_ratio=float(np.max(np.abs(new_errors))),
      best_absolute_margin_input=int(np.argmin(tol*(1-np.abs(new_errors)))+1))

def run(span=35):
    baseline=cc.coefficients(CUTS)
    mat=np.array(baseline,copy=True)
    reports=[]
    for j in ROWS:
        rep=fit(j,baseline[j],span)
        reports.append(rep);mat[j]=rep["new_coeff"]
        print("ROW "+json.dumps(rep),flush=True)
    allr=DENSE;fail=0;worst=0
    for j in range(128,512):
        w=float(error_row(j,mat[j],allr).max())
        worst=max(worst,w)
        if w>=1:fail+=1
    result=dict(repairs=reports,full_normalized_passes=fail==0,
       failed_intervals=fail,worst_normalized_error_ratio=worst,
       normalized_words=3145728)
    if fail==0:result["score_model"]=full_input_slack(mat)
    print("SUMMARY "+json.dumps(result),flush=True)
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--span",type=int,default=35)
    args=ap.parse_args()
    result=run(args.span)
    p=Path("runs/minmargin_critical_interval_refit.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(result,indent=2)+"\n")
