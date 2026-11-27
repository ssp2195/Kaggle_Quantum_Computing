"""Exhaustive greedy deletion of negligible dynamic phase rotations.

Exactly evaluate all 3,133,452 legal-input phases of the currently
full-verifier-pinned v9 circuit, including cyclic phase word wrapping.
For each of its original 3 Z_POW phase stages, try deleting a rotation
on one bit at a time; accept only if EVERY input still stays below
the competition's relative 2^-24 / sqrt(x) tolerance.

At each step pick the deletion that minimizes the new maximum
relative error, repeat until no further single deletion passes.
The constant/slope/curvature bit values are evaluated from precisely
the same integer omitted-column products and current 384-entry
coefficient table used by the native circuit.

This is an independently developed *global phase output pruning*
method. A mathematical screen is a candidate filter only: any edited
native circuit must pass the full pinned 200,540,928-trajectory verifier
and improve actual QxC before production promotion.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
from .benchmark1_phase_support_prune import support_cert
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product

CUTS=(10,13,12)
MAX=3133452
DENS=(36,39,37)
SIGNS=(+1,-1,+1)
PADS=(9,5,5)

def data_full_domain(chunk=200000):
    mat=cc.coefficients(CUTS)
    A=mat[:,0];B=mat[:,1]>>cc.BZ;C=mat[:,2]>>cc.CZ
    wb=max(1,int(B.max()).bit_length())
    wc=max(1,int(C.max()).bit_length())
    widths=(max(1,int(A.max()).bit_length())+9,
            wb+cc.LO-(CUTS[0]-cc.BZ)+5,
            wc+(2*cc.LO-CUTS[1])-(CUTS[2]-cc.CZ)+5)
    out=[np.empty(MAX,dtype=np.int64) for _ in range(3)]
    errors=np.empty(MAX,dtype=np.float64)
    tolerance=np.empty(MAX,dtype=np.float64)
    for start in range(1,MAX+1,chunk):
        z=np.arange(start,min(start+chunk,MAX+1),dtype=np.int64)
        s=10-(np.floor(np.log2(z)).astype(np.int64)//2)
        norm=z<<(2*s)
        j=norm>>cc.LO
        res=norm&((1<<cc.LO)-1)
        square=truncated_product(res,res,CUTS[1])
        vals=(A[j],truncated_product(B[j],res,CUTS[0]-cc.BZ),
                     truncated_product(C[j],square,CUTS[2]-cc.CZ))
        sl=slice(start-1,start-1+len(z))
        for t,(word,width) in enumerate(zip(vals,widths)):
            mask=(1<<width)-1
            out[t][sl]=((word<<s)&mask)|(word>>(width-s))
        phase=out[0][sl]/(1<<DENS[0])-out[1][sl]/(1<<DENS[1])+out[2][sl]/(1<<DENS[2])
        target=1/np.sqrt(z)
        tolerance[sl]=2**-24*target
        errors[sl]=((phase-target+0.5)%1-0.5)/tolerance[sl]
    return out,errors,tolerance,widths

def optimize():
    out,err,tol,widths=data_full_domain()
    base_worst=float(np.max(np.abs(err)))
    original_margin=float(np.min((1-np.abs(err))*tol))
    def estimate_qc(r,margin):
        eps=math.pi*margin/r
        return 85*(2232+(1+r*(.56*math.log2(1/eps)+5.3))/4)
    original_predicted_qc=estimate_qc(76,original_margin)
    print(json.dumps(dict(phase_model="original cyclic-fixed-point native v9",
       checked=MAX,baseline_worst_error_ratio=base_worst,
       original_absolute_min_margin_turns=original_margin,
       original_predicted_qc=original_predicted_qc,
       widths=widths)),flush=True)
    deleted=[]
    positions=[(t,j) for t,w in enumerate(widths) for j in range(w)]
    while True:
        best=None
        for t,j in positions:
            if (t,j) in deleted:continue
            angle=SIGNS[t]*(2.0**(j-DENS[t]))
            if angle==0:continue
            bit=(out[t]>>j)&1
            if not np.any(bit):candidate=(base_worst,(t,j),err.copy())
            else:
                delta=bit*(angle/tol)
                # No periodic branch can change for candidates with
                # max error <1: the correction is << 1/2 turn.
                adj=err-delta
                worst=float(np.max(np.abs(adj)))
                candidate=(worst,(t,j),adj)
            if candidate[0]<1 and (best is None or candidate[0]<best[0]):
                best=candidate
        if best is None:break
        newworst,(t,j),err=best
        deleted.append((t,j))
        base_worst=newworst
        print(json.dumps(dict(deleted_stage=t,deleted_bit=j,
            cumulative_deleted=len(deleted),
            worst_error_ratio=newworst)),flush=True)
        if len(deleted)>=20:break
    new_margin=float(np.min((1-np.abs(err))*tol))
    # A skipped gate is a fine Z rotation only if its exponent is not
    # Clifford (integer quarter turns are modeled separately).
    fine_removed=sum(1 for stage,j in deleted
                     if 0<2.0**(j-DENS[stage])%0.25)
    new_rotations=76-fine_removed
    new_predicted_qc=estimate_qc(new_rotations,new_margin)
    return dict(legal_inputs=MAX,baseline_model=base_worst if not deleted else None,
        original_min_margin_turns=original_margin,
        final_min_margin_turns=new_margin,
        previous_verified_qc=237684.49222389946,
        previous_predicted_qc=original_predicted_qc,
        new_predicted_qc=new_predicted_qc,
        predicted_improvement=new_predicted_qc<237684.49222389946,
        new_rotations=new_rotations,
        deleted_phase_bits=[list(t) for t in deleted],
        deleted_count=len(deleted),
        final_worst_error_ratio=base_worst,
        passes=base_worst<1,original_phase_register_widths=list(widths))

def main():
    r=optimize()
    p=Path("runs/exhaustive_diagonal_phase_pruning.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(r,indent=2)+"\n")
    print("SUMMARY "+json.dumps(r),flush=True)

if __name__=="__main__":main()
