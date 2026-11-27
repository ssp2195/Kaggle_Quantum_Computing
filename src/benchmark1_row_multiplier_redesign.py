"""Native cost comparisons of independently scheduled exact partial-product rows.

Uses the same integer arithmetic and independently certified p7/q26
coefficients as v5. Changes only multiplier control/dataflow.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from . import benchmark1_cost_cooptim as cc
from .benchmark1_zero_residual_pad import build as build_pad
from .priority1_truncated import _mul_cut
from .hybrid_kickmix import resource_counts

REFERENCE=2076880.0
INCUMBENT=272300.78036771005

def make(*,linear_row=False,quadratic_row=False):
    original=cc.mul_controlled
    call_counter=[0]
    def wrapped(builder,a,b,dest,drop,inverse=False):
        is_linear=(drop==8 and len(a)>len(b))
        # Linear is B(18) x r(13), quadratic is sq(13) x gamma(10).
        # Both use drop=8 for this precise saved circuit. Determine
        # identity by the operand widths, not by drop.
        is_linear=len(a)>len(b) and len(a)>13
        use=(linear_row and is_linear) or (quadratic_row and not is_linear)
        call_counter[0]+=1
        if use:
            return _mul_cut(builder,a,b,dest,drop,inverse=inverse)
        return original(builder,a,b,dest,drop,inverse=inverse)
    cc.mul_controlled=wrapped
    try:
        circuit,calls=build_pad()
    finally:
        cc.mul_controlled=original
    if call_counter[0]!=4:raise ValueError("Unexpected multiplier call count: "+str(call_counter))
    return circuit

def main():
 rows=[]
 for L,Q in [(False,False),(True,False),(False,True),(True,True)]:
    try:
     c=make(linear_row=L,quadratic_row=Q)
     costs=resource_counts(c)
     record=dict(linear_row=L,quadratic_row=Q,status="NATIVE_BUILT",
       **costs,can_beat_incumbent_by_toffoli_floor=costs["lower_bound_qc"]<INCUMBENT,
       can_reach_20x_by_toffoli_floor=costs["lower_bound_qc"]<REFERENCE/20)
     if record["can_beat_incumbent_by_toffoli_floor"] and (L or Q):
      from .hybrid_verify import verify_samples,boundary_inputs
      rec=verify_samples(c,boundary_inputs(random_count=128))
      record["coherent_sample"]=rec
    except Exception as e:
     record=dict(linear_row=L,quadratic_row=Q,status="FAILED",
        exception=type(e).__name__,message=str(e)[:500])
    print(json.dumps(record),flush=True)
    rows.append(record)
 out=Path('runs/multiplier_stage_comparison.json')
 out.parent.mkdir(parents=True,exist_ok=True)
 out.write_text(json.dumps(rows,indent=2)+'\n')

if __name__=="__main__":main()
