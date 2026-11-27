"""Native-only qubit high-water profiling with reversible-stage ablation.

Build variants of the original full-verified v9 constructor that omit
exactly one stage. These incomplete circuits are NOT submissions and
MUST NOT be used for any coherent accuracy or Q*C improvement claims.

Only purpose: identify the largest native qubit high-water contributor.
Measured with the precise pinned KickMix circuit.num_qubits and
operation_statistics, not an estimate of Python live registers.

Stage switches include all three phase pad sites, two controlled
multiply stages, symmetric square, and coefficient ROM. Stage skip
compiles a semantically incomplete native circuit intentionally.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from . import benchmark1_pruned_anf_rom as rom
from . import benchmark1_first_row_specialization as row
from .hybrid_kickmix import resource_counts

SHAPES={'linear':(16,13),'quad':(13,7)}
CASES=('baseline','no_intercept_phase','no_linear_phase',
       'no_quadratic_phase','no_all_phases',
       'no_linear_multiply','no_quadratic_multiply',
       'no_square','no_all_arithmetic','no_anf_rom',
       'only_ROM','only_intercept','only_linear','only_quadratic')

def build(case):
    if case not in CASES:raise ValueError(case)
    original_p=rom.phase_with_pad
    original_m=row.mul_first_row
    original_s=row.square_first_row
    original_r=rom.binary_pruned_rom
    phase_calls=[0]
    def phase_fn(b,word,shift,frac,sign,pad):
        n=phase_calls[0];phase_calls[0]+=1
        if case in ('no_all_phases','no_all_arithmetic','only_ROM'):return
        if case in ('no_intercept_phase','only_linear','only_quadratic') and n==0:return
        if case in ('no_linear_phase','only_intercept','only_quadratic') and n==1:return
        if case in ('no_quadratic_phase','only_intercept','only_linear') and n==2:return
        return original_p(b,word,shift,frac,sign,pad)
    def mul_fn(b,a,other,dest,drop,inverse=False):
        if case in ('no_all_arithmetic','only_ROM'):return
        if case in ('no_linear_multiply','only_intercept','only_quadratic') and (len(a),len(other))==SHAPES['linear']:return
        if case in ('no_quadratic_multiply','only_intercept','only_linear') and (len(a),len(other))==SHAPES['quad']:return
        return original_m(b,a,other,dest,drop,inverse)
    def sq_fn(b,r,dest,drop,inverse=False):
        if case in ('no_square','no_all_arithmetic','only_ROM','only_intercept','only_linear'):return
        return original_s(b,r,dest,drop,inverse)
    def decoder(b,addr,word,table):
        if case=='no_anf_rom':return
        return original_r(b,addr,word,table)
    rom.phase_with_pad=phase_fn
    row.mul_first_row=mul_fn
    row.square_first_row=sq_fn
    rom.binary_pruned_rom=decoder
    try:c=row.build(True,True)
    finally:
        rom.phase_with_pad=original_p
        row.mul_first_row=original_m
        row.square_first_row=original_s
        rom.binary_pruned_rom=original_r
    return c

def main():
    records=[]
    for case in CASES:
        try:
            c=build(case)
            result=dict(case=case,profiling_only_invalid_oracle=case!='baseline',
                   native=resource_counts(c))
        except Exception as e:
            result=dict(case=case,error=type(e).__name__+":"+str(e)[:500])
        records.append(result)
        print(json.dumps(result),flush=True)
    p=Path("runs/native-stage-qubit-peak-profiles.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(records,indent=2)+'\n')
if __name__=="__main__":main()
