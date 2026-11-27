"""Original four-stage zero-block binary search for reversible base-four normalization.

The binary normalization exponent is obtained from four *greedy* block tests:
   if top16 bits are zero, rotate 16 bits (8 base-four pairs);
   if top8 bits are zero, rotate 8 bits (4 pairs);
   if top4 bits are zero, rotate 4 bits (2 pairs);
   if top2 bits are zero, rotate 2 bits (1 pair).
Each flag becomes one binary digit of the exponent, with NO quantum additions.

Zero-window logic is a Boolean prefix-AND with measurement-based uncomputation
of the finished intermediate prefixes while original controls remain live.
After the phase, reverse rotations and recompute the prefixes to erase flags
causally. Reversible, coherent and exact on every nonzero 22-bit input.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from . import benchmark1_cost_cooptim as cc
from .benchmark1_zero_residual_pad import build as build_phase_v5
from .hybrid_kickmix import resource_counts

STAGES=(8,4,2,1)
BASELINE=272300.78036771005

def _lanes(x):
    import kickmix as km
    return km.array([x[i] for i in range(0,22,2)]),km.array([x[i] for i in range(1,22,2)])

def _compute_window_flag(b,high):
    b.x(high)
    intermediates=[]
    previous=high[0]
    for r in high[1:]:
        target=b.alloc_qubits(1)[0]
        b.init_and(previous,r,target)
        intermediates.append((previous,r,target))
        previous=target
    flag=intermediates[-1][2]
    for left,right,q in reversed(intermediates[:-1]):
        b.del_and(left,right,q)
        b.free(q)
    b.x(high)
    return flag

def _erase_window_flag(b,high,flag):
    b.x(high)
    intermediates=[]
    previous=high[0]
    for r in high[1:-1]:
        target=b.alloc_qubits(1)[0]
        b.init_and(previous,r,target)
        intermediates.append((previous,r,target))
        previous=target
    b.del_and(previous,high[-1],flag)
    b.free(flag)
    for left,right,q in reversed(intermediates):
        b.del_and(left,right,q)
        b.free(q)
    b.x(high)

def normalize_binary_block(b,x):
    if len(x)!=22:raise ValueError("Expected the 22-bit input register")
    even,odd=_lanes(x)
    msb_first=[]
    for pairs in STAGES:
        high=x[22-2*pairs:]
        flag=_compute_window_flag(b,high)
        b.left_rotate(even,pairs,control=flag)
        b.left_rotate(odd,pairs,control=flag)
        msb_first.append(flag)
    import kickmix as km
    exponent=km.array(list(reversed(msb_first)))
    return msb_first,exponent

def denormalize_binary_block(b,x,flags,exponent):
    if len(flags)!=4:raise ValueError("Expected four block flags")
    even,odd=_lanes(x)
    for pairs,flag in reversed(list(zip(STAGES,flags))):
        b.right_rotate(even,pairs,control=flag)
        b.right_rotate(odd,pairs,control=flag)
        _erase_window_flag(b,x[22-2*pairs:],flag)

def build():
    orig_n=cc.normalize_compact
    orig_u=cc.denormalize_compact
    cc.normalize_compact=normalize_binary_block
    cc.denormalize_compact=denormalize_binary_block
    try:
        circuit,calls=build_phase_v5()
    finally:
        cc.normalize_compact=orig_n
        cc.denormalize_compact=orig_u
    return circuit,calls

def numeric_normalization_check():
    # This is the exact classical map implied by the four reversible tests.
    N=3133452
    max_shift=0
    for first in range(1,N+1,200000):
        import numpy as np
        x=np.arange(first,min(N+1,first+200000),dtype=np.int64)
        v=x.copy();shift=np.zeros_like(x)
        for p in STAGES:
            condition=(v>>(22-2*p))==0
            v=np.where(condition,v<<(2*p),v)
            shift+=condition*p
        expected=(21-np.floor(np.log2(x)).astype(np.int64))//2
        if not np.array_equal(shift,expected):
            bad=int(x[np.nonzero(shift!=expected)[0][0]])
            raise ValueError("Binary block normalization mismatch x="+str(bad))
        if not np.array_equal(v,x<<(2*expected)):
            raise ValueError("Normalized mantissa mismatch")
        max_shift=max(max_shift,int(shift.max()))
    return dict(passes=True,checked=N,max_shift=max_shift,
       count_of_quantum_additions_for_exponent=0)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/binary_block_normalized.kmb"))
    p.add_argument("--math-only",action="store_true")
    p.add_argument("--verify-sample",type=int,default=256)
    a=p.parse_args()
    r=dict(classical_normalization=numeric_normalization_check())
    if not a.math_only:
        c,calls=build()
        r["native"]=resource_counts(c)
        r["phase_calls"]=calls
        r["toffoli_floor_below_current"]=r["native"]["lower_bound_qc"]<BASELINE
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        r["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
        if a.verify_sample:
            from .hybrid_verify import verify_samples,boundary_inputs
            sample=verify_samples(c,boundary_inputs(random_count=a.verify_sample))
            r["sample"]=sample
            if sample["failures"]:raise RuntimeError("Coherence failed: "+str(sample))
    print(json.dumps(r,indent=2))
if __name__=="__main__":main()
