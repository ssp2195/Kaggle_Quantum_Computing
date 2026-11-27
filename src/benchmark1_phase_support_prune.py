"""Full-domain quantum phase word bit-support pruning on verified v9.

The previous phase_with_pad emits Z_POW on every fixed pad/register bit,
even if that bit can never be 1 at its execution point for any legal x.
This fully exhaustive analysis computes OR-support bit masks on all
3,133,452 legal inputs for the actual three exact integer phase terms
A(j), T8(B(j),r), T8(C(j),T13(r,r)), after applying the four-stage
normalization exponent s and cyclic phase register rotation.

If bit j is zero for EVERY valid x, applying its Z_POW is physically
identity on the entire allowed coherent input subspace. We may omit it
without changing any phase, measurement outcome, or workspace cleanup.
All 3,133,452 input values are checked; no probabilistic pruning.

Preserve original 85q v9 carry-safe square, first-row adders, ANF ROM,
and pinned phase error certificate. No competitor code or circuit used.
"""
from __future__ import annotations
import argparse,json,hashlib
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product
from . import benchmark1_pruned_anf_rom as rom
from . import benchmark1_first_row_specialization as first
from .hybrid_kickmix import resource_counts

CUTS=(10,13,12)
MAX=3133452
PADS=(9,5,5)
PRIOR_QC=237684.49222389946
REF=2076880.0

@lru_cache(maxsize=1)
def support_cert(chunk=100000):
    mat=cc.coefficients(CUTS)
    A=mat[:,0]
    B=mat[:,1]>>cc.BZ
    C=mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(t.max()).bit_length()) for t in (A,B,C)]
    dl,ds,dq=CUTS
    widths=[wa,wb+cc.LO-(dl-cc.BZ),wc+(2*cc.LO-ds)-(dq-cc.CZ)]
    W=[w+p for w,p in zip(widths,PADS)]
    masks=[0,0,0]
    counts=[0]*3
    for start in range(1,MAX+1,chunk):
        x=np.arange(start,min(MAX+1,start+chunk),dtype=np.int64)
        s=10-(np.floor(np.log2(x)).astype(np.int64)//2)
        norm=x<<(2*s)
        j=norm>>cc.LO
        r=norm&((1<<cc.LO)-1)
        linear=truncated_product(B[j],r,dl-cc.BZ)
        square=truncated_product(r,r,ds)
        quadratic=truncated_product(C[j],square,dq-cc.CZ)
        for idx,word in enumerate((A[j],linear,quadratic)):
            width=W[idx]
            mod=(1<<width)-1
            shifted=((word<<s)&mod) | (word>>(width-s))
            masks[idx]|=int(np.bitwise_or.reduce(shifted))
            counts[idx]+=len(x)
    for i in range(3):
        if counts[i]!=MAX:raise RuntimeError("Incomplete domain coverage")
        if masks[i]>>W[i]:raise RuntimeError("Phase mask exceeds word width")
    keep=[[j for j in range(w) if (mask>>j)&1]
          for w,mask in zip(W,masks)]
    absent=[[j for j in range(w) if j not in keep[i]]
            for i,w in enumerate(W)]
    return dict(passes=True,legal_inputs=MAX,
       phase_words=widths,pad_words=PADS,physical_widths=W,
       support_masks=[hex(i) for i in masks],
       used_positions=keep,unused_positions=absent,
       all_positions=sum(W),
       suppressed_rotations=sum(map(len,absent)))

def build():
    certificate=support_cert()
    if not certificate["passes"]:raise RuntimeError("Support proof failed")
    previous=rom.phase_with_pad
    calls=[]
    def optimized(b,word,shift,frac,sign,pad):
        idx=len(calls)
        if idx>=3:raise RuntimeError("Unexpected phase call")
        assert len(word)==certificate["phase_words"][idx]
        assert pad==PADS[idx]
        flags=set(certificate["used_positions"][idx])
        filler=b.alloc_qubits(pad);b.reset(filler)
        full=word+filler
        b.left_rotate(full,shift=shift)
        den=1<<(frac+10)
        for j in range(len(full)):
            exponent=Fraction(sign*(1<<(j+1)),den)
            if (exponent%2) and j in flags:
                b.z_pow(full[j],exponent)
        b.right_rotate(full,shift=shift)
        b.free(filler)
        calls.append(idx)
    rom.phase_with_pad=optimized
    try:c=first.build(True,True)
    finally:rom.phase_with_pad=previous
    if len(calls)!=3:raise RuntimeError("Missing phase stage")
    return c

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--native",action="store_true")
    p.add_argument("--samples",type=int,default=256)
    p.add_argument("--output",type=Path,default=Path("out/support_pruned_phase.kmb"))
    args=p.parse_args()
    proof=support_cert()
    report=dict(support_certificate=proof)
    if args.native:
        circuit=build()
        s=resource_counts(circuit)
        report["native"]=s
        report["could_beat_v9_by_toffoli_floor"]=s["lower_bound_qc"]<PRIOR_QC
        report["could_beat_20x_by_toffoli_floor"]=s["lower_bound_qc"]<REF/20
        args.output.parent.mkdir(parents=True,exist_ok=True)
        circuit.to_file(args.output,format="kmb")
        report["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
        if args.samples and report["could_beat_v9_by_toffoli_floor"]:
            from .hybrid_verify import boundary_inputs,verify_samples
            report["coherent"]=verify_samples(circuit,boundary_inputs(random_count=args.samples))
            if report["coherent"]["failures"]:raise RuntimeError("Suppressed necessary phase bit")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
