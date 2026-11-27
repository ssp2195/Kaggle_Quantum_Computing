"""Native follow-through for independently co-refitted nonuniform exact truncation.

The certified v9 architecture is held fixed. For every candidate cut
triple (linear, square, quadratic), refit A/B/C against EXACT dropped
binary partial products, choose 128 unreachable table rows to cancel
top ANF interactions, validate all 3,145,728 normalized values, and
prove that the physical left-rotated phase words do not wrap for any
of the 3,133,452 legal x values.

Only after those proofs do we replace the v9 resource-gate circuit
and check the native Q/Toff/rotation cost and sampled coherent behavior.
No borrowed competition code or external circuit is incorporated.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim_exp as exp
from . import benchmark1_pruned_anf_rom as rom
from . import benchmark1_dontcare_boolean_rom as dc
from . import benchmark1_first_row_specialization as first
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts

CONFIGS=((11,13,12),(10,14,12),(10,13,13),
         (11,14,13),(11,14,12),(11,13,13),(10,14,13),(12,13,12))
CURRENT=241424.49222389946
REF=2076880.

def encoded_coefficients(cuts):
    mat=exp.coefficients(cuts)
    A,B,C=mat[:,0],mat[:,1]>>exp.BZ,mat[:,2]>>exp.CZ
    wa,wb,wc=[max(1,int(x.max()).bit_length()) for x in (A,B,C)]
    table=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb)) for a,b,c in zip(A,B,C)]
    encoded=dc.extend_packed(table,"affine")
    m=np.array(mat,copy=True)
    for j in range(128):
        v=int(encoded[j])
        a=v&((1<<wa)-1)
        b=(v>>wa)&((1<<wb)-1)
        c=(v>>(wa+wb))&((1<<wc)-1)
        m[j]=a,b<<exp.BZ,c<<exp.CZ
    assert np.array_equal(m[128:],mat[128:])
    m.setflags(write=False)
    return m,encoded

def pad_certificate(cuts,m,chunk=200000):
    dl,ds,dq=cuts
    beta=m[:,1]>>exp.BZ
    gamma=m[:,2]>>exp.CZ
    wb=max(1,int(beta.max()).bit_length())
    wc=max(1,int(gamma.max()).bit_length())
    wl=wb+exp.LO-(dl-exp.BZ)
    wsq=2*exp.LO-ds
    wq=wc+wsq-(dq-exp.CZ)
    ml=[0]*11;mq=[0]*11
    for start in range(1,3133453,chunk):
        x=np.arange(start,min(3133453,start+chunk),dtype=np.int64)
        shift=10-(np.floor(np.log2(x)).astype(np.int64)//2)
        norm=x<<(2*shift)
        addr=norm>>exp.LO
        rem=norm & ((1<<exp.LO)-1)
        linear=truncated_product(beta[addr],rem,dl-exp.BZ)
        sq=truncated_product(rem,rem,ds)
        quad=truncated_product(gamma[addr],sq,dq-exp.CZ)
        for s in range(11):
            sel=shift==s
            if sel.any():
                ml[s]=max(ml[s],int(linear[sel].max()))
                mq[s]=max(mq[s],int(quad[sel].max()))
    pl=max(0,max(i+v.bit_length()-wl for i,v in enumerate(ml)))
    pq=max(0,max(i+v.bit_length()-wq for i,v in enumerate(mq)))
    return dict(passes=pl<=5 and pq<=5,
          linear_pad=pl,quadratic_pad=pq,
          linear_word_bits=wl,quadratic_word_bits=wq,
          max_linear_by_shift=ml,max_quadratic_by_shift=mq)

def build(cuts,table,encoded):
    oldexp,oldcut=rom.cc,rom.CUTS
    oldopt=dc.optimized_coefficients
    rom.cc=exp
    rom.CUTS=cuts
    dc.optimized_coefficients=lambda policy:(table,encoded)
    try:
        circuit=first.build(True,True)
    finally:
        rom.cc=oldexp
        rom.CUTS=oldcut
        dc.optimized_coefficients=oldopt
    return circuit

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--max',type=int,default=8)
    p.add_argument('--samples',type=int,default=128)
    a=p.parse_args()
    rows=[]
    for cuts in CONFIGS[:a.max]:
        try:
            proof=exp.full_screen(cuts)
            report=dict(cuts=cuts,math=proof)
            if not proof["passes"]:
                report["stage"]="REJECTED_INEXACT_MATH"
            else:
                m,encoded=encoded_coefficients(cuts)
                pads=pad_certificate(cuts,m)
                report["pads"]=pads
                if not pads["passes"]:
                    report["stage"]="REJECTED_PHASE_REGISTER_OVERFLOW"
                else:
                    c=build(cuts,m,encoded)
                    native=resource_counts(c)
                    report["native"]=native
                    report["can_beat_v8_by_toffoli_floor"]=native["lower_bound_qc"]<CURRENT
                    report["can_reach_20x_by_floor"]=native["lower_bound_qc"]<REF/20
                    report["stage"]="NATIVE_BUILT"
                    if report["can_beat_v8_by_toffoli_floor"] and a.samples:
                        from .hybrid_verify import boundary_inputs,verify_samples
                        report["sample"]=verify_samples(c,boundary_inputs(random_count=a.samples))
        except Exception as e:
            report=dict(cuts=cuts,stage="FAILED",error=type(e).__name__+":"+str(e)[:450])
        rows.append(report)
        print(json.dumps(report),flush=True)
    out=Path('runs/alt_cut_native_screen.json')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=="__main__":main()
