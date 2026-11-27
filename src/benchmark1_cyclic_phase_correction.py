"""Original cyclic phase-word synthesis with exact overflow compensation.

For a w-bit integer A and binary s, 0<=s<=10<w:
    R = rotl_w(A,s) = (A*2**s mod 2**w) + floor(A*2**s/2**w)
    A*2**s = R + (2**w-1)*floor(A*2**s/2**w).

Thus one can apply the phase of A*2**s/2**(Q+10) WITHOUT a
10-qubit high-zero phase-padding bank: cyclically rotate A in place, phase
the rotated bits, and correct the low s wrapped bits by their exact
coefficients (2**w-1)*2**j/2**(Q+10). Each correction is controlled on
(s>j) and A_rot[j], implemented via a temporary AND and measurement-based
Clifford feedback. Undo the cyclic rotation. The operation exactly equals
the previous padded phase for EVERY valid pair (A,s), not just samples.

This is an independent, bit-exact architecture; zero third-party circuit
reuse. The native compiler and complete verifier determine cost and validity.
"""
from __future__ import annotations
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path
from . import benchmark1_cost_cooptim as cc
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .hybrid_kickmix import resource_counts

Q=26;PAD=10;CUTS=(10,13,12)
INCUMBENT=268172.78036771005
REFERENCE=2076880

def integer_identity():
    # Exhaust exact actual intercept coefficients x every s=0,...,10.
    co=cc.coefficients(CUTS)
    wa=max(1,int(co[:,0].max()).bit_length())
    mask=(1<<wa)-1
    worst=0
    for A in map(int,co[:,0]):
        for s in range(11):
            r=((A<<s)&mask)|(A>>(wa-s) if s else 0)
            h=sum(((r>>j)&1)<<j for j in range(s))
            corrected=r+((1<<wa)-1)*h
            if corrected!=(A<<s):
                return dict(passes=False,A=A,s=s,actual=corrected,expected=A<<s)
            worst=max(worst,h)
    return dict(passes=True,coefficient_values_checked=co.shape[0]*11,
        source_word_bits=wa,max_wrapped_word=worst,
        identity="A*2^s=rotl_w(A,s)+(2^w-1)*low_s(rotl_w(A,s))")

def phase_no_padding(builder,word,shift,frac_bits,sign):
    w=len(word)
    if w<11:raise ValueError("Intercept word too narrow for s<=10")
    builder.left_rotate(word,shift=shift)
    denom=1<<(frac_bits+10)
    for i in range(w):
        exp=Fraction(sign*(1<<(i+1)),denom)
        if exp%2:builder.z_pow(word[i],exp)
    for j in range(10):
        flag=builder.alloc_qubits(1)[0]
        builder.reset(flag)
        builder.flip_if_greater_than(shift,j,target=flag)
        control=builder.alloc_qubits(1)[0]
        builder.reset(control)
        builder.init_and(flag,word[j],control)
        correction=Fraction(sign*((1<<w)-1)*(1<<(j+1)),denom)
        if correction%2:builder.z_pow(control,correction)
        builder.del_and(flag,word[j],control)
        builder.free(control)
        builder.flip_if_greater_than(shift,j,target=flag)
        builder.free(flag)
    builder.right_rotate(word,shift=shift)

def build():
    ident=integer_identity()
    if not ident["passes"]:raise ValueError("Exact cyclic identity failure")
    oldphase,oldn,oldu=cc._phase,cc.normalize_compact,cc.denormalize_compact
    calls=[]
    def newphase(b,word,shift,frac,sign):
        idx=len(calls)
        calls.append((len(word),frac,sign))
        if idx==0:
            phase_no_padding(b,word,shift,frac,sign)
        elif idx in (1,2):
            phase_with_pad(b,word,shift,frac,sign,5)
        else:raise ValueError("Unexpected phase")
    cc._phase=newphase
    cc.normalize_compact=normalize_binary_block
    cc.denormalize_compact=denormalize_binary_block
    try:
        c=cc.build(CUTS,swap_quadratic=True)
    finally:
        cc._phase=oldphase
        cc.normalize_compact=oldn
        cc.denormalize_compact=oldu
    return c,calls

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--native',action='store_true')
    p.add_argument('--sample',type=int,default=128)
    p.add_argument('--output',type=Path,default=Path('out/cyclic_compensated_phase.kmb'))
    a=p.parse_args()
    r=dict(identity=integer_identity())
    if a.native:
        c,calls=build()
        stats=resource_counts(c)
        r["resources"]=stats
        r["calls"]=calls
        r["can_beat_v6_floor"]=stats["lower_bound_qc"]<INCUMBENT
        r["can_reach_20x_floor"]=stats["lower_bound_qc"]<REFERENCE/20
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        r["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
        if a.sample and stats["lower_bound_qc"]<INCUMBENT:
            from .hybrid_verify import verify_samples,boundary_inputs
            r["coherent"]=verify_samples(c,boundary_inputs(random_count=a.sample))
            if r["coherent"]["failures"]:raise RuntimeError("Cyclic phase coherence failed")
    print(json.dumps(r,indent=2))

if __name__=="__main__":main()
