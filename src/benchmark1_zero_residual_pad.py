"""Exact input-dependent padding elimination for Benchmark I.

Original contribution: the residual-dependent phase terms vanish whenever
normalization left-shifts the 22-bit input by >=14 positions. Therefore
their phase word is identically 0 for shifts s>=7, and five *upper* zero
padding qubits suffice for rotations at every s where the word is nonzero.

The constant phase keeps all ten padding qubits. The denominator of Z_POW
is NOT changed when trimming the register: it remains 2**(fraction+10).
This ensures exactly the same diagonal oracle on all valid input states.

No Kaggle participant circuit, source or table is used.
"""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import argparse, hashlib, json

from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts
from .hybrid_model import MAX_X

PADS=(10,5,5)
INCUMBENT_QC=283167.97477128694
REFERENCE_QC=2076880.0
CUTS=(10,13,12)


def residual_zero_in_short_pad_region(x: int) -> bool:
    """Mathematical implication checked independently of any coefficient."""
    if not 1<=x<=MAX_X:
        raise ValueError("Input outside Benchmark I range")
    s=(22-x.bit_length())//2
    assert 0<=s<=10
    return s<7 or ((x << (2*s)) & ((1<<13)-1))==0


def full_pad_certificate():
    """Exact full-input word-height proof for the chosen shortened pads.

    Strictly stronger than the zero-residual s>=7 implication: at s=6,
    the residual words have a spare leading zero bit. This is measured
    from the exact truncated integer products at all legal inputs.
    """
    from .benchmark1_min_pad_proof import analyze
    cert=analyze()
    feasible=(cert["checked"]==MAX_X and
        cert["min_linear_pad"]<=PADS[1] and
        cert["min_quadratic_pad"]<=PADS[2])
    return dict(cert,passes=bool(feasible),
                phase_pad_qubits=list(PADS),
                proof="Every shifted exact product fits without cyclic wrap for every x")

def phase_with_pad(builder, word, shift, frac_bits, sign, pad_bits):
    pad=builder.alloc_qubits(pad_bits)
    builder.reset(pad)
    extended=word+pad
    builder.left_rotate(extended, shift=shift)
    # Keep the original 10-bit scale factor even when the physical pad is
    # shorter. This is critical: changing it would change the oracle.
    denominator=1<<(frac_bits+10)
    for j in range(len(extended)):
        exp=Fraction(sign*(1<<(j+1)),denominator)
        if exp%2:
            builder.z_pow(extended[j],exp)
    builder.right_rotate(extended,shift=shift)
    builder.free(pad)


def build():
    original=cc._phase
    calls=[]
    def replacement(builder,word,shift,frac_bits,sign):
        index=len(calls)
        if index>=3:
            raise RuntimeError("Unexpected fourth phase in fixed architecture")
        calls.append((len(word),frac_bits,sign))
        phase_with_pad(builder,word,shift,frac_bits,sign,PADS[index])
    cc._phase=replacement
    try:
        circuit=cc.build(CUTS,swap_quadratic=True)
    finally:
        cc._phase=original
    if len(calls)!=3 or [v[2] for v in calls]!=[1,-1,1]:
        raise RuntimeError("Phase schedule changed: "+str(calls))
    return circuit,calls


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/phase_pad5_exact.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    p.add_argument("--only-math",action="store_true")
    a=p.parse_args()
    certificate=full_pad_certificate()
    if not certificate["passes"]:
        raise RuntimeError(str(certificate))
    result=dict(padding=certificate,normalized_approximation=cc.full_screen(CUTS))
    if not a.only_math:
        circuit,calls=build()
        s=resource_counts(circuit)
        result["native"]=s
        result["phase_calls"]=calls
        result["may_beat_incumbent"]=s["lower_bound_qc"]<INCUMBENT_QC
        a.output.parent.mkdir(parents=True,exist_ok=True)
        circuit.to_file(a.output,format="kmb")
        result["output"]=str(a.output)
        result["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
        if a.verify_sample:
            from .hybrid_verify import verify_samples,boundary_inputs
            result["native_sample"]=verify_samples(
                circuit,boundary_inputs(random_count=a.verify_sample))
            if result["native_sample"]["failures"]:
                print(json.dumps(result,indent=2))
                raise SystemExit("Sampled quantum coherence verification failed")
    print(json.dumps(result,indent=2))
if __name__=="__main__":
    main()
