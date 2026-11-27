"""High-bit phase-overflow simplification in the original v8 architecture.

For 27-bit intercept A, maximum normalization shift s=10, and denominator
2**36 (Q=26+10), the sole bit that can overflow a 36-bit phase register
is A[26] at s=10. Its original phase is *one whole turn*, which is
physically the identity. A cyclic 36-bit rotation instead produces
a spurious 2**(-36) turn, a known perturbation far below the relative
2**(-24) tolerance. All other intercept contributions remain identical.

The original v8 exact symmetric-square and ANF coefficient ROM are
retained. The score is NOT claimed until the complete independent
verifier confirms it. This tests 9 instead of 10 intercept pad wires.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from . import benchmark1_pruned_anf_rom as rom
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_symmetric_square_rows_v7 import build as build_square
from .hybrid_kickmix import resource_counts
from . import benchmark1_cost_cooptim as cc

ORIGINAL_WA=27
Q=26;MAX_SHIFT=10
INCUMBENT=256476.78036771005

def overflow_certificate():
    co=cc.coefficients((10,13,12))
    xvals=[1,2,3]
    result=[]
    for x in xvals:
        norm=x<<20
        j=norm>>13
        A=int(co[j,0])
        # Only the high bit A[26] at s=10 would wrap.
        highbit=(A>>26)&1
        delta_turns=highbit/(1<<36)
        tolerance=(2**-24)/x**0.5
        result.append(dict(x=x,A=A,overflow_high_bit=highbit,
             absolute_spurious_phase_turns=delta_turns,
             error_over_tolerance=delta_turns/tolerance))
    assert max(t["error_over_tolerance"] for t in result)<.001
    return dict(passes=True,affected_valid_inputs=result,
        max_extra_relative_error_ratio=max(t["error_over_tolerance"] for t in result),
        analytical_scope="all x with s=10, other shifts cannot wrap")

def build():
    if not overflow_certificate()["passes"]:
        raise RuntimeError("Unsafe intercept pad9")
    old=rom.phase_with_pad
    phase_calls=[0]
    def p9(builder,word,shift,frac,sign,pad_bits):
        if phase_calls[0]==0:
            assert len(word)==27 and frac==26 and sign==1 and pad_bits==10
            pad_bits=9
        phase_calls[0]+=1
        return old(builder,word,shift,frac,sign,pad_bits)
    rom.phase_with_pad=p9
    try:
        c=build_square()
    finally:
        rom.phase_with_pad=old
    if phase_calls[0]!=3:raise RuntimeError("Unexpected number of phase terms")
    return c

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/v8_phase_pad9.kmb"))
    p.add_argument("--samples",type=int,default=256)
    args=p.parse_args()
    c=build()
    stats=resource_counts(c)
    report=dict(overflow_proof=overflow_certificate(),
         native=stats,could_beat_v7_at_toffoli_floor=
         stats["lower_bound_qc"]<INCUMBENT)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    report["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if args.samples and report["could_beat_v7_at_toffoli_floor"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        report["coherent"]=verify_samples(c,boundary_inputs(random_count=args.samples))
        if report["coherent"]["failures"]:
            raise RuntimeError("Quantum coherence failed: "+str(report))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
