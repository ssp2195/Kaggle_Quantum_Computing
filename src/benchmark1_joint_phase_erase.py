"""Original jointly optimized coefficient lifetime and conditional intercept phase.

Previous isolated interventions failed to lower peak width:
  early B erasure alone: 85 qubits, 2317 Toff (old arithmetic);
  phase shortening alone: 85 qubits, 2290+ Toff (old arithmetic).
This combines *both* with full-verified v9 carry-free first-row products
and symmetric squaring. The coefficient B bank is measured and erased
immediately after its final use, while the intercept borrows a shortened
physical high-zero bank and corrects the omitted high-bit phase whenever
its intended contribution is not a whole turn.

No change to mathematical A/B/C table, exact integer polynomial or
all-domain normalization. All candidates undergo native resource gates,
sampled quantum coherent cleanup, and would require full 200.5M shot
independent verification before promotion.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from . import benchmark1_early_coefficient_erasure as early
from . import benchmark1_first_row_specialization as row
from . import benchmark1_phase_conditioned_wrap as wrapped
from .benchmark1_zero_residual_pad import phase_with_pad as original_phase
from .hybrid_kickmix import resource_counts

PRIOR_FULL=237684.49222389946
REF=2076880.

def build(p=6,erase_b="qrom",erase_c="qrom"):
    cert=wrapped.certified_wrap_bound(p)
    if not cert["passes"]:raise RuntimeError("Wrap bound fails")
    old_phase=early.phase_with_pad
    old_mul=early.mul_controlled
    old_sq=early.square_symmetric_rows
    phases=[]
    def phase_fn(b,word,shift,frac,sign,pad):
        if not phases:
            if len(word)!=27 or frac!=26 or sign!=1:
                raise RuntimeError("Unexpected intercept phase")
            wrapped.correction(b,word,shift,frac,p)
        else:
            old_phase(b,word,shift,frac,sign,pad)
        phases.append((len(word),frac,sign))
    early.phase_with_pad=phase_fn
    early.mul_controlled=row.mul_first_row
    early.square_symmetric_rows=row.square_first_row
    try:c=early.build(erase_b=erase_b,erase_c=erase_c)
    finally:
        early.phase_with_pad=old_phase
        early.mul_controlled=old_mul
        early.square_symmetric_rows=old_sq
    if len(phases)!=3:raise RuntimeError("Unexpected phase stage count")
    return c,cert

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--pads",type=int,nargs="+",default=[8,7,6,5,4])
    parser.add_argument("--samples",type=int,default=256)
    parser.add_argument("--erase",choices=("qrom","anf"),default="qrom")
    a=parser.parse_args()
    rows=[]
    for p in a.pads:
        try:
            circuit,bound=build(p,a.erase,a.erase)
            stats=resource_counts(circuit)
            r=dict(pad=p,erase=a.erase,wrap=bound,
                   native=stats,
                   beats_verified_v9_at_toffoli_floor=stats["lower_bound_qc"]<PRIOR_FULL,
                   reaches_20x_by_toffoli_floor=stats["lower_bound_qc"]<REF/20)
            if r["beats_verified_v9_at_toffoli_floor"] and a.samples:
                from .hybrid_verify import boundary_inputs,verify_samples
                r["coherent"]=verify_samples(circuit,boundary_inputs(random_count=a.samples))
        except Exception as e:
            r=dict(pad=p,erase=a.erase,error=type(e).__name__+":"+str(e)[:450])
        rows.append(r)
        print(json.dumps(r),flush=True)
    path=Path('runs/early_erase_and_shortphase_native.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,indent=2)+"\n")
if __name__=="__main__":main()
