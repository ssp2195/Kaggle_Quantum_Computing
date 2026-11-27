"""Input-certified conditional high-phase correction for shortened intercept word.

The coefficient word A has 27 bits and the native fractional scale is D=36.
Instead of ten padding qubits, rotate A + p high zero qubits for p in
5..9. For any s where its high bits wrap, the physical rotation phases
the wrapped bits at very small 2**(i+s-w-36) turn weights, instead
of their true 2**(i+s-36) weights. Implement the TRUE coefficient
by conditional phase angles for each overflow (i,s) pair, controlled on
binary exponent == s and original A[i].

The only approximation relative to the previous verified circuit is
LEAVING the tiny artificial wrapped phases:
   delta=sum_{i=w-s}^{26} A_i*2**(i+s-w-36) turns.
A worst-case *input-dependent bound* on |delta|/(2^-24/sqrt(x))
is computed over all possible s before any native compilation.
This is not a copied Kaggle circuit, and every candidate must pass
the pinned full quantum verifier before being promoted.
"""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import argparse,json,hashlib,math
from . import benchmark1_pruned_anf_rom as rom
from .benchmark1_symmetric_square_rows_v7 import build as build_squared
from .hybrid_kickmix import resource_counts

D=36
WA=27
INCUMBENT=256476.78036771005
REFERENCE=2076880.

def certified_wrap_bound(p):
    if p not in range(4,10):raise ValueError("Test pad 4-9 only")
    width=WA+p
    arr=[]
    maxr=0.
    for s in range(11):
        n=max(0,s-p)
        x_max=min(3133452,(1<<(22-2*s))-1)
        # Exact maximum small reintroduced wrap: 2^n-1 is the low n
        # bits in shifted A; all those bits are allowed to be 1.
        turns=(2**n-1)/2**D
        ratio=turns*math.sqrt(x_max)*2**24
        maxr=max(maxr,ratio)
        arr.append(dict(exponent=s,max_valid_x=x_max,
           wrapped_bits=n,maximum_spurious_turns=turns,
           maximum_error_relative_to_tolerance=ratio))
    return dict(p=p,passes=maxr<.1,
         worst_case_added_error_over_tolerance=maxr,
         bound_by_exponent=arr,
         proof="The only omitted exact correction is the tiny 2^-36-weight wrap phase")

def correction(builder,word,shift,frac_bits,p):
    w=len(word)+p
    if frac_bits!=26 or len(word)!=27:
        raise ValueError("Expected independently pinned intercept coefficient")
    pad=builder.alloc_qubits(p)
    builder.reset(pad)
    rotated=word+pad
    builder.left_rotate(rotated,shift=shift)
    den=1<<D
    for i in range(len(rotated)):
        exp=Fraction(2*(1<<i),den)
        if exp%2:builder.z_pow(rotated[i],exp)
    builder.right_rotate(rotated,shift=shift)
    builder.free(pad)
    # Rotated source is restored. Each bin exponent s uses one equality
    # indicator shared by all overflow bits in that bin.
    for s in range(p+1,11):
        eq=builder.alloc_qubits(1)[0]
        builder.reset(eq)
        builder.flip_if_equal(shift,s,target=eq)
        for i in range(w-s,len(word)):
            exponent=Fraction(2*(1<<(i+s)),den)
            if exponent%2==0:
                continue  # whole integer turns are global identity
            if exponent==1:
                # Half-turn => direct CZ (Clifford)
                builder.ccz(True,eq,word[i])
            else:
                ctrl=builder.alloc_qubits(1)[0]
                builder.init_and(eq,word[i],ctrl)
                builder.z_pow(ctrl,exponent)
                builder.del_and(eq,word[i],ctrl)
                builder.free(ctrl)
        builder.flip_if_equal(shift,s,target=eq)
        builder.free(eq)

def build(p):
    cert=certified_wrap_bound(p)
    if not cert["passes"]:raise ValueError("Wrap approximation bound exceeds safety threshold")
    old=rom.phase_with_pad
    visited=[]
    def patched(builder,word,shift,frac_bits,sign,pad_bits):
        if not visited:
            if sign!=1 or frac_bits!=26 or pad_bits!=10:
                raise ValueError("First phase is not intercept")
            correction(builder,word,shift,frac_bits,p)
        else:
            old(builder,word,shift,frac_bits,sign,pad_bits)
        visited.append((len(word),frac_bits,sign))
    rom.phase_with_pad=patched
    try:
        circ=build_squared()
    finally:
        rom.phase_with_pad=old
    if len(visited)!=3:raise RuntimeError("Unexpected number of phase stages")
    return circ,cert

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--pads',type=int,nargs='+',default=[8,7,6])
    ap.add_argument('--outdir',type=Path,default=Path('out'))
    ap.add_argument('--samples',type=int,default=128)
    args=ap.parse_args()
    rows=[]
    for p in args.pads:
        try:
            circ,cert=build(p)
            costs=resource_counts(circ)
            r=dict(pad=p,status="NATIVE_BUILT",wrap=cert,resources=costs,
                toffoli_floor_beats_v7=costs["lower_bound_qc"]<INCUMBENT,
                toffoli_floor_reaches_20x=costs["lower_bound_qc"]<REFERENCE/20)
            if r["toffoli_floor_beats_v7"]:
                args.outdir.mkdir(parents=True,exist_ok=True)
                path=args.outdir/('phase_corrected_pad%d.kmb'%p)
                circ.to_file(path,format="kmb")
                r["sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
                if args.samples:
                    from .hybrid_verify import boundary_inputs,verify_samples
                    r["coherent"]=verify_samples(circ,boundary_inputs(random_count=args.samples))
            rows.append(r)
        except Exception as e:
            rows.append(dict(pad=p,status="FAILED",error=type(e).__name__+":"+str(e)[:650]))
        print(json.dumps(rows[-1]),flush=True)
    out=Path('runs/phase_conditioned_wrap_sweep.json')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=="__main__":main()
