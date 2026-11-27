"""Original parity-compressed half-square direct phase oracle.

For 13-bit r, define S=T_13(r,r), where T_d is the exact sum of binary
partial-product columns at/above d (NOT floor of a rounded product).
Pair the off-diagonal terms (i,j),(j,i): each contributes an even
integer. Every surviving diagonal i=j has 2i>=14, hence contributes
2^(2i-13), also even. Thus S is ALWAYS EVEN and H=S/2 is an exact
12-bit word.

Crucially, for every nonnegative coefficient g:
       T_8(g,S) = T_7(g,H).
No arithmetic approximation or target phase is changed.

The original v9 circuit unnecessarily allocated 13 square bits and
multiplied gamma by the permanently zero least significant square bit.
This implementation directly constructs H as a carry-safe symmetric
half-square in 12 qubits and uses an exact 7-column-cut product,
with the original first-row carry-free shortcut. The initial ANF ROM,
phase coefficients, input normalization and measured cleanup are
unchanged. Independent development, no competing Kaggle circuit.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from . import benchmark1_cost_cooptim as cc
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from .benchmark1_pruned_anf_rom import binary_pruned_rom
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_first_row_specialization import mul_first_row
from .hybrid_kickmix import resource_counts

CUTS=(10,13,12)
INCUMBENT=237684.49222389946
REFERENCE=2076880.0

def half_square_certificate():
    r=np.arange(8192,dtype=np.int64)
    square=truncated_product(r,r,13)
    if np.any(square&1):
        raise RuntimeError("Odd truncated symmetric square! identity failed")
    half=square>>1
    assert half.max()<1<<12
    # Check ALL possible 7-bit nonnegative coefficients and all residuals.
    wrong=0
    for g in range(128):
        before=truncated_product(np.int64(g),square,8)
        after=truncated_product(np.int64(g),half,7)
        wrong+=int(np.count_nonzero(before!=after))
    return dict(passes=wrong==0,checked_residuals=8192,
        checked_coefficient_values=128,
        checked_product_pairs=8192*128,
        mismatched_products=wrong,max_half_square=int(half.max()),
        identity="T8(g,T13(r,r))=T7(g,T13(r,r)/2)")

def half_square_rows(builder,x,dest,inverse=False):
    if len(x)!=13 or len(dest)!=12:raise ValueError("Expected 13bit input/12bit half square")
    d=13
    allrows=[]
    for i in range(13):
        terms=[]
        if 2*i>=d:
            terms.append((2*i-d-1,None))
        for j in range(i+1,13):
            if i+j>=d:
                terms.append((i+j-d,j))
        if terms:
            start=min(p for p,_ in terms)
            end=max(p for p,_ in terms)
            if start<0 or end+1>=len(dest):
                raise RuntimeError("Missing carry or invalid half-square bit")
            allrows.append((i,start,end,terms))
    if not allrows:raise RuntimeError("No live rows")
    special=allrows[0]
    order=allrows if not inverse else list(reversed(allrows))
    for i,start,end,terms in order:
        row=builder.alloc_qubits(end-start+2)
        builder.reset(row)
        for pos,j in terms:
            if j is None:builder.cx(x[i],row[pos-start])
            else:builder.init_and(x[i],x[j],row[pos-start])
        if i==special[0]:
            for k in range(len(row)-1):
                builder.cx(row[k],dest[start+k])
        elif inverse:
            builder.isub(row,target=dest[start:start+len(row)])
        else:
            builder.iadd(row,target=dest[start:start+len(row)])
        for pos,j in reversed(terms):
            if j is None:builder.cx(x[i],row[pos-start])
            else:builder.del_and(x[i],x[j],row[pos-start])
        builder.free(row)

def build():
    if not half_square_certificate()["passes"]:raise RuntimeError("Half-square proof failure")
    if not cc.full_screen(CUTS)["passes"]:raise RuntimeError("P7 phase numerical certificate failure")
    import kickmix as km
    mat,packed=optimized_coefficients("affine")
    assert proxy(packed)["top_interaction_monomials"]==0
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    tA=[int(a) for a in A]
    tBC=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(b,x)
    word=b.alloc_qubits(wa+wb+wc)
    b.reset(word)
    binary_pruned_rom(b,x[cc.LO:],word,packed)
    a=word[:wa];beta=word[wa:wa+wb];gamma=word[wa+wb:]
    with b.mark("intercept_phase"):
        phase_with_pad(b,a,shift,cc.Q,+1,9)
        b.del_lookup(table=tA,address=x[cc.LO:],target=a)
        b.free(a)
    dl,ds,dq=CUTS
    with b.mark("linear_phase"):
        drop=dl-cc.BZ
        out=b.alloc_qubits(wb+cc.LO-drop);b.reset(out)
        mul_first_row(b,beta,x[:cc.LO],out,drop)
        phase_with_pad(b,out,shift,cc.Q+cc.LO-dl,-1,5)
        mul_first_row(b,beta,x[:cc.LO],out,drop,inverse=True)
        b.free(out)
    with b.mark("half_square_quadratic_phase"):
        half=b.alloc_qubits(12);b.reset(half)
        half_square_rows(b,x[:cc.LO],half)
        drop=dq-cc.CZ-1
        out=b.alloc_qubits(wc+len(half)-drop);b.reset(out)
        mul_first_row(b,half,gamma,out,drop)
        phase_with_pad(b,out,shift,cc.Q+2*cc.LO-ds-dq,+1,5)
        mul_first_row(b,half,gamma,out,drop,inverse=True)
        b.free(out)
        half_square_rows(b,x[:cc.LO],half,inverse=True)
        b.free(half)
    b.del_lookup(table=tBC,address=x[cc.LO:],target=word[wa:])
    b.free(word[wa:])
    denormalize_binary_block(b,x,flags,shift)
    return b.finish_circuit()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,default=Path("out/half_square_v10.kmb"))
    ap.add_argument("--samples",type=int,default=256)
    args=ap.parse_args()
    p=half_square_certificate()
    c=build()
    s=resource_counts(c)
    r=dict(math=p,native=s,
      can_beat_full_v9_at_toffoli_floor=s["lower_bound_qc"]<INCUMBENT,
      can_reach_20x_at_toffoli_floor=s["lower_bound_qc"]<REFERENCE/20)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    r["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if args.samples and r["can_beat_full_v9_at_toffoli_floor"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        r["coherent"]=verify_samples(c,boundary_inputs(random_count=args.samples))
        if r["coherent"]["failures"]:raise RuntimeError("Half-square circuit coherent FAIL: "+str(r))
    print(json.dumps(r,indent=2))

if __name__=="__main__":main()
