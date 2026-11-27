"""Original carry-safe row-compressed symmetric square, not per-term increments.

For n-bit x, x² = sum_i x_i*2^(2i) +
    sum_{i<j} (x_i AND x_j)*2^(i+j+1).
The exact *partial product* truncation at d discards ordered terms
whose i+j<d, including both halves of each cross term.

Group all terms sharing i into a single contiguous clean row, then
add that row once (one carry pass). Compared with naive symmetric
squares that increment for EACH bit product, this saves repeated
carries. Rows are processed by increasing i, so a single extra
high row bit captures all possible overflow at each stage.

Each temporary AND is X-basis measured and corrected by KickMix; output
is uncomputed by reverse row subtraction. No competing code is reused.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts
from . import benchmark1_pruned_anf_rom as native
from .benchmark1_dontcare_boolean_rom import build as build_v7
from .benchmark1_best_verified_v7 import VERIFIED_QC

def square_symmetric_rows(builder,x,dest,drop,inverse=False):
    n=len(x)
    if not 0<=drop<2*n or len(dest)!=2*n-drop:
        raise ValueError("Square output width mismatch")
    allrows=[]
    for i in range(n):
        contrib=[]
        if 2*i>=drop:
            contrib.append((2*i-drop,None))
        for j in range(i+1,n):
            if i+j>=drop:
                contrib.append((i+j+1-drop,j))
        if contrib:
            start=min(p for p,j in contrib)
            end=max(p for p,j in contrib)
            if end+1>=len(dest):
                raise ValueError("No carry-safe output bit")
            allrows.append((i,start,end,contrib))
    for i,start,end,contrib in (reversed(allrows) if inverse else allrows):
        length=end-start+2
        if start+length>len(dest):
            raise ValueError("Missing single carry bit in square output")
        row=builder.alloc_qubits(length)
        builder.reset(row)
        for pos,j in contrib:
            if j is None:
                builder.cx(x[i],row[pos-start])
            else:
                builder.init_and(x[i],x[j],row[pos-start])
        if inverse:
            builder.isub(row,target=dest[start:start+length])
        else:
            builder.iadd(row,target=dest[start:start+length])
        for pos,j in reversed(contrib):
            if j is None:
                builder.cx(x[i],row[pos-start])
            else:
                builder.del_and(x[i],x[j],row[pos-start])
        builder.free(row)

def exact_math_check():
    n=13;d=13
    x=np.arange(1<<n,dtype=np.int64)
    rows=[]
    for i in range(n):
        terms=[]
        if 2*i>=d: terms.append((2*i-d,i,None))
        for j in range(i+1,n):
            if i+j>=d:terms.append((i+j+1-d,i,j))
        if terms:rows.append(terms)
    val=np.zeros_like(x)
    for terms in rows:
        row=np.zeros_like(x)
        for pos,i,j in terms:
            term=(x>>i)&1 if j is None else ((x>>i)&1)&((x>>j)&1)
            row+=term<<pos
        val+=row
    ref=truncated_product(x,x,d)
    bad=int(np.count_nonzero(val!=ref))
    return dict(passes=bad==0,checked=1<<n,failures=bad,
                dropped_partial_product_columns=d,
                rows=len(rows),full_square_output_width=2*n-d)

def build():
    if not exact_math_check()["passes"]:raise RuntimeError("Symmetric square integer model failed")
    old=native._mul_cut
    def replacement(builder,a,b,dest,drop,inverse=False):
        if len(a)==len(b) and all(x==y for x,y in zip(a,b)):
            return square_symmetric_rows(builder,a,dest,drop,inverse)
        return old(builder,a,b,dest,drop,inverse)
    native._mul_cut=replacement
    try:
        c,_=build_v7("affine")
    finally:
        native._mul_cut=old
    return c

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/symmetric_rows_v7.kmb"))
    p.add_argument("--samples",type=int,default=256)
    args=p.parse_args()
    math=exact_math_check()
    c=build()
    stats=resource_counts(c)
    report=dict(math=math,native=stats,
                toffoli_floor_beats_v7=stats["lower_bound_qc"]<VERIFIED_QC,
                twenty_x_toffoli_floor=2076880/20)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    report["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if args.samples and report["toffoli_floor_beats_v7"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        report["coherent"]=verify_samples(c,boundary_inputs(random_count=args.samples))
        if report["coherent"]["failures"]:
            raise RuntimeError("Coherent sampled phase failed: "+str(report))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
