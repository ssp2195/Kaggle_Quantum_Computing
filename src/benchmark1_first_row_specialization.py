"""Original carry-free first-row specialization for reversible cut multipliers.

In the verified oracle, each product accumulator starts at |0>.
Its first selected source row can therefore be assigned by CCX into
known-zero output bits, with no ripple carry. On exact uncomputation,
the same row is last, after all other additions have been reversed,
so the same CCX gates return it to zero. This is exact and preserves
all intermediate phase values. Use the longest selected row to
maximize avoided carry work.

Also test an initial-row optimization of the symmetric-square stage:
the first clean AND row is copied with Clifford CX, rather than
ripple-added. The remaining row order is arbitrary, as integer
addition commutes and the output has guaranteed sufficient width.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from .benchmark1_phase_pad9_symmetric import build as v8_build
from . import benchmark1_pruned_anf_rom as rom
from . import benchmark1_symmetric_square_rows_v7 as rows
from .hybrid_kickmix import resource_counts

VERIFIED_V8=241424.49222389946
REFERENCE=2076880

def mul_first_row(builder,a,b,dest,drop,inverse=False):
    na,nb=len(a),len(b)
    if len(dest)!=na+nb-drop:raise ValueError("Product output width mismatch")
    js=[j for j in range(nb) if max(0,drop-j)<na]
    if not js:return
    special=max(js,key=lambda j:(na-max(0,drop-j),-j))
    order=[special]+[j for j in js if j!=special]
    if inverse:order=list(reversed(order))
    for j in order:
        first=max(0,drop-j)
        offset=j+first-drop
        source=a[first:]
        if j==special:
            # The destination is exactly zero when the first addend
            # arrives, or exactly this row when the last one is erased.
            for i,q in enumerate(source):
                builder.ccx(q,b[j],dest[offset+i])
        elif inverse:
            builder.isub(source,target=dest[offset:offset+len(source)+1],control=b[j])
        else:
            builder.iadd(source,target=dest[offset:offset+len(source)+1],control=b[j])

def square_first_row(builder,x,dest,drop,inverse=False):
    n=len(x)
    if len(dest)!=2*n-drop:raise ValueError("Square output width mismatch")
    allrows=[]
    for i in range(n):
        terms=[]
        if 2*i>=drop:terms.append((2*i-drop,None))
        for j in range(i+1,n):
            if i+j>=drop:terms.append((i+j+1-drop,j))
        if terms:
            start=min(z for z,_ in terms)
            end=max(z for z,_ in terms)
            allrows.append((i,start,end,terms))
    special=max(allrows,key=lambda row:(row[2]-row[1]+1,-row[0]))
    order=[special]+[x for x in allrows if x!=special]
    if inverse:order=list(reversed(order))
    for i,start,end,terms in order:
        row=builder.alloc_qubits(end-start+2)
        builder.reset(row)
        for pos,j in terms:
            if j is None:builder.cx(x[i],row[pos-start])
            else:builder.init_and(x[i],x[j],row[pos-start])
        if (i,start,end,terms)==special:
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

def build(linear=False,square=False):
    previous_l=rom.mul_controlled
    previous_s=rows.square_symmetric_rows
    rom.mul_controlled=mul_first_row if linear else previous_l
    rows.square_symmetric_rows=square_first_row if square else previous_s
    try: c=v8_build()
    finally:
        rom.mul_controlled=previous_l
        rows.square_symmetric_rows=previous_s
    return c

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=256)
    ap.add_argument("--outputdir",type=Path,default=Path("out"))
    args=ap.parse_args()
    reports=[]
    for first,sq in [(False,False),(True,False),(False,True),(True,True)]:
        try:
            c=build(first,sq)
            s=resource_counts(c)
            result=dict(linear_first=first,square_first=sq,native=s,
              beats_verified_v8_at_floor=s["lower_bound_qc"]<VERIFIED_V8,
              potentially_twenty_x=s["lower_bound_qc"]<REFERENCE/20)
            if first or sq:
                args.outputdir.mkdir(parents=True,exist_ok=True)
                p=args.outputdir/(f"firstrow_{int(first)}_{int(sq)}.kmb")
                c.to_file(p,format="kmb")
                result["sha256"]=hashlib.sha256(p.read_bytes()).hexdigest()
            if (first or sq) and result["beats_verified_v8_at_floor"] and args.samples:
                from .hybrid_verify import boundary_inputs,verify_samples
                result["coherent"]=verify_samples(c,boundary_inputs(random_count=args.samples))
        except Exception as e:
            result=dict(linear_first=first,square_first=sq,
                        error=type(e).__name__+":"+str(e)[:500])
        print(json.dumps(result),flush=True)
        reports.append(result)
    p=Path("runs/firstrow_shortcuts.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(reports,indent=2)+"\n")

if __name__=="__main__":main()
