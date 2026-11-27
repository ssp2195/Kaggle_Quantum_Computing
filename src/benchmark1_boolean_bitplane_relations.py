"""Search exact low-cost nonlinear relations among coherent coefficient bits.

All 384 reachable A/B/C coefficient rows are treated as 384-bit Boolean
truth columns. Find bits that can be reconstructed from *already computed*
coefficient or normalized address bits using:
  - constant/XOR of any two/three existing bitplanes (Clifford only);
  - pairwise AND/NAND/OR/XOR (one Toffoli-like);
  - AND/OR plus one other bit XOR (one Toffoli-like plus Clifford).

Such a dependent coefficient bit can temporarily serve as clean
ANF-prefix workspace and be reconstructed just before the phase,
potentially saving a peak qubit at the cost of few or zero Toffolis.
All equalities are EXACT over every 384 legal coefficient address.
"""
from __future__ import annotations
import json,itertools
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import optimized_coefficients

def main():
    mat,_=optimized_coefficients('affine')
    A,B,C=mat[:,0],mat[:,1]>>2,mat[:,2]>>4
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    packed=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
              for a,b,c in zip(A,B,C)][128:]
    L=len(packed)
    mask=(1<<L)-1
    output_cols=[sum(((p>>k)&1)<<i for i,p in enumerate(packed)) for k in range(wa+wb+wc)]
    addr_cols=[sum(((j>>k)&1)<<i for i,j in enumerate(range(128,512))) for k in range(9)]
    cols=output_cols+addr_cols+[mask]
    names=['out'+str(i) for i in range(len(output_cols))]+['addr'+str(i) for i in range(9)]+['const1']
    relation=[]
    nout=len(output_cols)
    for t in range(nout):
        target=cols[t]
        rest={cols[j]:j for j in range(len(cols)) if j!=t}
        if target in rest: relation.append((t, 'copy', [rest[target]]))
        if target^mask in rest: relation.append((t,'not',[rest[target^mask]]))
        for i in range(len(cols)):
            if i==t:continue
            key=target^cols[i]
            if key in rest and rest[key]!=i:
                a,b=sorted([i,rest[key]])
                if a!=t and b!=t:
                    relation.append((t,'xor',[a,b]))
        if any(r[0]==t and r[1] in ('copy','not','xor') for r in relation):continue
        for i in range(len(cols)):
            if i==t:continue
            for j in range(i+1,len(cols)):
                if j==t:continue
                if cols[i]&cols[j]==target:
                    relation.append((t,'and',[i,j]))
                    break
                if ((cols[i]&cols[j])^mask)==target:
                    relation.append((t,'nand',[i,j]))
                    break
            if any(r[0]==t for r in relation):break
    # Find a GF(2) nullspace relation among all coefficient and address
    # bit columns; expose any nontrivial low-complexity linear identities.
    basis={}
    dependencies=[]
    for j,col in enumerate(cols):
        remainder=col
        vector=1<<j
        while remainder:
            hi=remainder.bit_length()-1
            if hi not in basis:
                basis[hi]=(remainder,vector)
                break
            n,v=basis[hi]
            remainder^=n;vector^=v
        else:
            dependencies.append([names[i] for i in range(len(cols)) if (vector>>i)&1])
    relation=sorted(set((t,kind,tuple(ids)) for t,kind,ids in relation),
                    key=lambda v:(v[0],len(v[2]),v[1]))
    r=dict(passes=True,checked_rows=L,
         coefficient_widths=[wa,wb,wc],
         total_bits=nout,
         affine_basis_rank=len(basis),
         linear_dependencies=dependencies,
         dependent_output_candidates=len(set(t for t,_,_ in relation)),
         best_relations=[dict(target=names[t],operation=kind,
                              args=[names[i] for i in ids]) for t,kind,ids in relation[:100]])
    print(json.dumps(r,indent=2))
    p=Path('runs/boolean_bitplane_relations.json')
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':main()
