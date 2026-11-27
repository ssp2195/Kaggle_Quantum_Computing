"""Original coefficient-output-bit borrowing in a coherent multioutput ANF ROM.

For all 384 legal normalized addresses j=128,...,511, the integer
intercept coefficient A has:
   A[25] XOR A[26] = 1.
This is fully checked against the EXACT existing coefficient matrix.

Instead of computing both bits of A in the main multioutput GF(2) ROM,
we temporarily use the allocated A[25] (initially clean zero) as a
recyclable AND-prefix scratch qubit. The decoder computes every
coefficient bit except A[25]. All AND measurements and Clifford
feedforward leave the scratch clean. Finally, reconstruct A[25]
with X(A[25]) and CX(A[26],A[25]). The full 50-bit coefficient
register now equals the original exactly on every legal input.

This eliminates one clean ancillary bit from the recursion peak with
ZERO extra Toffoli-like gates and ZERO extra fine rotations. QROM
unlookup and all subsequent arithmetic/phase are unchanged. The
conjugation holds coherently in superpositions of all legal addresses.

The 128 impossible table rows need not satisfy the complement relation.
No competitor source or circuit is imported.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from . import benchmark1_pruned_anf_rom as rom
from . import benchmark1_dontcare_boolean_rom as dc
from .benchmark1_first_row_specialization import build as v9_build
from .hybrid_kickmix import resource_counts

BASE=237684.49222389946
REF=2076880.
BORROWED=25
COMPLEMENT=26

def complement_certificate():
    mat,_=dc.optimized_coefficients("affine")
    A=mat[:,0]
    vals=[int(q) for q in A]
    invalid=[j for j in range(128,512)
             if (((vals[j]>>BORROWED)^(vals[j]>>COMPLEMENT))&1)!=1]
    return dict(passes=not invalid,checked_rows=384,
       borrowed_output_bit=BORROWED,source_output_bit=COMPLEMENT,
       failed_rows=invalid,
       relation="A[25] XOR A[26] is constantly 1 over allowed addresses")

def decode_borrowing(b,address,target,table):
    import kickmix as km
    if len(target)!=50 or len(address)!=9:
        raise ValueError("Only independently certified p7 50-bit coefficient table")
    coeff=rom.transform(table,len(address))
    if any(int(w).bit_length()>len(target) for w in coeff):
        raise RuntimeError("ANF output width wrong")
    required={0}
    for m,w in enumerate(coeff):
        if not (int(w)&~(1<<BORROWED)):continue
        t=m
        while True:
            required.add(t)
            if not t:break
            t^=1<<(t.bit_length()-1)
    freepool=[target[BORROWED]]
    borrowed_uses=[0]
    def emit(control,value):
        v=int(value)&~(1<<BORROWED)
        i=0
        while v:
            if v&1:
                if control is None:b.x(target[i])
                else:b.cx(control,target[i])
            i+=1
            v>>=1
    emit(None,coeff[0])
    def visit(start,parent,mask):
        for j in range(start,len(address)):
            child=mask|(1<<j)
            if child not in required:continue
            if parent is None:
                q=address[j];kind="input"
            elif freepool:
                q=freepool.pop()
                b.init_and(parent,address[j],q)
                borrowed_uses[0]+=1
                kind="borrowed"
            else:
                q=b.alloc_qubits(1)[0]
                b.init_and(parent,address[j],q)
                kind="owned"
            emit(q,coeff[child])
            visit(j+1,q,child)
            if kind in ("owned","borrowed"):
                b.del_and(parent,address[j],q)
                if kind=="owned":b.free(q)
                else:freepool.append(q)
    visit(0,None,0)
    if len(freepool)!=1 or freepool[0]!=target[BORROWED]:
        raise RuntimeError("Borrowed clean output qubit not restored")
    b.x(target[BORROWED])
    b.cx(target[COMPLEMENT],target[BORROWED])
    return dict(borrowed_ands=borrowed_uses[0],
        all_coefficient_bits_reconstructed=True)

def build():
    if not complement_certificate()["passes"]:
        raise RuntimeError("Borrow relation not valid")
    orig=rom.binary_pruned_rom
    rom.binary_pruned_rom=decode_borrowing
    try:c=v9_build(True,True)
    finally:rom.binary_pruned_rom=orig
    return c

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/output_borrowed_anf.kmb"))
    p.add_argument("--samples",type=int,default=256)
    a=p.parse_args()
    report=dict(math=complement_certificate())
    c=build();cost=resource_counts(c)
    report["native"]=cost
    report["can_beat_current_at_toffoli_floor"]=cost["lower_bound_qc"]<BASE
    report["can_reach_20x_at_toffoli_floor"]=cost["lower_bound_qc"]<REF/20
    a.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(a.output,format="kmb")
    report["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
    if a.samples and report["can_beat_current_at_toffoli_floor"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        report["sample"]=verify_samples(c,boundary_inputs(random_count=a.samples))
        if report["sample"]["failures"]:
            raise RuntimeError("Coherent scratch reuse failed: "+str(report["sample"]))
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
