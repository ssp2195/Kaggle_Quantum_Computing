"""Independent don't-care-optimized reversible coefficient Boolean network.

The normalized j=x_norm>>13 is ALWAYS in [128,511]. Hence table entries
for j=0..127 are unreachable don't-cares. Rather than padding them with
zeros, independently choose those 128 multi-output table values to
reduce nonlinear ANF monomial structure. The valid 384 entries (and all
output phases) remain bit-for-bit unchanged.

For each 7-bit low monomial m, let a,b,c be ANF coefficients of the
three reachable high-block values f01, f10, f11. A free coefficient v
for f00 changes four ANF blocks to:
    [v, v^a, v^b, v^a^b^c].
Choosing v=a^b^c kills the highest-degree h7*h8 block exactly.
We also search other candidates under a weighted shared-AND cost proxy.
The inverse GF(2) Möbius transform recovers an exact invalid table.

No outside Kaggle participant's code, coefficient tables, or circuit
is read. Every generated candidate is tested by the native oracle
verifier before any promotion.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import json,hashlib,argparse
import numpy as np
from . import benchmark1_cost_cooptim as cc
from . import benchmark1_pruned_anf_rom as rom
from .hybrid_kickmix import resource_counts

VALID_START=1<<cc.P
TARGET=2076880.0/20
INCUMBENT=268172.78036771005
CUTS=(10,13,12)


def extend_packed(table, policy):
    if policy not in ("zero","affine","weighted","majority"):
        raise ValueError(policy)
    source=list(map(int,table))
    if len(source)!=512 or any(v!=0 for v in source[:128]):
        raise ValueError("Expected 512 rows with first 128 unreachable zeros")
    if policy=="zero":return source
    # High address blocks: f00 invalid, f01 j=128..255,
    # f10 j=256..383, f11 j=384..511.
    A=rom.transform(source[128:256],7)
    B=rom.transform(source[256:384],7)
    C=rom.transform(source[384:512],7)
    low=[]
    for m,(a,b,c) in enumerate(zip(A,B,C)):
        fourth=a^b^c
        if policy=="affine":v=fourth
        else:
            candidates=set((0,a,b,fourth))
            if policy=="majority":
                # Minimize number of nonzero Boolean coefficient vectors.
                score=lambda v:(bool(v)+bool(v^a)+bool(v^b)+bool(v^fourth),bool(v^fourth))
            else:
                # The highest-order h7h8 interaction needs an extra shared
                # AND, so prioritize removing it.
                score=lambda v:(.35*bool(v)+1.*bool(v^a)+
                                   1.*bool(v^b)+2.25*bool(v^fourth),
                                 bool(v^fourth),v.bit_count())
            v=min(candidates,key=score)
        low.append(v)
    # ANF is its own inverse (mobius transform on the binary cube).
    values=rom.transform(low,7)
    return values+source[128:]


def optimized_coefficients(policy):
    orig=cc.coefficients(CUTS)
    if orig.shape!=(512,3):raise ValueError("Unexpected coefficient shape")
    matrix=np.array(orig,copy=True)
    A,B,C=matrix[:,0],matrix[:,1]>>cc.BZ,matrix[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    ta=[int(v) for v in A]
    tbc=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    table=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
             for a,b,c in zip(A,B,C)]
    new_table=extend_packed(table,policy)
    allmask=(1<<(wa+wb+wc))-1
    for i in range(128):
        packed=int(new_table[i])&allmask
        a=packed&((1<<wa)-1)
        b=(packed>>wa)&((1<<wb)-1)
        c=(packed>>(wa+wb))&((1<<wc)-1)
        matrix[i]=[a,b<<cc.BZ,c<<cc.CZ]
    if not np.array_equal(matrix[128:],orig[128:]):
        raise RuntimeError("Reachable coefficient row changed!")
    # Check all 512 packed entries reconstruct intended extension.
    AA,BB,CC=(matrix[:,0],matrix[:,1]>>cc.BZ,matrix[:,2]>>cc.CZ)
    observed=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb)) for a,b,c in zip(AA,BB,CC)]
    if observed!=new_table:
        raise RuntimeError("Table reconstruction disagrees with GF2 padding")
    matrix.setflags(write=False)
    return matrix,new_table

def proxy(table):
    coeff=rom.transform(table,9)
    active={m for m,c in enumerate(coeff) if c!=0}
    prefix=set()
    for mask in active:
        while mask:
            prefix.add(mask)
            mask^=1<<(mask.bit_length()-1)
    return dict(active_monomials=len(active),active_nontrivial_prefixes=
               sum(1 for m in prefix if m.bit_count()>1),
        top_interaction_monomials=sum(1 for m in active if (m&384)==384))

def build(policy):
    org=cc.coefficients
    table,new_packed=optimized_coefficients(policy)
    cc.coefficients=lambda cuts=CUTS:table
    try:
        circuit=rom.build(init="anf",erase_a="qrom",erase_bc="qrom")
    finally:
        cc.coefficients=org
    return circuit,new_packed

def main():
    rows=[]
    for policy in ("zero","affine","weighted","majority"):
        try:
            table,t=optimized_coefficients(policy)
            v=proxy(t)
            circuit,p=build(policy)
            s=resource_counts(circuit)
            out=dict(policy=policy,boolean_network=v,native=s,
                 could_beat_incumbent_floor=s["lower_bound_qc"]<INCUMBENT,
                 could_reach_20x_floor=s["lower_bound_qc"]<TARGET)
            if s["lower_bound_qc"]<INCUMBENT and policy!="zero":
                from .hybrid_verify import verify_samples,boundary_inputs
                out["coherent"]=verify_samples(circuit,boundary_inputs(random_count=128))
        except Exception as e:
            out=dict(policy=policy,error=type(e).__name__+":"+str(e)[:300])
        rows.append(out);print(json.dumps(out),flush=True)
    path=Path("runs/dontcare_boolean_rom.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,indent=2)+"\n")
if __name__=="__main__":main()
