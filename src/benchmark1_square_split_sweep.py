"""Original exact quantum square split-parameter search.

For 13-bit r=L+2^k H and 4<=k<=7, dropped-column (not
rounded) partial-product arithmetic obeys:
  T13(r,r) = T_{13-2k}(H,H) + 2*T_{13-k}(L,H),
when 13-2k>=0, else 2^{2k-13} H^2 + 2*T_{13-k}(L,H).
Low L^2 contributes nothing because 2k-2<13.

Native builder computes H² in carry-safe rows and one strongly cut
cross-product, adds the cross output once, and reverses every operation.
No dependency on competitor implementation. Full 8192-word certificates
and exact quantum coherence checks guard every candidate.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .priority1_truncated import _mul_cut as mul_cut
from .benchmark1_symmetric_square_rows_v7 import square_symmetric_rows
from .benchmark1_phase_pad9_symmetric import build as v8_build
from . import benchmark1_pruned_anf_rom as rom
from .hybrid_kickmix import resource_counts

REFERENCE=2076880.0
V8=241424.49222389946

def math_certificate(k):
    if k not in (4,5,6,7):raise ValueError(k)
    r=np.arange(8192,dtype=np.int64);low=r&((1<<k)-1);hi=r>>k
    exact=truncated_product(r,r,13)
    dH=13-2*k
    high=truncated_product(hi,hi,dH) if dH>=0 else (hi*hi)<<(-dH)
    cross=2*truncated_product(low,hi,13-k)
    bad=int(np.count_nonzero(exact!=high+cross))
    return dict(split=k,passes=bad==0,failures=bad,
      checked=8192,max_output=int(exact.max()),
      high_word_bits=13-k,high_square_drop=dH,cross_drop=13-k)

def square_split(b, r, dest, drop, k=6, inverse=False):
    if len(r)!=13 or len(dest)!=13 or drop!=13:
        raise ValueError("This family is certified only for T13 of 13bits")
    L,H=r[:k],r[k:]
    dH=13-2*k
    high_target=dest if dH>=0 else dest[-dH:]
    high_drop=max(0,dH)
    if not inverse:
        square_symmetric_rows(b,H,high_target,high_drop)
    cross=b.alloc_qubits(k)
    b.reset(cross)
    mul_cut(b,L,H,cross,13-k)
    if inverse:b.isub(cross,target=dest[1:])
    else:b.iadd(cross,target=dest[1:])
    mul_cut(b,L,H,cross,13-k,inverse=True)
    b.free(cross)
    if inverse:
        square_symmetric_rows(b,H,high_target,high_drop,inverse=True)

def make(k):
    if not math_certificate(k)["passes"]:
        raise RuntimeError("Algebraic split certificate failed")
    original=rom._mul_cut
    def patched(b,a,a2,dest,drop,inverse=False):
        if len(a)==13 and len(a2)==13 and all(u==v for u,v in zip(a,a2)):
            return square_split(b,a,dest,drop,k,inverse)
        return original(b,a,a2,dest,drop,inverse)
    rom._mul_cut=patched
    try:c=v8_build()
    finally:rom._mul_cut=original
    return c

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--splits",nargs="+",type=int,default=[4,5,6,7])
    parser.add_argument("--samples",type=int,default=128)
    args=parser.parse_args()
    result=[]
    for k in args.splits:
        try:
            cert=math_certificate(k)
            circ=make(k);s=resource_counts(circ)
            record=dict(split=k,math=cert,native=s,
                could_beat_v8_by_toffoli_floor=s["lower_bound_qc"]<V8,
                reaches_20x_by_floor=s["lower_bound_qc"]<REFERENCE/20)
            if args.samples and record["could_beat_v8_by_toffoli_floor"]:
                from .hybrid_verify import boundary_inputs,verify_samples
                record["sample"]=verify_samples(circ,boundary_inputs(random_count=args.samples))
        except Exception as e:
            record=dict(split=k,error=type(e).__name__+":"+str(e)[:400])
        print(json.dumps(record),flush=True)
        result.append(record)
    path=Path("runs/symmetric-square-split-sweep.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+"\n")

if __name__=="__main__":main()
