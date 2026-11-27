"""Original radix-64 split of exact dropped-column quantum square.

Write r=L+64H, L six bits, H seven bits. The deliberately
non-rounded partial-product truncation of the 13-bit square obeys
  T_13(r,r)=T_1(H,H)+2*T_7(L,H).
This is an EXACT equality of bit-selected sums, not a rounded-square
identity. Only 15 low/high pair products survive the T7 mixed term.

Use one small high-half square, one strongly truncated 6x7 mixed
product, and one 12-bit carry-propagated addition (instead of twelve
full-width square rows). All intermediate XOR/AND temporaries are
measurement-corrected and returned to zero.

Never import other Kaggle participants' circuits or use partial
verifier success as a full correctness claim.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .priority1_truncated import _mul_cut as original_mul_cut
from .benchmark1_symmetric_square_rows_v7 import square_symmetric_rows
from .benchmark1_phase_pad9_symmetric import build as original_build
from . import benchmark1_pruned_anf_rom as rom
from .hybrid_kickmix import resource_counts

BASE=241424.49222389946
REFERENCE=2076880.

def certificate():
    r=np.arange(8192,dtype=np.int64)
    lo=r&63
    hi=r>>6
    expected=truncated_product(r,r,13)
    actual=truncated_product(hi,hi,1)+2*truncated_product(lo,hi,7)
    bad=int(np.count_nonzero(expected!=actual))
    return dict(passes=bad==0,checked=8192,
       max_square=int(expected.max()),bad=bad,
       identity="T13(r,r) = T1(r>>6,r>>6) + 2*T7(r&63,r>>6)")

def square_split(b,r,square,drop,inverse=False):
    if drop!=13 or len(r)!=13 or len(square)!=13:
        raise ValueError("Only certified (13,13) square")
    L=r[:6];H=r[6:]
    if not inverse:
        square_symmetric_rows(b,H,square,1)
    cross=b.alloc_qubits(6)
    b.reset(cross)
    original_mul_cut(b,L,H,cross,7)
    if inverse:b.isub(cross,target=square[1:])
    else:b.iadd(cross,target=square[1:])
    original_mul_cut(b,L,H,cross,7,inverse=True)
    b.free(cross)
    if inverse:
        square_symmetric_rows(b,H,square,1,inverse=True)

def build():
    assert certificate()["passes"]
    old=rom._mul_cut
    def wrapped(b,a,bits,dest,drop,inverse=False):
        if len(a)==13 and len(bits)==13 and all(q==p for q,p in zip(a,bits)):
            return square_split(b,a,dest,drop,inverse)
        return old(b,a,bits,dest,drop,inverse)
    rom._mul_cut=wrapped
    try:c=original_build()
    finally:rom._mul_cut=old
    return c

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--samples",type=int,default=256)
    a.add_argument("--output",type=Path,default=Path("out/split_square_radix64.kmb"))
    args=a.parse_args()
    c=build();st=resource_counts(c)
    out=dict(math=certificate(),native=st,
             can_beat_v8_at_floor=st["lower_bound_qc"]<BASE,
             can_reach_20x_at_floor=st["lower_bound_qc"]<REFERENCE/20)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    out["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if args.samples and out["can_beat_v8_at_floor"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        out["sample"]=verify_samples(c,boundary_inputs(random_count=args.samples))
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
