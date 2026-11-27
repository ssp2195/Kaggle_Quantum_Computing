"""Replace general truncated square with algebraically symmetric square circuit.

Original square microarchitecture test: diagonal bits are copied, off-diagonal
terms are doubled once. Exact dropped-column integer arithmetic must agree
on every normalized residual; only full quantum verification can promote.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .benchmark1_zero_residual_pad import build as build_reduced_pad
from .benchmark1_symmetric_square import square_symmetric
from .hybrid_kickmix import resource_counts
from .priority1_exact import truncated_product

INCUMBENT=272300.78036771005

def mathematical_square_certificate():
    # Certify all 8192 squares, all 13 output positions and partial cuts.
    n=cc.LO;d=13
    r=np.arange(1<<n,dtype=np.int64)
    full=(r*r).astype(np.int64)
    selected=np.zeros_like(r)
    for i in range(n):
        if 2*i>=d:selected+=(((r>>i)&1)<<(2*i-d))
        for j in range(i+1,n):
            if i+j>=d:selected+=((((r>>i)&1)&((r>>j)&1))<<(i+j+1-d))
    ref=truncated_product(r,r,d)
    failures=int(np.count_nonzero(selected!=ref))
    return dict(passes=not failures,residuals=1<<n,failures=failures,
        sample_max=int(np.max(selected)),drop=d)

def build():
    cert=mathematical_square_certificate()
    if not cert["passes"]:raise ValueError(str(cert))
    old=cc._mul_cut
    def replacement(builder,a,b,dest,drop,inverse=False):
        if len(a)==len(b) and all(q==z for q,z in zip(a,b)):
            return square_symmetric(builder,a,dest,drop,inverse)
        return old(builder,a,b,dest,drop,inverse)
    cc._mul_cut=replacement
    try:
        circuit,_calls=build_reduced_pad()
    finally:
        cc._mul_cut=old
    return circuit

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=Path('out/symmetric_square_phase5.kmb'))
    ap.add_argument('--samples',type=int,default=256)
    args=ap.parse_args()
    cert=mathematical_square_certificate()
    c=build()
    stats=resource_counts(c)
    result=dict(math=cert,resources=stats,could_beat_incumbent=stats["lower_bound_qc"]<INCUMBENT)
    if result["could_beat_incumbent"]:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(args.output,format='kmb')
        result["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
        from .hybrid_verify import verify_samples,boundary_inputs
        result["coherent_sample"]=verify_samples(c,boundary_inputs(random_count=args.samples))
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
