"""Experimental direct-controlled-adder multiplier: no AND row register.

For each input bit b[j], directly add the live a bit-slice to the
target under that bit's control. This implements exactly the same
truncated partial-product sum as src.priority1_truncated._mul_cut.
Square operands alias; in that case use the reference implementation
until the native backend has passed an independent alias validation.
"""
from __future__ import annotations
from .priority1_truncated import _mul_cut as _reference

def mul_controlled(builder,a,b,dest,drop,inverse=False):
    na,nb=len(a),len(b)
    if len(dest)!=na+nb-drop:
        raise ValueError('Incorrect destination width')
    # Native controlled adder may not support a control equal to an addend
    # qubit; keep the established square implementation.
    if any(ai==bj for ai in a for bj in b):
        return _reference(builder,a,b,dest,drop,inverse=inverse)
    for j in (range(nb-1,-1,-1) if inverse else range(nb)):
        first=max(0,drop-j)
        if first>=na:continue
        offset=j+first-drop
        source=a[first:]
        target=dest[offset:offset+len(source)+1]
        if inverse:
            builder.isub(source,target=target,control=b[j])
        else:
            builder.iadd(source,target=target,control=b[j])


def build_sparse_controlled():
    from . import priority1_truncated
    from .benchmark1_sparse_coeff import build
    previous=priority1_truncated._mul_cut
    priority1_truncated._mul_cut=mul_controlled
    try:
        return build()
    finally:
        priority1_truncated._mul_cut=previous


def main():
    import argparse,json
    from pathlib import Path
    from .hybrid_kickmix import resource_counts
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=Path('out/controlled_coeff_invsqrt.kmb'))
    p.add_argument('--verify',action='store_true')
    p.add_argument('--random-count',type=int,default=1024)
    a=p.parse_args()
    circuit=build_sparse_controlled()
    costs=resource_counts(circuit)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(a.output,format='kmb')
    report={'resources':costs,
            'toffoli_floor':costs['lower_bound_qc'],
            'may_beat_reference':costs['lower_bound_qc']<430259.61862786254,
            'may_exceed_10x':costs['lower_bound_qc']<207688.}
    if a.verify:
        from .hybrid_verify import verify_samples,boundary_inputs
        report['coherent']=verify_samples(circuit,boundary_inputs(random_count=a.random_count))
        if report['coherent']['failures']:
            print(json.dumps(report,indent=2))
            raise SystemExit(1)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
