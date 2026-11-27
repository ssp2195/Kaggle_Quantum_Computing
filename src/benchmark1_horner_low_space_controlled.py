"""Low-space exact Horner: no explicit AND row registers.

The previous in-place Horner synthesizer allocated an entire n-bit row
for each product bit, pushing peak width above 150. We instead perform
native controlled arithmetic into the *entire* accumulator suffix, so
all borrow/carry propagation is retained exactly.

Mathematical coefficients are independently generated in the existing
p6/p7 Chebyshev modules. No published competitor source is incorporated.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .hybrid_kickmix import resource_counts
from . import benchmark1_horner_xag as hx
from . import benchmark1_centered_horner as cx
from .benchmark1_compact_normalize import normalize_compact,denormalize_compact

CURRENT=272300.78036771005
TARGET=2076880/20

def controlled_row_update(builder,source,operand,dest,drop,*,cleanup="mbuc",inverse=False):
    n,m=len(source),len(operand)
    for j in (range(m-1,-1,-1) if inverse else range(m)):
        first=max(0,drop-j)
        if first>=n:continue
        start=j+first-drop
        if start+n-first>len(dest):
            raise ValueError("Accumulator too narrow")
        if inverse:
            builder.iadd(source[first:],target=dest[start:],control=operand[j])
        else:
            builder.isub(source[first:],target=dest[start:],control=operand[j])

def signed_controlled_update(builder,source,operand,dest,drop,*,cleanup="mbuc",inverse=False):
    n,m=len(source),len(operand)
    for j in (range(m-1,-1,-1) if inverse else range(m)):
        first=max(0,drop-j)
        if first>=n:continue
        start=j+first-drop
        if start+n-first>len(dest):
            raise ValueError("Accumulator too narrow")
        positive=j==m-1
        subtract=(not positive)!=inverse
        if subtract:
            builder.isub(source[first:],target=dest[start:],control=operand[j])
        else:
            builder.iadd(source[first:],target=dest[start:],control=operand[j])

def builder(family,segment_bits,guard,normalization):
    if family not in ("signed","unsigned"):raise ValueError(family)
    if normalization not in ("compact","barrel"):raise ValueError(normalization)
    module=hx if family=="unsigned" else cx
    old=(module._inplace_product if family=="unsigned"
           else module.signed_product_update)
    if family=="unsigned":module._inplace_product=controlled_row_update
    else:module.signed_product_update=signed_controlled_update
    old_n=module.normalize_barrel;old_u=module.denormalize_barrel
    if normalization=="compact":
        module.normalize_barrel=normalize_compact
        module.denormalize_barrel=denormalize_compact
    try:
        if family=="unsigned":
            if segment_bits!=6:raise ValueError("Unsigned only p6")
            cfg=hx.XAGHornerConfig(guard=guard,lookup="qrom",cleanup="mbuc")
            math=hx.screen(cfg)
            c=module.build(cfg)
        else:
            cfg=cx.CenteredConfig(segment_bits=segment_bits,guard=guard,
                 lookup="qrom",cleanup="mbuc")
            math=cx.screen(cfg)
            c=module.build(cfg)
    finally:
        if family=="unsigned":module._inplace_product=old
        else:module.signed_product_update=old
        module.normalize_barrel=old_n
        module.denormalize_barrel=old_u
    return c,math

def main():
    records=[]
    for family,p,guards in [("unsigned",6,(4,5)),("signed",6,(3,4)),("signed",7,(3,4))]:
        for guard in guards:
            for norm in ("compact","barrel"):
                try:
                    c,math=builder(family,p,guard,norm)
                    r=resource_counts(c)
                    out=dict(family=family,segment_bits=p,guard=guard,
                      normalization=norm,math=math,resources=r,
                      could_beat_incumbent=r["lower_bound_qc"]<CURRENT,
                      could_reach_20x=r["lower_bound_qc"]<TARGET)
                except Exception as e:
                    out=dict(family=family,segment_bits=p,guard=guard,
                       normalization=norm,error=f"{type(e).__name__}: {str(e)[:450]}")
                records.append(out)
                print(json.dumps(out),flush=True)
    path=Path("runs/horner_controlled_low_space.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(records,indent=2)+"\n")
if __name__=="__main__":main()
