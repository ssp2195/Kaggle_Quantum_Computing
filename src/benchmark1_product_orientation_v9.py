"""Independent commutative product orientation optimization on certified v9.

For partial-column truncated products T_d(a,b), the bit-product
multiset is commutative: T_d(a,b)=T_d(b,a), even when lower columns are
omitted. Native controlled adders, however, are asymmetrical: cost
depends on which operand supplies row controls. Test both orientations
for linear B(11)*r(13) and quadratic sq(13)*gamma(12), using the
certified carry-free first actual row and its exact reverse erasure.

Every alternative must pass native quantum coherence before being
considered for the independent full-domain verifier.
"""
from __future__ import annotations
import json
from pathlib import Path
from . import benchmark1_first_row_specialization as f
from .hybrid_kickmix import resource_counts

REF=2076880.;V8=241424.49222389946

def build(swap_linear=False,swap_quad=False):
    old=f.mul_first_row
    from .benchmark1_cost_cooptim import coefficients,BZ,CZ
    mat=coefficients((10,13,12))
    bw=max(1,int((mat[:,1]>>BZ).max()).bit_length())
    cw=max(1,int((mat[:,2]>>CZ).max()).bit_length())
    if bw==cw:raise RuntimeError("Ambiguous product shape")
    def orient(b,a,other,dest,drop,inverse=False):
        if (len(a),len(other))==(bw,13) and swap_linear:
            return old(b,other,a,dest,drop,inverse)
        if (len(a),len(other))==(13,cw) and swap_quad:
            return old(b,other,a,dest,drop,inverse)
        return old(b,a,other,dest,drop,inverse)
    f.mul_first_row=orient
    try:c=f.build(linear=True,square=True)
    finally:f.mul_first_row=old
    return c

def main():
    rows=[]
    for l,q in ((False,False),(True,False),(False,True),(True,True)):
        try:
            c=build(l,q)
            cost=resource_counts(c)
            out=dict(swap_linear=l,swap_quadratic=q,native=cost,
                could_improve_previous_at_floor=cost["lower_bound_qc"]<V8,
                can_reach_20x_at_floor=cost["lower_bound_qc"]<REF/20)
            if out["could_improve_previous_at_floor"]:
                from .hybrid_verify import verify_samples,boundary_inputs
                out["coherent"]=verify_samples(c,boundary_inputs(random_count=256))
        except Exception as e:
            out=dict(swap_linear=l,swap_quadratic=q,error=type(e).__name__+":"+str(e)[:420])
        rows.append(out);print(json.dumps(out),flush=True)
    p=Path('runs/firstrow-product-orientation.json')
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=="__main__":main()
