"""Original early coefficient erasure in a single lookup with shared ANF input.

The verified v8 loads A(27),B(11),C(12) together, then erases A after
its phase, but holds B through the quadratic stage even though its
linear phase is complete. Quantum-controlled B unlookup immediately
after the linear phase frees B before the large square/C-product
workspace, reducing peak width by up to eleven qubits.

After the quadratic phase, C is likewise erased separately. Explore
measurement-assisted native QROM with padded clean targets, versus the
coherent self-inverse shared ANF Boolean decoder. The quantum table
is unchanged on every legal input; no competitors' circuits are read.
"""
from __future__ import annotations
from pathlib import Path
import argparse,json
from .benchmark1_dontcare_boolean_rom import optimized_coefficients
from .benchmark1_pruned_anf_rom import binary_pruned_rom
from .benchmark1_binary_block_normalizer import (
    normalize_binary_block,denormalize_binary_block)
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled
from .benchmark1_symmetric_square_rows_v7 import square_symmetric_rows
from .hybrid_kickmix import resource_counts
from . import benchmark1_cost_cooptim as cc

CUTS=(10,13,12);Q=26;REFERENCE=2076880.
V8=241424.49222389946

def build(erase_b="qrom",erase_c="qrom"):
    if erase_b not in ("qrom","anf") or erase_c not in ("qrom","anf"):
        raise ValueError("Invalid decoder")
    import kickmix as km
    mat,_=optimized_coefficients("affine")
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    tA=list(map(int,A));tB=list(map(int,B));tC=list(map(int,C))
    packed=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
            for a,b,c in zip(A,B,C)]
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(22,name="x")
    flags,scale=normalize_binary_block(builder,x)
    address=x[cc.LO:]
    coeff=builder.alloc_qubits(wa+wb+wc)
    with builder.mark("coefficient_decode"):
        builder.reset(coeff)
        binary_pruned_rom(builder,address,coeff,packed)
    a=coeff[:wa];beta=coeff[wa:wa+wb];gamma=coeff[wa+wb:]
    with builder.mark("phase_intercept"):
        phase_with_pad(builder,a,scale,Q,+1,9)
        builder.del_lookup(table=tA,address=address,target=a)
        builder.free(a)
    dl,ds,dq=CUTS
    with builder.mark("linear_phase"):
        prod=builder.alloc_qubits(wb+cc.LO-(dl-cc.BZ));builder.reset(prod)
        mul_controlled(builder,beta,x[:cc.LO],prod,dl-cc.BZ)
        phase_with_pad(builder,prod,scale,Q+cc.LO-dl,-1,5)
        mul_controlled(builder,beta,x[:cc.LO],prod,dl-cc.BZ,inverse=True)
        builder.free(prod)
    def erase(table,reg,width,mode):
        if mode=="anf":
            binary_pruned_rom(builder,address,reg,table)
        else:
            # The native measured QROM unlookup requires enough target
            # qubits for its unary address decomposition. Append clean
            # zero-valued output bits without changing table semantics.
            out_width=max(width,16)
            extra=builder.alloc_qubits(out_width-width)
            builder.reset(extra)
            union=km.array(list(reg)+list(extra))
            builder.del_lookup(table=table,address=address,target=union)
            builder.free(extra)
        builder.free(reg)
    with builder.mark("early_beta_erasure"):
        erase(tB,beta,wb,erase_b)
    with builder.mark("quadratic_phase"):
        square=builder.alloc_qubits(2*cc.LO-ds);builder.reset(square)
        square_symmetric_rows(builder,x[:cc.LO],square,ds)
        quadratic=builder.alloc_qubits(wc+len(square)-(dq-cc.CZ))
        builder.reset(quadratic)
        mul_controlled(builder,square,gamma,quadratic,dq-cc.CZ)
        phase_with_pad(builder,quadratic,scale,Q+2*cc.LO-ds-dq,+1,5)
        mul_controlled(builder,square,gamma,quadratic,dq-cc.CZ,inverse=True)
        builder.free(quadratic)
        square_symmetric_rows(builder,x[:cc.LO],square,ds,inverse=True)
        builder.free(square)
    with builder.mark("gamma_erasure"):
        erase(tC,gamma,wc,erase_c)
    denormalize_binary_block(builder,x,flags,scale)
    return builder.finish_circuit()

def main():
    r=[]
    for mb,mc in (("qrom","qrom"),("anf","qrom"),("qrom","anf"),("anf","anf")):
        try:
            c=build(mb,mc)
            s=resource_counts(c)
            record=dict(erase_b=mb,erase_c=mc,native=s,
                beats_v8_by_toffoli_floor=s["lower_bound_qc"]<V8,
                reaches_20x_by_floor=s["lower_bound_qc"]<REFERENCE/20)
            if record["beats_v8_by_toffoli_floor"]:
                from .hybrid_verify import boundary_inputs,verify_samples
                record["coherent"]=verify_samples(c,boundary_inputs(random_count=256))
        except Exception as e:
            record=dict(erase_b=mb,erase_c=mc,error=type(e).__name__+":"+str(e)[:600])
        print(json.dumps(record),flush=True)
        r.append(record)
    path=Path('runs/early_coefficient_erasure.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2)+'\n')
if __name__=="__main__":main()
