"""Original low-peak-workspace coefficient-streamed phase oracle.

The verified v8 keeps A(27), B(11), C(12) coefficients in one 50-bit
bank. But their use is sequential: intercept, linear correction, then
quadratic correction. Here the coefficient bank is split into three
quantum lifetimes; each coefficient is loaded only for its stage and
measured/uncomputed before the next one starts.

Search native QROM versus original don’t-care GF(2) ANF loaders and
two-stage A+BC versus three-stage A+B+C pipelines. All variants
retain EXACT same 384 reachable coefficients and exact symmetric-square
arithmetic and binary-normalized exponent as verified v8. No third-party
competition source is used.

The changed circuit must pass the independent quantum verifier and actual
Q*C before promotion; Toffoli floors alone are not evidence of improvement.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from .benchmark1_pruned_anf_rom import binary_pruned_rom
from .benchmark1_binary_block_normalizer import (
    normalize_binary_block,denormalize_binary_block)
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled
from .benchmark1_symmetric_square_rows_v7 import square_symmetric_rows,exact_math_check
from .hybrid_kickmix import resource_counts
from . import benchmark1_cost_cooptim as cc

CUTS=(10,13,12)
INCUMBENT=241424.49222389946
REFERENCE=2076880.

def build(mode="split",first="anf",other="qrom"):
    if mode not in ("two","split"):raise ValueError(mode)
    if first not in ("qrom","anf") or other not in ("qrom","anf"):
        raise ValueError("Unsupported table synthesis mode")
    if not cc.full_screen(CUTS)["passes"]:raise ValueError("Polynomial screen failed")
    if not exact_math_check()["passes"]:raise ValueError("Square screen failed")
    import kickmix as km
    mat,_=optimized_coefficients("affine")
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    ta=[int(v) for v in A]
    tb=[int(v) for v in B]
    tc=[int(v) for v in C]
    tbc=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(builder,x)
    address=x[cc.LO:]
    def load(table,width,loader):
        reg=builder.alloc_qubits(width)
        if loader=="qrom":
            builder.init_lookup(table=table,address=address,target=reg)
        else:
            # XAG uses a multioutput GF(2) Möbius representation with
            # shared conditional ANDs, reversible at every used address.
            builder.reset(reg)
            binary_pruned_rom(builder,address,reg,table)
        return reg
    def erase(table,reg):
        builder.del_lookup(table=table,address=address,target=reg)
        builder.free(reg)
    with builder.mark("intercept_coefficient_lifetime"):
        reg=load(ta,wa,first)
        phase_with_pad(builder,reg,shift,cc.Q,+1,9)
        erase(ta,reg)
    dl,ds,dq=CUTS
    def linear(breg):
        prod=builder.alloc_qubits(wb+cc.LO-(dl-cc.BZ))
        builder.reset(prod)
        mul_controlled(builder,breg,x[:cc.LO],prod,dl-cc.BZ)
        phase_with_pad(builder,prod,shift,cc.Q+cc.LO-dl,-1,5)
        mul_controlled(builder,breg,x[:cc.LO],prod,dl-cc.BZ,inverse=True)
        builder.free(prod)
    def quadratic(creg):
        sq=builder.alloc_qubits(2*cc.LO-ds)
        builder.reset(sq)
        square_symmetric_rows(builder,x[:cc.LO],sq,ds)
        prod=builder.alloc_qubits(wc+len(sq)-(dq-cc.CZ))
        builder.reset(prod)
        mul_controlled(builder,sq,creg,prod,dq-cc.CZ)
        phase_with_pad(builder,prod,shift,cc.Q+2*cc.LO-ds-dq,+1,5)
        mul_controlled(builder,sq,creg,prod,dq-cc.CZ,inverse=True)
        builder.free(prod)
        square_symmetric_rows(builder,x[:cc.LO],sq,ds,inverse=True)
        builder.free(sq)
    if mode=="two":
        with builder.mark("shared_beta_gamma_lifetime"):
            bc=load(tbc,wb+wc,other)
            linear(bc[:wb])
            quadratic(bc[wb:])
            erase(tbc,bc)
    else:
        with builder.mark("linear_coefficient_lifetime"):
            b=load(tb,wb,other)
            linear(b)
            erase(tb,b)
        with builder.mark("quadratic_coefficient_lifetime"):
            c=load(tc,wc,other)
            quadratic(c)
            erase(tc,c)
    denormalize_binary_block(builder,x,flags,shift)
    return builder.finish_circuit()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--variants",nargs="+",default=
        ["two-anf-qrom","two-anf-anf","split-anf-qrom","split-qrom-qrom","split-anf-anf"])
    ap.add_argument("--samples",type=int,default=128)
    args=ap.parse_args()
    rows=[]
    for v in args.variants:
        mode,first,other=v.split("-")
        try:
            c=build(mode,first,other)
            stats=resource_counts(c)
            result=dict(variant=v,status="NATIVE_BUILT",native=stats,
                reaches_20x_by_toffoli_floor=stats["lower_bound_qc"]<REFERENCE/20,
                beats_v8_by_toffoli_floor=stats["lower_bound_qc"]<INCUMBENT)
            if stats["lower_bound_qc"]<INCUMBENT and args.samples:
                from .hybrid_verify import boundary_inputs,verify_samples
                result["sample"]=verify_samples(c,boundary_inputs(random_count=args.samples))
                if result["sample"]["failures"]:
                    result["status"]="COHERENT_SAMPLE_FAILED"
        except Exception as e:
            result=dict(variant=v,status="FAILED",
                        error=type(e).__name__+":"+str(e)[:450])
        rows.append(result)
        print(json.dumps(result),flush=True)
    p=Path("runs/coeff_stream_native_screen.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(rows,indent=2)+"\n")
if __name__=="__main__":main()
