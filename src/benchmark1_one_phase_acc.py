"""Original single fixed-point phase accumulator for quadratic oracle.

Instead of 3 phase-register rotations (76 non-Clifford Z_POW), form
a single signed fixed-point phase integer in the *existing intercept
coefficient* word. With the verified (dl,ds,dq)=(10,13,12):
  frac(A)=26, frac(L)=Q+LO-dl=29,
  frac(P)=Q+2LO-ds-dq=27.
Therefore form a coherent 30-bit accumulation:
  acc = 8*A - L + 4*P
where L and P are the exact independently verified omitted-column
partial products. One phase call uses frac_bits=29 and the same
normalization exponent shift s.

After phasing, invert P and L additions to recover exactly original A,
then use unchanged measured coefficient unlookups. Preserves v9
measurement-assisted cleanup and carry-free first rows.

This is a different circuit-construction algorithm from competitor
submissions. Full verification is required for any score.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import optimized_coefficients
from .benchmark1_pruned_anf_rom import binary_pruned_rom
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_first_row_specialization import mul_first_row,square_first_row
from .hybrid_kickmix import resource_counts
from . import benchmark1_cost_cooptim as cc

CUTS=(10,13,12)
REFERENCE=2076880.
V9=237684.49222389946
F=29

def build():
    assert cc.full_screen(CUTS)["passes"]
    import kickmix as km
    mat,_=optimized_coefficients("affine")
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    if wa+3!=30:raise RuntimeError("Unexpected 30bit phase accumulator")
    tA=[int(v) for v in A]
    tBC=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    tALL=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
          for a,b,c in zip(A,B,C)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(b,x)
    addr=x[cc.LO:]
    bank=b.alloc_qubits(wa+wb+wc)
    b.reset(bank)
    binary_pruned_rom(b,addr,bank,tALL)
    a=bank[:wa];beta=bank[wa:wa+wb];gamma=bank[wa+wb:]
    low=b.alloc_qubits(3);b.reset(low)
    acc=km.array(list(low)+list(a))
    if len(acc)!=30:raise RuntimeError("bad phase bank")
    dl,ds,dq=CUTS
    def linear(inverse=False):
        L=b.alloc_qubits(wb+cc.LO-(dl-cc.BZ))
        b.reset(L)
        mul_first_row(b,beta,x[:cc.LO],L,dl-cc.BZ)
        if inverse:b.iadd(L,target=acc)
        else:b.isub(L,target=acc)
        mul_first_row(b,beta,x[:cc.LO],L,dl-cc.BZ,inverse=True)
        b.free(L)
    def quadratic(inverse=False):
        sq=b.alloc_qubits(2*cc.LO-ds);b.reset(sq)
        square_first_row(b,x[:cc.LO],sq,ds)
        P=b.alloc_qubits(wc+len(sq)-(dq-cc.CZ))
        b.reset(P)
        mul_first_row(b,sq,gamma,P,dq-cc.CZ)
        if inverse:b.isub(P,target=acc[2:])
        else:b.iadd(P,target=acc[2:])
        mul_first_row(b,sq,gamma,P,dq-cc.CZ,inverse=True)
        b.free(P)
        square_first_row(b,x[:cc.LO],sq,ds,inverse=True)
        b.free(sq)
    with b.mark("single_phase_form_accumulator"):
        linear()
        quadratic()
    with b.mark("single_phase_30bit"):
        phase_with_pad(b,acc,shift,F,+1,10)
    with b.mark("single_phase_restore_A"):
        quadratic(inverse=True)
        linear(inverse=True)
    b.free(low)
    b.del_lookup(table=tA,address=addr,target=a)
    b.free(a)
    b.del_lookup(table=tBC,address=addr,target=bank[wa:])
    b.free(bank[wa:])
    denormalize_binary_block(b,x,flags,shift)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/one_phase_accumulator.kmb"))
    p.add_argument("--samples",type=int,default=256)
    args=p.parse_args()
    c=build()
    stats=resource_counts(c)
    out=dict(native=stats,verified_v9=V9,
        toffoli_floor_below_v9=stats["lower_bound_qc"]<V9,
        toffoli_floor_can_reach_20x=stats["lower_bound_qc"]<REFERENCE/20)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    out["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if args.samples and out["toffoli_floor_below_v9"]:
        from .hybrid_verify import boundary_inputs,verify_samples
        out["coherent"]=verify_samples(c,boundary_inputs(random_count=args.samples))
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
