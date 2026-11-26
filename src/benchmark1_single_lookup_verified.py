"""Load A/B/C once; early measurement-corrected cleanup of phase-consumed A.

Unlike XOR-buffer-morph, init_lookup only acts on a clean fresh register.
"""
from __future__ import annotations
from pathlib import Path
import argparse,json
from .priority1_model import Approximation,coefficient_matrix
from .priority1_exact import full_cut_screen
from .priority1_truncated import _mul_cut
from .benchmark1_compact_normalize import normalize_compact,denormalize_compact
from .hybrid_kickmix import _phase,resource_counts
CONFIG=Approximation("chebyshev",2,7,26)
DROPS=(10,13,10)
REFERENCE_QC=430259.61862786254
VERIFIED_SHA256="1726e2a438e4296590210c2422004ed1feda9153c09f28a9edc31dce930e0be6"
VERIFIED_COST=367309.6152203282

def build(*,drops=DROPS):
    import kickmix as km
    check=full_cut_screen(CONFIG,tuple(drops))
    if not check["passes"] or check["worst_error_over_tolerance"]>=.85:
        raise ValueError("Full-domain mathematical approximation invalid")
    mat=coefficient_matrix(CONFIG)
    wa,wb,wc=[max(1,int(mat[:,i].max()).bit_length()) for i in range(3)]
    table_a=[int(row[0]) for row in mat]
    table_bc=[int(row[1])|(int(row[2])<<wb) for row in mat]
    table_abc=[a|(bc<<wa) for a,bc in zip(table_a,table_bc)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_compact(b,x)
    with b.mark("joint_coefficient_load"):
        word=b.alloc_qubits(wa+wb+wc)
        b.init_lookup(table=table_abc,address=x[13:],target=word)
    a=word[:wa];bc=word[wa:]
    beta=bc[:wb];gamma=bc[wb:]
    with b.mark("phase_A_and_early_A_erasure"):
        _phase(b,a,shift,26,+1)
        b.del_lookup(table=table_a,address=x[13:],target=a)
        b.free(a)
    dl,ds,dq=drops
    mul=_mul_cut
    with b.mark("linear_phase"):
        linear=b.alloc_qubits(len(beta)+13-dl)
        b.reset(linear)
        mul(b,beta,x[:13],linear,dl)
        _phase(b,linear,shift,26+13-dl,-1)
        mul(b,beta,x[:13],linear,dl,inverse=True)
        b.free(linear)
    with b.mark("quadratic_phase"):
        sq=b.alloc_qubits(26-ds);b.reset(sq)
        _mul_cut(b,x[:13],x[:13],sq,ds)
        quad=b.alloc_qubits(len(gamma)+len(sq)-dq);b.reset(quad)
        mul(b,gamma,sq,quad,dq)
        _phase(b,quad,shift,26+26-ds-dq,+1)
        mul(b,gamma,sq,quad,dq,inverse=True);b.free(quad)
        _mul_cut(b,x[:13],x[:13],sq,ds,inverse=True);b.free(sq)
    with b.mark("BC_erase_and_restore"):
        b.del_lookup(table=table_bc,address=x[13:],target=bc)
        b.free(bc)
    denormalize_compact(b,x,flags,shift)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/joint_coefficients.kmb"))
    p.add_argument("--drops",default="10,13,10",help="linear,square,quadratic partial-product omissions")
    a=p.parse_args()
    cuts=tuple(int(v) for v in a.drops.split(","))
    if len(cuts)!=3: raise ValueError("Exactly 3 cuts required")
    circuit=build(drops=cuts)
    stats=resource_counts(circuit)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(a.output,format="kmb")
    import hashlib
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    if cuts==DROPS and digest!=VERIFIED_SHA256:
        raise RuntimeError("Pinned verified circuit hash mismatch: "+digest)
    print(json.dumps(dict(**stats,cuts=cuts,reference_qc=REFERENCE_QC,verified_qc=VERIFIED_COST,circuit_sha256=digest,file=str(a.output)),indent=2))

if __name__=="__main__":main()
