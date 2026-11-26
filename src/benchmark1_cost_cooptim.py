"""Cost-aware coefficient refit under EXACT native partial-product cuts.

Candidate family is the already verified p7 q26 piecewise quadratic
phase oracle, but integer coefficient synthesis is changed from fitting
a smooth polynomial to optimizing the actual dropped-column quantum
arithmetic. Fitting samples a deterministic low-bit-stratified grid,
then repairs individual failing intervals against ALL 8192 mantissa
residuals. The resulting table is certified on every normalized
interval/residual pair before constructing any native circuit.

Experimental; no score until native full pinned verification.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import argparse,json
import numpy as np
from .priority1_model import Approximation,coefficient_matrix
from .priority1_exact import truncated_product
from .hybrid_model import MAX_X,TOL
from .benchmark1_compact_normalize import normalize_compact,denormalize_compact
from .hybrid_kickmix import _phase,resource_counts
from .benchmark1_controlled_cut import mul_controlled
from .priority1_truncated import _mul_cut

P=7
LO=13
Q=26
BZ=2
CZ=4
CURRENT=302641.04424652964
H=1<<LO

def _mesh():
    return np.unique(np.r_[np.linspace(0,H-1,121,dtype=int),
                  np.arange(0,H,257),
                  np.random.default_rng(10).integers(0,H,120)]).astype(np.int64)

def _error_for_row(j,r,A,B,C,cuts,square=None):
    dl,ds,dq=cuts
    z=(j+r/H)/(1<<P)
    lin=truncated_product(B>>BZ,r,dl-BZ)
    if square is None:square=truncated_product(r,r,ds)
    quad=truncated_product(C>>CZ,square,dq-CZ)
    phase=A/(1<<Q)-lin/2**(Q+LO-dl)+quad/2**(Q+2*LO-ds-dq)
    return np.abs(phase*np.sqrt(z)-1)/TOL

def _local_refit(j,center,cuts,rs,*,span_b,span_c):
    dl,ds,dq=cuts
    desired=(1<<Q)/np.sqrt((j+rs/H)/(1<<P))
    weight=np.sqrt((j+rs/H)/(1<<P))/((1<<Q)*TOL)
    square=truncated_product(rs,rs,ds)
    _,b0,c0=map(int,center)
    best=(float('inf'),None)
    for dc in range(-span_c,span_c+1):
        c=c0+dc*(1<<CZ)
        if c<0:continue
        quad=truncated_product(c>>CZ,square,dq-CZ)/2**(2*LO-ds-dq)
        for db in range(-span_b,span_b+1):
            b=b0+db*(1<<BZ)
            if b<0:continue
            lin=truncated_product(b>>BZ,rs,dl-BZ)/2**(LO-dl)
            wanted=desired+lin-quad
            mid=(float(wanted.max())+float(wanted.min()))/2
            for a in (int(np.floor(mid)),int(np.ceil(mid))):
                if a<0:continue
                loss=float(np.max(np.abs((a-wanted)*weight)))
                if loss<best[0]:
                    best=(loss,(a,b,c))
    if best[1] is None:raise ValueError('No nonnegative coefficient fit')
    return best

@lru_cache(maxsize=4)
def coefficients(cuts=(10,12,13)):
    if cuts not in ((10,12,13),(10,13,12)):
        raise ValueError('Only screened cut settings supported')
    mesh=_mesh()
    seed=coefficient_matrix(Approximation("chebyshev",2,P,Q))
    out=np.array(seed,copy=True)
    dl,ds,dq=cuts
    square=truncated_product(mesh,mesh,ds)
    for j in range(1<<P,4<<P):
        z=(j+mesh/H)/(1<<P)
        desired=(1<<Q)/np.sqrt(z)
        weight=np.sqrt(z)/((1<<Q)*TOL)
        nodes=.5+.5*np.cos(np.pi*(np.arange(3)+.5)/3)
        v=np.linalg.solve(
            np.polynomial.polynomial.polyvander(nodes,2),
            1/np.sqrt((j+nodes)/(1<<P)))*(1<<Q)
        b0=round(-float(v[1])/(1<<BZ))
        c0=round(float(v[2])/(1<<CZ))
        best=(float('inf'),None)
        for dc in range(-6,7):
            c=(c0+dc)*(1<<CZ)
            if c<0:continue
            quad=truncated_product(c>>CZ,square,dq-CZ)/2**(2*LO-ds-dq)
            for db in range(-6,7):
                b=(b0+db)*(1<<BZ)
                if b<0:continue
                lin=truncated_product(b>>BZ,mesh,dl-BZ)/2**(LO-dl)
                wanted=desired+lin-quad
                mid=(float(wanted.max())+float(wanted.min()))/2
                for a in (int(np.floor(mid)),int(np.ceil(mid))):
                    if a<0:continue
                    loss=float(np.max(np.abs((a-wanted)*weight)))
                    if loss<best[0]:best=(loss,(a,b,c))
        if best[1] is None:raise ValueError('No fit at '+str(j))
        out[j]=best[1]
    # Every possible 13-bit residual is checked in every normalized
    # interval (strictly stronger than the competition's restricted x).
    dense=np.arange(H,dtype=np.int64)
    for j in range(1<<P,4<<P):
        A,B,C=map(int,out[j])
        err=_error_for_row(j,dense,A,B,C,cuts)
        if float(np.max(err))>1:
            _,best=_local_refit(j,out[j],cuts,dense,span_b=15,span_c=15)
            out[j]=best
    out.setflags(write=False)
    return out

def full_screen(cuts=(10,12,13)):
    c=coefficients(cuts)
    worst=0.;worst_j=0;worst_r=0;failures=0
    residual=np.arange(H,dtype=np.int64)
    for j in range(1<<P,4<<P):
        A,B,C=map(int,c[j])
        e=_error_for_row(j,residual,A,B,C,cuts)
        n=int(np.count_nonzero(e>1));failures+=n
        ix=int(np.argmax(e))
        if float(e[ix])>worst:worst=float(e[ix]);worst_j=j;worst_r=int(residual[ix])
    return dict(passes=failures==0,failures=failures,
        all_normalized_inputs=(3*(1<<P))*H,
        worst_error_over_tolerance=worst,worst_j=worst_j,
        worst_residual=worst_r,cuts=list(cuts))

def build(cuts=(10,12,13),*,swap_quadratic=False):
    import kickmix as km
    cert=full_screen(cuts)
    if not cert['passes'] or cert['worst_error_over_tolerance']>=1:
        raise ValueError('Exhaustive normalized phase certificate rejected: '+str(cert))
    mat=coefficients(cuts)
    A,B,C=mat[:,0],mat[:,1]>>BZ,mat[:,2]>>CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    tA=[int(a) for a in A]
    tBC=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    table=[int(a)|(int(bc)<<wa) for a,bc in zip(tA,tBC)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,k=normalize_compact(b,x)
    word=b.alloc_qubits(wa+wb+wc)
    b.init_lookup(table=table,address=x[LO:],target=word)
    a=word[:wa];bc=word[wa:];beta,gamma=bc[:wb],bc[wb:]
    with b.mark("constant_A_phase"):
        _phase(b,a,k,Q,+1)
        b.del_lookup(table=tA,address=x[LO:],target=a)
        b.free(a)
    dl,ds,dq=cuts
    with b.mark("linear_phase"):
        drop=dl-BZ
        product=b.alloc_qubits(wb+LO-drop)
        b.reset(product)
        mul_controlled(b,beta,x[:LO],product,drop)
        _phase(b,product,k,Q+LO-dl,-1)
        mul_controlled(b,beta,x[:LO],product,drop,inverse=True)
        b.free(product)
    with b.mark("quadratic_phase"):
        sq=b.alloc_qubits(2*LO-ds)
        b.reset(sq)
        _mul_cut(b,x[:LO],x[:LO],sq,ds)
        drop=dq-CZ
        product=b.alloc_qubits(wc+len(sq)-drop)
        b.reset(product)
        first,second=(sq,gamma) if swap_quadratic else (gamma,sq)
        mul_controlled(b,first,second,product,drop)
        _phase(b,product,k,Q+2*LO-ds-dq,+1)
        mul_controlled(b,first,second,product,drop,inverse=True)
        b.free(product)
        _mul_cut(b,x[:LO],x[:LO],sq,ds,inverse=True)
        b.free(sq)
    b.del_lookup(table=tBC,address=x[LO:],target=bc)
    b.free(bc)
    denormalize_compact(b,x,flags,k)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("command",choices=("screen","build","verify"))
    p.add_argument("--cuts",default="10,12,13")
    p.add_argument("--swap-quadratic",action="store_true")
    p.add_argument("--output",type=Path)
    p.add_argument("--circuit",type=Path)
    p.add_argument("--random-count",type=int,default=2048)
    a=p.parse_args()
    cuts=tuple(map(int,a.cuts.split(',')))
    if a.command=="screen":
        result={"math":full_screen(cuts)}
    elif a.command=="build":
        if a.output is None:p.error("--output required")
        c=build(cuts,swap_quadratic=a.swap_quadratic)
        stats=resource_counts(c)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        result=dict(math=full_screen(cuts),resources=stats,
            toffoli_floor_beats_incumbent=stats["lower_bound_qc"]<CURRENT,
            toffoli_floor_can_reach_ten_x=stats["lower_bound_qc"]<207688,
            output=str(a.output))
    else:
        if a.circuit is None:p.error("--circuit required")
        from .hybrid_verify import verify_samples,boundary_inputs
        import kickmix as km
        c=km.Circuit.from_file(a.circuit,format="kmb")
        result=verify_samples(c,boundary_inputs(random_count=a.random_count))
    print(json.dumps(result,indent=2))
    if result.get("failures",0) or result.get("math",{}).get("failures",0):
        raise SystemExit(1)

if __name__=="__main__":main()
