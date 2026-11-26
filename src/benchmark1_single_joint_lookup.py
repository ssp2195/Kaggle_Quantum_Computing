"""Single QROM A/B/C lookup with early measurement-based A-only erasure.

Native kickmix must verify that del_lookup of A (originally packed with BC)
preserves coherent phase and the still-live BC bits. Not a score claim.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .benchmark1_sparse_coeff import coefficients,math_screen
from .benchmark1_compact_normalize import normalize_compact,denormalize_compact
from .hybrid_kickmix import _phase,resource_counts

DROPS=(10,12,12)
VERIFIED_QC=430259.61862786254

def build(*,controlled=False):
    import kickmix as km
    from .priority1_truncated import _mul_cut
    if controlled:
        from .benchmark1_controlled_cut import mul_controlled as multiplier
    else:multiplier=_mul_cut
    cert=math_screen()
    if not cert['passes'] or cert['worst_error_over_tolerance']>=.95:
        raise ValueError(str(cert))
    data=coefficients()
    A=data[:,0]
    B=data[:,1]>>2
    C=data[:,2]>>4
    wa,wb,wc=[max(1,int(s.max()).bit_length()) for s in (A,B,C)]
    ta=[int(a) for a in A]
    tbc=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    table=[int(a)|(int(bc)<<wa) for a,bc in zip(ta,tbc)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name='x')
    unused,k=normalize_compact(b,x)
    storage=b.alloc_qubits(wa+wb+wc)
    with b.mark('single_ABC_lookup'):
        b.init_lookup(table=table,address=x[13:],target=storage)
    a=storage[:wa]
    bc=storage[wa:]
    beta,gamma=bc[:wb],bc[wb:]
    with b.mark('A_phase_and_early_measurement_erasure'):
        _phase(b,a,k,26,+1)
        b.del_lookup(table=ta,address=x[13:],target=a)
        b.free(a)
    dl,ds,dq=DROPS
    with b.mark('linear_phase'):
        d=dl-2
        product=b.alloc_qubits(len(beta)+13-d);b.reset(product)
        multiplier(b,beta,x[:13],product,d)
        _phase(b,product,k,26+13-dl,-1)
        multiplier(b,beta,x[:13],product,d,inverse=True)
        b.free(product)
    with b.mark('quadratic_phase'):
        square=b.alloc_qubits(26-ds);b.reset(square)
        multiplier(b,x[:13],x[:13],square,ds)
        d=dq-4
        product=b.alloc_qubits(len(gamma)+len(square)-d);b.reset(product)
        multiplier(b,gamma,square,product,d)
        _phase(b,product,k,26+26-ds-dq,+1)
        multiplier(b,gamma,square,product,d,inverse=True)
        b.free(product)
        multiplier(b,x[:13],x[:13],square,ds,inverse=True)
        b.free(square)
    with b.mark('BC_cleanup'):
        b.del_lookup(table=tbc,address=x[13:],target=bc)
        b.free(bc)
    denormalize_compact(b,x,unused,k)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('command',choices=('build','verify'))
    p.add_argument('--controlled',action='store_true')
    p.add_argument('--output',type=Path)
    p.add_argument('--circuit',type=Path)
    p.add_argument('--random-count',type=int,default=2048)
    a=p.parse_args()
    if a.command=='build':
        if a.output is None:p.error('--output required')
        c=build(controlled=a.controlled)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format='kmb')
        r=dict(resources=resource_counts(c),math=math_screen(),output=str(a.output))
        r['can_beat_reference_native_floor']=r['resources']['lower_bound_qc']<VERIFIED_QC
        r['can_reach_10x_native_floor']=r['resources']['lower_bound_qc']<207688
    else:
        if a.circuit is None:p.error('--circuit required')
        from .hybrid_verify import verify_samples,boundary_inputs
        import kickmix as km
        c=km.Circuit.from_file(a.circuit,format='kmb')
        r=verify_samples(c,boundary_inputs(random_count=a.random_count))
    print(json.dumps(r,indent=2))
    if r.get('failures',0):raise SystemExit(1)
if __name__=='__main__':main()
