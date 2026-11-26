"""Hardware-aware sparse-coefficient Benchmark-I candidate (EXPERIMENTAL).

No claim of a score until native kickmix and full independent verifier PASS.
Fits interval coefficients with B divisible by 4 and C divisible by 16;
uses exact shifted partial products and omits their known-zero input bits.
"""
from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import argparse
import json
import numpy as np

Q=26
P=7
DROPS=(10,12,12)
B_ZERO=2
C_ZERO=4
MAX_X=3_133_452
TOL=2.**-24
VERIFIED_QC=430259.61862786254

@lru_cache(None)
def coefficients():
    from .priority1_model import Approximation, coefficient_matrix
    seed=coefficient_matrix(Approximation('chebyshev',2,P,Q))
    u=np.unique(np.r_[np.linspace(0.,1.,65), .5+.5*np.cos(np.pi*np.arange(65)/64)])
    table=np.array(seed,copy=True)
    for j in range(1<<P,4<<P):
        z=(j+u)/(1<<P)
        desired=(1<<Q)/np.sqrt(z)
        weight=np.sqrt(z)/((1<<Q)*TOL)
        a0,b0,c0=(int(v) for v in seed[j])
        best=(float('inf'),None)
        for dc in range(-3,4):
            c=(round(c0/16)+dc)*16
            if c<0: continue
            value=desired-c*u*u
            secant=value[0]-value[-1]
            for db in range(-4,5):
                b=(round(secant/4)+db)*4
                if b<0: continue
                wanted=desired+b*u-c*u*u
                amid=(float(wanted.min())+float(wanted.max()))/2
                for a in (int(np.floor(amid)),int(np.ceil(amid))):
                    if a<0:continue
                    worst=float(np.max(np.abs((a-b*u+c*u*u-desired)*weight)))
                    if worst<best[0]:
                        best=(worst,(a,b,c))
        if best[1] is None:raise RuntimeError('coefficient fit failed')
        table[j]=best[1]
    assert np.all((table[:,1]&3)==0) and np.all((table[:,2]&15)==0)
    table.setflags(write=False)
    return table

def math_screen(chunk=150000):
    from .priority1_exact import truncated_product
    mat=coefficients()
    worst=0.;worst_x=1;failures=0
    dl,ds,dq=DROPS
    for start in range(1,MAX_X+1,chunk):
        x=np.arange(start,min(start+chunk,MAX_X+1),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        norm=x<<(20-2*k)
        A,B,C=mat[norm>>13].T
        r=norm&8191
        linear=truncated_product(B,r,dl)
        sq=truncated_product(r,r,ds)
        quad=truncated_product(C,sq,dq)
        phase=np.ldexp(A/(1<<Q)-linear/2**(Q+13-dl)+
                       quad/2**(Q+26-ds-dq),-k)
        error=np.abs(phase*np.sqrt(x)-1)/TOL
        i=int(np.argmax(error))
        if float(error[i])>worst:worst=float(error[i]);worst_x=int(x[i])
        failures+=int(np.count_nonzero(error>1))
    return dict(passes=not failures,failures=failures,
                worst_error_over_tolerance=worst,worst_x=worst_x,
                checked_inputs=MAX_X)

def build():
    import kickmix as km
    from .priority1_truncated import _mul_cut
    from .benchmark1_compact_normalize import normalize_compact, denormalize_compact
    from .hybrid_kickmix import _phase
    numerical=math_screen()
    if not numerical['passes'] or numerical['worst_error_over_tolerance']>=.95:
        raise ValueError('Exhaustive phase gate failed: '+str(numerical))
    mat=coefficients()
    A=mat[:,0]
    B=mat[:,1]>>B_ZERO
    C=mat[:,2]>>C_ZERO
    wa=int(A.max()).bit_length()
    wb=int(B.max()).bit_length()
    wc=int(C.max()).bit_length()
    ta=[int(a) for a in A]
    tbc=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(22,name='x')
    unused,k=normalize_compact(builder,x)
    with builder.mark('constant_A'):
        a=builder.alloc_qubits(wa)
        builder.init_lookup(table=ta,address=x[13:],target=a)
        _phase(builder,a,k,Q,+1)
        builder.del_lookup(table=ta,address=x[13:],target=a)
        builder.free(a)
    with builder.mark('compressed_BC_lookup'):
        bc=builder.alloc_qubits(wb+wc)
        builder.init_lookup(table=tbc,address=x[13:],target=bc)
    beta,gamma=bc[:wb],bc[wb:]
    dl,ds,dq=DROPS
    with builder.mark('linear_phase'):
        d=dl-B_ZERO
        result=builder.alloc_qubits(len(beta)+13-d)
        builder.reset(result)
        _mul_cut(builder,beta,x[:13],result,d)
        _phase(builder,result,k,Q+13-dl,-1)
        _mul_cut(builder,beta,x[:13],result,d,inverse=True)
        builder.free(result)
    with builder.mark('quadratic_phase'):
        square=builder.alloc_qubits(26-ds)
        builder.reset(square)
        _mul_cut(builder,x[:13],x[:13],square,ds)
        d=dq-C_ZERO
        result=builder.alloc_qubits(len(gamma)+len(square)-d)
        builder.reset(result)
        _mul_cut(builder,gamma,square,result,d)
        _phase(builder,result,k,Q+26-ds-dq,+1)
        _mul_cut(builder,gamma,square,result,d,inverse=True)
        builder.free(result)
        _mul_cut(builder,x[:13],x[:13],square,ds,inverse=True)
        builder.free(square)
    with builder.mark('cleanup'):
        builder.del_lookup(table=tbc,address=x[13:],target=bc)
        builder.free(bc)
    denormalize_compact(builder,x,unused,k)
    return builder.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('command',choices=('screen','build','verify'))
    p.add_argument('--output',type=Path)
    p.add_argument('--circuit',type=Path)
    p.add_argument('--random-count',type=int,default=4096)
    a=p.parse_args()
    if a.command=='screen':
        report=dict(math=math_screen())
    elif a.command=='build':
        if a.output is None:p.error('--output required')
        from .hybrid_kickmix import resource_counts
        circuit=build()
        resources=resource_counts(circuit)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        circuit.to_file(a.output,format='kmb')
        report=dict(math=math_screen(),resources=resources,file=str(a.output),
          toffoli_floor_beats_verified=resources['lower_bound_qc']<VERIFIED_QC,
          toffoli_floor_can_reach_10x=resources['lower_bound_qc']<207688.)
    else:
        if a.circuit is None:p.error('--circuit required')
        import kickmix as km
        from .hybrid_verify import verify_samples,boundary_inputs
        circuit=km.Circuit.from_file(a.circuit,format='kmb')
        report=verify_samples(circuit,boundary_inputs(random_count=a.random_count))
    print(json.dumps(report,indent=2,sort_keys=True))
    if report.get('failures',0):raise SystemExit(1)

if __name__=='__main__':main()
