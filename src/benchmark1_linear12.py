"""Independent p12 piecewise-linear inverse-sqrt phase oracle.

A *completely different approximation* from v5: 12 high fractional mantissa
bits index a linear-interpolation table; remaining eight low bits receive
one truncated variable-coefficient product. No square or quadratic stage.

Across ALL 3*2^20 normalized inputs, exact omitted-column math with
cut=5 has worst error below 0.68 x the allowed tolerance. The native
lookup implementation is screened for large table costs before considering
full quantum verification. No competitor source/table/circuit is imported.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts
from .benchmark1_controlled_cut import mul_controlled
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block

P=12;LO=20-P;Q=26;DROP=5;H=1<<LO;N=1<<P
TOL=2**-24
REFERENCE=2076880.
INCUMBENT=272300.78036771005

def make_table():
    j=np.arange(N,4*N,dtype=np.float64)
    A=np.zeros(4*N,dtype=np.int64)
    B=np.zeros(4*N,dtype=np.int64)
    A[N:]=np.rint((1<<Q)/np.sqrt(j/N)).astype(np.int64)
    B[N:]=np.rint((1<<Q)*(1/np.sqrt(j/N)-1/np.sqrt((j+1)/N))).astype(np.int64)
    return A,B

def screen():
    A,B=make_table();j=np.arange(N,4*N,dtype=np.float64)
    a=A[N:];b=B[N:]
    worst=0.;worst_index=0;failures=0
    for r in range(H):
        term=truncated_product(b,r,DROP)
        z=(j+r/H)/N
        phase=a/(1<<Q)-term/(1<<(Q+LO-DROP))
        error=np.abs(phase*np.sqrt(z)-1)/TOL
        failures+=int(np.count_nonzero(error>=1))
        ix=int(np.argmax(error))
        if float(error[ix])>worst:
            worst=float(error[ix]);worst_index=(int(j[ix])*H)+r
    return dict(passes=failures==0,failures=failures,
        normalized_inputs=3*N*H,worst_error_ratio=worst,
        worst_normalized_word=worst_index,
        segment_bits=P,low_bits=LO,precision=Q,partial_product_cut=DROP,
        coefficient_widths=[max(1,int(A.max()).bit_length()),
                            max(1,int(B.max()).bit_length())])

def build(*,native_lookup_only=False):
    import kickmix as km
    report=screen()
    if not report["passes"]:raise ValueError(str(report))
    A,B=make_table()
    wa,wb=report["coefficient_widths"]
    table=[int(av)|(int(bv)<<wa) for av,bv in zip(A,B)]
    ta=[int(av) for av in A]
    tb=[int(bv) for bv in B]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(b,x)
    c=b.alloc_qubits(wa+wb)
    with b.mark("linear_coefficients_lookup"):
        b.init_lookup(table=table,address=x[LO:],target=c)
    a,sl=c[:wa],c[wa:]
    if native_lookup_only:
        # Deliberately incomplete and unverified: measuring the size of the
        # initial 14-address-bit QROM *alone* gives a rigorous native
        # lower bound on every full oracle that uses this architecture.
        return b.finish_circuit()
    with b.mark("intercept_phase"):
        phase_with_pad(b,a,shift,Q,+1,10)
        b.del_lookup(table=ta,address=x[LO:],target=a)
        b.free(a)
    with b.mark("one_residual_multiplier"):
        product=b.alloc_qubits(wb+LO-DROP);b.reset(product)
        mul_controlled(b,sl,x[:LO],product,DROP)
        phase_with_pad(b,product,shift,Q+LO-DROP,-1,3)
        mul_controlled(b,sl,x[:LO],product,DROP,inverse=True)
        b.free(product)
    with b.mark("cleanup"):
        b.del_lookup(table=tb,address=x[LO:],target=sl)
        b.free(sl)
        denormalize_binary_block(b,x,flags,shift)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--math-only",action="store_true")
    p.add_argument("--lookup-only",action="store_true")
    p.add_argument("--output",type=Path,default=Path("out/linear_p12.kmb"))
    p.add_argument("--sample",type=int,default=0)
    a=p.parse_args()
    result=dict(math=screen())
    if not a.math_only:
        c=build(native_lookup_only=a.lookup_only)
        stats=resource_counts(c)
        result["resources"]=stats
        result["toffoli_floor_can_reach_20x"]=stats["lower_bound_qc"]<REFERENCE/20
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        result["sha256"]=hashlib.sha256(a.output.read_bytes()).hexdigest()
        if a.sample and not a.lookup_only:
            from .hybrid_verify import verify_samples,boundary_inputs
            result["coherent"]=verify_samples(c,boundary_inputs(random_count=a.sample))
            if result["coherent"]["failures"]:raise ValueError(str(result))
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
