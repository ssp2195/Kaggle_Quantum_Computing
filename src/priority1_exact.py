"""Full finite-domain screen for deliberately omitted partial products.

This exactly models the integer arithmetic of _mul_cut, not a floating
approximation to floor(a*b/2**drop). Every specified Benchmark I input is
evaluated; this does not substitute for kickmix's quantum verifier.
"""
import numpy as np
from .hybrid_model import MAX_X, TOL
from .priority1_model import Approximation, coefficient_matrix


def truncated_product(a, b, drop: int):
    a=np.asarray(a,dtype=np.int64)
    b=np.asarray(b,dtype=np.int64)
    omitted=np.zeros(np.broadcast_shapes(a.shape,b.shape),dtype=np.int64)
    for bit in range(drop):
        omitted += ((b>>bit)&1)*(a&((1<<(drop-bit))-1))*(1<<bit)
    return (a*b-omitted)>>drop


def full_cut_screen(config: Approximation, drops: tuple[int,int,int],
                    max_x: int = MAX_X):
    if config.degree!=2 or len(drops)!=3 or any(not isinstance(d,int) or d<0 for d in drops):
        raise ValueError("Only quadratic nonnegative integer drop triples supported")
    linear_drop,square_drop,quadratic_drop=drops
    lo,q=config.low_bits,config.precision
    coeffs=coefficient_matrix(config)
    if (linear_drop >= lo or square_drop >= lo*2 or
            quadratic_drop >= lo*2-square_drop):
        raise ValueError("Drop count exceeds supported phase-register ranges")
    worst,failures,worst_x=0.,0,1
    for start in range(1,max_x+1,200_000):
        x=np.arange(start,min(start+200_000,max_x+1),dtype=np.int64)
        exponent=np.floor(np.log2(x)).astype(np.int64)//2
        normalized=x<<(20-2*exponent)
        idx=normalized>>lo
        mantissa=normalized&((1<<lo)-1)
        a,b,c=(coeffs[idx,i] for i in range(3))
        lin=truncated_product(b,mantissa,linear_drop)
        square=truncated_product(mantissa,mantissa,square_drop)
        quad=truncated_product(c,square,quadratic_drop)
        phase=np.ldexp(
            a/2**q-lin/2**(q+lo-linear_drop)+
            quad/2**(q+2*lo-square_drop-quadratic_drop),
            -exponent)
        reference=1/np.sqrt(x.astype(float))
        ratio=np.abs(phase-reference)/(TOL*reference)
        ix=int(np.argmax(ratio))
        if ratio[ix]>worst:
            worst,worst_x=float(ratio[ix]),int(x[ix])
        failures+=int(np.count_nonzero(ratio>1))
    return dict(family=config.family,degree=config.degree,
                segment_bits=config.segment_bits,precision=q,
                drops=list(drops),checked_inputs=max_x,
                worst_error_over_tolerance=worst,worst_x=worst_x,
                failed_inputs=failures,passes=failures==0)
