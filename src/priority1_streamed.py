"""Low-workspace normalized direct phase: stream A, B, C from QROM separately.

Each coefficient is loaded, used, uncomputed and freed before the next.
This trades multiple passes through a coefficient table for lower peak
ancilla use. Uses the exact same phase function and dropped-product
semantics as priority1_truncated.build_truncated.
"""
from __future__ import annotations
from .hybrid_kickmix import _normalize, _denormalize, _phase, resource_counts
from .hybrid_model import NORMAL_BITS
from .priority1_model import Approximation, coefficient_matrix
from .priority1_exact import full_cut_screen
from .priority1_truncated import _mul_cut


def build_streamed(config: Approximation, drops=(10, 12, 11)):
    import kickmix as km
    if config.degree != 2:
        raise ValueError("Coefficient streaming currently requires a quadratic")
    dl, ds, dq = tuple(int(i) for i in drops)
    result = full_cut_screen(config, (dl,ds,dq))
    if result["failed_inputs"] or result["worst_error_over_tolerance"]>=.85:
        raise ValueError("Uncertified phase truncation")

    mat=coefficient_matrix(config)
    tables=[[int(v) for v in mat[:,degree]] for degree in range(3)]
    widths=[max(1,max(vals).bit_length()) for vals in tables]
    lo=config.low_bits
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(NORMAL_BITS,name="x")
    with builder.mark("normalize"):
        flags,shift=_normalize(builder,x)

    def lookup(degree):
        term=builder.alloc_qubits(widths[degree])
        builder.init_lookup(table=tables[degree],address=x[lo:],target=term)
        return term

    def release(degree,term):
        builder.del_lookup(table=tables[degree],address=x[lo:],target=term)
        builder.free(term)

    with builder.mark("constant"):
        a=lookup(0)
        _phase(builder,a,shift,config.precision,+1)
        release(0,a)

    with builder.mark("linear"):
        b=lookup(1)
        linear=builder.alloc_qubits(len(b)+lo-dl)
        builder.reset(linear)
        _mul_cut(builder,b,x[:lo],linear,dl)
        _phase(builder,linear,shift,config.precision+lo-dl,-1)
        _mul_cut(builder,b,x[:lo],linear,dl,inverse=True)
        builder.free(linear)
        release(1,b)

    with builder.mark("quadratic"):
        square=builder.alloc_qubits(2*lo-ds)
        builder.reset(square)
        _mul_cut(builder,x[:lo],x[:lo],square,ds)
        c=lookup(2)
        quad=builder.alloc_qubits(len(c)+len(square)-dq)
        builder.reset(quad)
        _mul_cut(builder,c,square,quad,dq)
        _phase(builder,quad,shift,config.precision+2*lo-ds-dq,+1)
        _mul_cut(builder,c,square,quad,dq,inverse=True)
        builder.free(quad)
        release(2,c)
        _mul_cut(builder,x[:lo],x[:lo],square,ds,inverse=True)
        builder.free(square)

    with builder.mark("restore_input"):
        _denormalize(builder,x,flags,shift)
    return builder.finish_circuit()
