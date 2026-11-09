"""Native kickmix implementations of quantized quadratic/cubic direct phases.

Uses the *same* normalization and low-level arithmetic as the verified
8/26 Taylor baseline so resource differences can be attributed to the
mathematical approximation. Alternative architectures can be added
without changing the exact verifier contract.
"""
from __future__ import annotations
from .priority1_model import Approximation, pack_coefficients
from .hybrid_kickmix import _normalize, _denormalize, _phase, _mul, resource_counts
from .hybrid_model import NORMAL_BITS


def build_priority1(config: Approximation):
    import kickmix as km
    table, widths = pack_coefficients(config)
    low_bits = config.low_bits
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(NORMAL_BITS,name="x")
    with builder.mark("normalize"):
        flags,shift=_normalize(builder,x)
    with builder.mark("coefficient_lookup"):
        packed=builder.alloc_qubits(sum(widths))
        builder.init_lookup(table=table,address=x[low_bits:],target=packed)
        coeff=[]
        start=0
        for w in widths:
            coeff.append(packed[start:start+w])
            start+=w

    with builder.mark("phase_constant"):
        _phase(builder,coeff[0],shift,config.precision,+1)

    # Successively compute u^2 and u^3; hold only one monomial power
    # at a time, then uncompute it after imprinting its phase.
    powers=[]
    for degree in range(1,config.degree+1):
        source = x[:low_bits] if degree==1 else powers[-1]
        if degree > 1:
            power=builder.alloc_qubits(degree*low_bits)
            builder.reset(power)
            _mul(builder,source,x[:low_bits],power)
            powers.append(power)
            source=power
        with builder.mark("phase_degree_"+str(degree)):
            result=builder.alloc_qubits(len(coeff[degree])+len(source))
            builder.reset(result)
            _mul(builder,coeff[degree],source,result)
            _phase(builder,result,shift,
                   config.precision+degree*low_bits,
                   (-1)**degree)
            _mul(builder,coeff[degree],source,result,inverse=True)
            builder.free(result)

    for degree in reversed(range(2,config.degree+1)):
        power=powers.pop()
        previous=x[:low_bits] if degree==2 else powers[-1]
        _mul(builder,previous,x[:low_bits],power,inverse=True)
        builder.free(power)
    with builder.mark("cleanup"):
        builder.del_lookup(table=table,address=x[low_bits:],target=packed)
        builder.free(packed)
        _denormalize(builder,x,flags,shift)
    return builder.finish_circuit()
