"""Direct-phase quadratic circuit with rigorously budgeted partial-product truncation.

Unlike truncating a finished product, the reversible multiplier skips low
partial products entirely. Their maximum possible numerical contribution
is bounded independently of x. Does not assume ignored carries are zero.
"""
from __future__ import annotations
from .priority1_model import Approximation, pack_coefficients
from .hybrid_model import NORMAL_BITS, TOL
from .hybrid_kickmix import _normalize, _denormalize, _phase


def omitted_partial_products(n_a: int, n_b: int, drop: int) -> int:
    return sum(1 << (i+j)
               for i in range(n_a) for j in range(n_b)
               if i+j < drop)


def choose_drop(n_a: int, n_b: int, fractional_bits: int,
                scale: float = 1., tolerance_fraction: float = .025) -> int:
    """Uniform relative phase-error bound: 2**25 * error at normalized z."""
    drop=0
    for proposed in range(1, min(n_a+n_b-1,32)+1):
        bound=(
            2**25 * scale *
            omitted_partial_products(n_a,n_b,proposed) /
            2**fractional_bits
        )
        if bound > tolerance_fraction:
            break
        drop=proposed
    return drop


def plan(config: Approximation, tolerance_fraction=.025):
    if config.degree != 2:
        raise ValueError("Truncated implementation is quadratic only")
    table, widths = pack_coefficients(config)
    lo, q = config.low_bits, config.precision
    # Decode maximal integer quadratic coefficient before selecting square drop.
    max_c = max(v >> (widths[0]+widths[1]) for v in table)
    linear = choose_drop(widths[1],lo,q+lo,
                         tolerance_fraction=tolerance_fraction)
    square = choose_drop(lo,lo,2*lo,scale=max_c/2**q,
                         tolerance_fraction=tolerance_fraction)
    quadratic = choose_drop(widths[2],2*lo-square,
                            q+2*lo-square,
                            tolerance_fraction=tolerance_fraction)
    def budget(na,nb,drop,denom,scale=1):
        return 2**25 * omitted_partial_products(na,nb,drop) * scale /2**denom
    bounds = dict(linear=budget(widths[1],lo,linear,q+lo),
                  square=budget(lo,lo,square,2*lo,max_c/2**q),
                  quadratic=budget(widths[2],2*lo-square,
                                    quadratic,q+2*lo-square))
    return dict(linear_drop=linear,square_drop=square,
                quadratic_drop=quadratic,addition_error_bound=sum(bounds.values()),
                components=bounds,coefficient_widths=widths)


def _mul_cut(builder,a,b,dest,drop,inverse=False):
    """Compute sum_{i+j>=drop} a_i b_j 2**(i+j-drop), reversibly.

    This is *not* floor(a*b/2**drop): discarded low products could have
    carried. The omission is deliberately tracked in the error certificate.
    """
    na,nb=len(a),len(b)
    if len(dest)!=na+nb-drop:
        raise ValueError("Incorrect destination width")
    js=range(nb-1,-1,-1) if inverse else range(nb)
    for j in js:
        first=max(0,drop-j)
        if first>=na:
            continue
        row=builder.alloc_qubits(na-first+1)
        builder.reset(row)
        for i in range(first,na):
            if a[i]==b[j]:
                builder.cx(a[i],row[i-first])
            else:
                builder.init_and(a[i],b[j],row[i-first])
        start=j+first-drop
        if inverse:
            builder.isub(row,target=dest[start:start+len(row)])
        else:
            builder.iadd(row,target=dest[start:start+len(row)])
        for i in reversed(range(first,na)):
            if a[i]==b[j]:
                builder.cx(a[i],row[i-first])
            else:
                builder.del_and(a[i],b[j],row[i-first])
        builder.free(row)


def build_truncated(config: Approximation, *, tolerance_fraction=.025):
    import kickmix as km
    parameters=plan(config,tolerance_fraction)
    dl,ds,dq=(parameters[k] for k in
              ("linear_drop","square_drop","quadratic_drop"))
    table,widths=pack_coefficients(config)
    lo=config.low_bits
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(NORMAL_BITS,name="x")
    with builder.mark("normalize"):
        flags,shift=_normalize(builder,x)
    with builder.mark("coefficient_lookup"):
        coeff=builder.alloc_qubits(sum(widths))
        builder.init_lookup(table=table,address=x[lo:],target=coeff)
        a=coeff[:widths[0]]
        b=coeff[widths[0]:widths[0]+widths[1]]
        c=coeff[widths[0]+widths[1]:]

    with builder.mark("constant_phase"):
        _phase(builder,a,shift,config.precision,+1)

    with builder.mark("truncated_linear"):
        linear=builder.alloc_qubits(len(b)+lo-dl)
        builder.reset(linear)
        _mul_cut(builder,b,x[:lo],linear,dl)
        _phase(builder,linear,shift,config.precision+lo-dl,-1)
        _mul_cut(builder,b,x[:lo],linear,dl,inverse=True)
        builder.free(linear)

    with builder.mark("truncated_quadratic"):
        square=builder.alloc_qubits(2*lo-ds)
        builder.reset(square)
        _mul_cut(builder,x[:lo],x[:lo],square,ds)
        quadratic=builder.alloc_qubits(len(c)+len(square)-dq)
        builder.reset(quadratic)
        _mul_cut(builder,c,square,quadratic,dq)
        _phase(builder,quadratic,shift,
               config.precision+2*lo-ds-dq,+1)
        _mul_cut(builder,c,square,quadratic,dq,inverse=True)
        builder.free(quadratic)
        _mul_cut(builder,x[:lo],x[:lo],square,ds,inverse=True)
        builder.free(square)

    with builder.mark("cleanup"):
        builder.del_lookup(table=table,address=x[lo:],target=coeff)
        builder.free(coeff)
        _denormalize(builder,x,flags,shift)
    return builder.finish_circuit()
