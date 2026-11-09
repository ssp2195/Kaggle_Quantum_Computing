"""Assess measured phase-gadget streaming against current kickmix QC budget.

This computes a *cost proxy*, not a validated 22-bit inverse-square-root
oracle: distributing the nested truncated square/product alters carry
semantics. A complete candidate would need full-domain accuracy screening
and coherent verifier acceptance before any score could be claimed.

The proxy counts all direct (coefficient bit, mantissa bit) monomials
above the existing linear cut, and all independent (coefficient bit,
mantissa pair) monomials above the two quadratic cuts. It does not
pretend that their sum exactly equals the nested reversible truncation
scheme, or that all their rotation angles must be distinct.
"""
import argparse
import json
from math import log2

from .benchmark1_cost_gate import TARGET_QC, CURRENT_VERIFIED_QC

SEGMENT_BITS = 7
LOW_BITS = 13
COEFFICIENT_WIDTHS = (27, 18, 11)
LINEAR_DROP = 10
SQUARE_DROP = 12
QUADRATIC_DROP = 11
NATIVE_ROTATIONS = 88
NATIVE_TOFFOLIS = 3218
NATIVE_Q = 112
VERIFIED_TOFF_EQUIVALENT = 3841.603737748773
REFERENCE = 2_076_880


def raw_phase_term_counts():
    a, b, c = COEFFICIENT_WIDTHS
    low = LOW_BITS
    linear = [(i, j) for i in range(b) for j in range(low)
              if i + j >= LINEAR_DROP]
    sq_pairs = [(j, k) for j in range(low) for k in range(j, low)
                if j + k >= SQUARE_DROP]
    quadratic = [(i, j, k)
                 for i in range(c)
                 for j, k in sq_pairs
                 if i + j + k - SQUARE_DROP >= QUADRATIC_DROP]
    cubic = [x for x in quadratic if x[1] != x[2]]
    repeated = [x for x in quadratic if x[1] == x[2]]
    return dict(
        intercept_bits=a,
        linear_bit_pair_terms=len(linear),
        square_bit_pairs=len(sq_pairs),
        quadratic_unique_terms=len(quadratic),
        quadratic_crossbit_degree3_terms=len(cubic),
        quadratic_squarebit_degree2_terms=len(repeated),
        total_nonconstant_raw_phase_terms=len(linear) + len(quadratic),
    )


def assess():
    terms = raw_phase_term_counts()
    n = terms["total_nonconstant_raw_phase_terms"]
    # Actual measured per-rotation contribution in the current verifier,
    # NOT a universal per-angle lower bound.
    eq_per_rotation = (VERIFIED_TOFF_EQUIVALENT-NATIVE_TOFFOLIS-.25)/NATIVE_ROTATIONS
    per_term = dict()
    for assumed_qubits in (22, 40, 50, 63, 112):
        per_term[str(assumed_qubits)] = dict(
            rotations_only_cost_if_one_fine_angle_per_raw_term=
                assumed_qubits*n*eq_per_rotation,
            allowable_equivalent_operations_at_30x=TARGET_QC/assumed_qubits,
            max_fine_rotations_if_zero_Toffolis_and_same_rate=
                int(TARGET_QC/(assumed_qubits*eq_per_rotation)),
        )
    return dict(
        status="cost proxy only; not a complete valid phase oracle",
        approximation_basis="p7 q26 quadratic, raw direct-product monomial counts",
        **terms,
        current_verified=dict(qubits=NATIVE_Q,toffolis=NATIVE_TOFFOLIS,
                              fine_rotations=NATIVE_ROTATIONS,
                              qc=430259.61862786254),
        current_effective_non_Toffoli_equivalent_per_fine_rotation=
            eq_per_rotation,
        estimated_rotation_cost_only=per_term,
        naive_exponent_scaling_rotations_if_11_per_raw_term=n*11,
        target_qc=TARGET_QC,
        important_caveats=[
            "Nested cut-products and distributed monomial cuts are not exactly identical.",
            "No common-subexpression sharing, lookups, normalization, or prefix resources counted.",
            "Rotation cost changes with precision budget and the number of angles.",
            "Some terms can be merged or omitted only after full-domain mathematical screening.",
            "Term counts alone do not establish a mathematical lower bound on all architectures.",
        ],
    )


def native_single_term(mode="inline"):
    """Native resource proxy for one variable-exponent quadratic phase term."""
    import kickmix as km
    from .hybrid_kickmix import _phase, resource_counts
    if mode not in ("inline","mbuc"):
        raise ValueError("mode must be inline or mbuc")
    b=km.CircuitBuilder()
    x=b.create_quantum_register(2,name="x")
    exponent=b.create_quantum_register(4,name="exponent")
    t=b.alloc_qubits(1)
    b.init_and(x[0],x[1],t[0])
    _phase(b,t,exponent,frac_bits=20,sign=+1)
    if mode=="mbuc":
        b.del_and(x[0],x[1],t[0])
    else:
        result=b.alloc_bits(1)[0]
        b.hmr(t[0],result)
        b.ccz(x[0],x[1],result)
        b.free(result)
    b.free(t)
    return resource_counts(b.finish_circuit())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native",action="store_true")
    p.add_argument("--output")
    args=p.parse_args()
    data=assess()
    if args.native:
        data["native_one_shifted_term"]={
            "inline":native_single_term("inline"),
            "mbuc":native_single_term("mbuc")
        }
    print(json.dumps(data,indent=2))
    if args.output:
        from pathlib import Path
        path=Path(args.output)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(data,indent=2)+"\n")


if __name__=="__main__":
    main()
