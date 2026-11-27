"""Best independently full-verified Benchmark 1 inverse-square-root phase oracle v9.

Original exact carry-free first-row arithmetic for both truncated
variable-coefficient multipliers and row-compressed symmetric squaring.
The first row is assigned into an initially-zero accumulator directly;
every following row keeps the previously certified carry-safe ordering.
The inverse undoes rows in precisely reverse order, erasing the first
directly when the accumulator again equals the initial first row.

Retains separately certified four-stage normalization, GF(2) unreachable
QROM don't-care ANF decoder, nine-qubit intercept phase, and
measurement-corrected cleanup. Every reachable coefficient and effective
phase is bit-for-bit identical to verified v8.

Full independent pinned verifier: PASS on 200,540,928 trajectories.
This is NOT an official Kaggle leaderboard submission.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_first_row_specialization import build
from .benchmark1_phase_pad9_symmetric import overflow_certificate
from .benchmark1_symmetric_square_rows_v7 import exact_math_check
from .benchmark1_binary_block_normalizer import numeric_normalization_check
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256="cb6c18f085c0c051cfa03856414c3f43037ad7abf8703a0a1be51ba4c44ac66c"
VERIFIED_QC=237684.49222389946
REFERENCE_QC=2076880.0
PRIOR_V8_QC=241424.49222389946
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46021681763"
EXPECTED_RESOURCES=(85,2232,76)
EXPECTED_VERIFIED_SHOTS=200540928
WORST_ERROR_RATIO=0.9687625254666724

def build_verified():
    if not overflow_certificate()["passes"]:raise RuntimeError("Intercept overflow proof failed")
    if not exact_math_check()["passes"]:raise RuntimeError("Square exact dropped-column equality failed")
    if not numeric_normalization_check()["passes"]:raise RuntimeError("Normalization full legal domain failed")
    if not cc.full_screen((10,13,12))["passes"]:raise RuntimeError("Phase math normalized domain failed")
    table,packed=optimized_coefficients("affine")
    if not (table[128:]==cc.coefficients((10,13,12))[128:]).all():raise RuntimeError("Reachable coeff changed")
    if proxy(packed)["top_interaction_monomials"]:raise RuntimeError("ANF decoder Boolean topology regressed")
    c=build(True,True)
    s=resource_counts(c)
    measured=(s["qubits"],s["toffoli_like"],s["fine_rotations"])
    if measured!=EXPECTED_RESOURCES:raise RuntimeError(f"Native resource regression: {measured}")
    return c

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,default=Path("out/best_verified_v9.kmb"))
    ap.add_argument("--verify-sample",type=int,default=0)
    args=ap.parse_args()
    c=build_verified()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    sha=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if sha!=VERIFIED_SHA256:raise RuntimeError(f"Pinned v9 SHA mismatch: {sha}")
    r=dict(status="INDEPENDENT_PUBLIC_FULL_VERIFIER_PASS_NOT_OFFICIAL_KAGGLE",
        file=str(args.output),sha256=sha,verified_qc=VERIFIED_QC,
        reference_improvement=REFERENCE_QC/VERIFIED_QC,
        previous_v8_qc=PRIOR_V8_QC,resources=resource_counts(c),
        verified_shots=EXPECTED_VERIFIED_SHOTS,
        worst_error_ratio=WORST_ERROR_RATIO,
        verifier_run=VERIFIED_RUN)
    if args.verify_sample:
        from .hybrid_verify import boundary_inputs,verify_samples
        r["sample"]=verify_samples(c,boundary_inputs(random_count=args.verify_sample))
        if r["sample"]["failures"]:raise RuntimeError("Coherent sample mismatch: "+str(r["sample"]))
    print(json.dumps(r,indent=2))

if __name__=="__main__":main()
