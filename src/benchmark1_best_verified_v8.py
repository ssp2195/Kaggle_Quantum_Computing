"""Independently full-verified 85-qubit Benchmark I winning circuit (v8).

Original carry-safe row-compressed symmetric-square arithmetic shares carry
propagation by building and adding one clean temporary row per source bit.
Combined with an original GF(2) don't-care coefficient ROM and the four-stage
greedy base-four normalizer, it has exact controlled phase behavior on all
valid inputs.

The intercept phase discards a single 2*pi whole-turn contribution at x=1,
permitting a 9-qubit rather than 10-qubit padding register. The omitted turn
is a global identity; accuracy otherwise unchanged.

Full pinned public verifier mirror PASS on all 3,133,452 legal inputs times
64 salted coherent trajectories, NOT an official Kaggle leaderboard score.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_phase_pad9_symmetric import build,overflow_certificate
from .benchmark1_symmetric_square_rows_v7 import exact_math_check
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from .benchmark1_binary_block_normalizer import numeric_normalization_check
from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256="dd3361737c5292f4eecf90ae79837659b3b434ff6beae589268d5b32a28ab774"
VERIFIED_QC=241424.49222389946
REFERENCE_QC=2076880.0
PRIOR_V7_QC=256476.78036771005
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46017040842"
EXPECTED_RESOURCES=(85,2276,76)
EXPECTED_VERIFIED_SHOTS=200540928
WORST_ERROR_RATIO=0.9687625254666724

def build_verified():
    if not overflow_certificate()["passes"]:raise RuntimeError("Whole-turn overflow is unsafe")
    if not exact_math_check()["passes"]:raise RuntimeError("Square exact integer equality failed")
    if not numeric_normalization_check()["passes"]:raise RuntimeError("Normalization full domain failed")
    if not cc.full_screen((10,13,12))["passes"]:raise RuntimeError("Polynomial full numerical domain failed")
    co,packed=optimized_coefficients("affine")
    if not (co[128:]==cc.coefficients((10,13,12))[128:]).all():
        raise RuntimeError("Reachable coefficient modifications detected")
    if proxy(packed)["top_interaction_monomials"]:raise RuntimeError("ANF sparsity regression")
    c=build()
    s=resource_counts(c)
    actual=(s["qubits"],s["toffoli_like"],s["fine_rotations"])
    if actual!=EXPECTED_RESOURCES:
        raise RuntimeError(f"Native resource regression: {actual} vs {EXPECTED_RESOURCES}")
    return c

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_verified_v8.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    args=p.parse_args()
    c=build_verified()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    digest=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError(f"SHA mismatch: {digest} vs {VERIFIED_SHA256}")
    report=dict(status="PINNED_FULL_MIRROR_PASS_NOT_OFFICIAL_KAGGLE",
       sha256=digest,output=str(args.output),
       native=resource_counts(c),verified_qc=VERIFIED_QC,
       improvement=REFERENCE_QC/VERIFIED_QC,
       previous_verified_qc=PRIOR_V7_QC,
       full_verifier_run=VERIFIED_RUN,
       worst_error_ratio=WORST_ERROR_RATIO,
       verified_shots=EXPECTED_VERIFIED_SHOTS)
    if args.verify_sample:
        from .hybrid_verify import boundary_inputs,verify_samples
        sample=verify_samples(c,boundary_inputs(random_count=args.verify_sample))
        if sample["failures"]:raise RuntimeError(str(sample))
        report["sample"]=sample
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
