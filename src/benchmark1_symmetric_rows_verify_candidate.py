"""Pinned full-verifier candidate: original carry-safe row-symmetric square.

Combines the independently full-verified v7 Boolean don't-care coefficient
ROM with independently authored row-compressed symmetric squares. All exact
partial products and phase gates are otherwise preserved.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_symmetric_square_rows_v7 import build,exact_math_check
from .hybrid_kickmix import resource_counts
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from . import benchmark1_cost_cooptim as cc
from .benchmark1_binary_block_normalizer import numeric_normalization_check

EXPECTED_RESOURCES=(86,2284,76)
BASELINE_QC=256476.78036771005
REFERENCE=2076880.0

def construct():
    assert exact_math_check()["passes"]
    assert numeric_normalization_check()["passes"]
    assert cc.full_screen((10,13,12))["passes"]
    mat,table=optimized_coefficients("affine")
    assert (mat[128:]==cc.coefficients((10,13,12))[128:]).all()
    assert proxy(table)["top_interaction_monomials"]==0
    c=build()
    s=resource_counts(c)
    assert (s["qubits"],s["toffoli_like"],s["fine_rotations"])==EXPECTED_RESOURCES,s
    return c,s

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/original_square_rows_v8.kmb"))
    p.add_argument("--samples",type=int,default=256)
    a=p.parse_args()
    c,s=construct()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(a.output,format="kmb")
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    report=dict(status="CANDIDATE_AWAITING_COMPLETE_INDEPENDENT_VERIFIER",
      native=s,sha256=digest,
      all_normalized_inputs=3145728,
      all_square_inputs=8192,
      legal_inputs=3133452,
      prior_verified_qc=BASELINE_QC,
      twenty_x_target_qc=REFERENCE/20)
    if a.samples:
        from .hybrid_verify import verify_samples,boundary_inputs
        report["coherent"]=verify_samples(c,boundary_inputs(random_count=a.samples))
        if report["coherent"]["failures"]:
            raise RuntimeError("Native quantum coherent sample failed: "+str(report))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
