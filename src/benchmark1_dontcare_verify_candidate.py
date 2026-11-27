"""Independent full-verifier candidate: Boolean ROM don't-care optimization.

Use 384 unchanged reachable p7 coefficients, but synthesize all 128
unreachable j<128 entries as XOR f01^f10^f11, canceling every
h7*h8 ANF interaction. The original optimized polynomial and binary
normalization (v6) remain bit-for-bit unchanged on all allowed inputs.

Sampled quantum phase was PASS. Full pinned mirror required to promote.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import build,optimized_coefficients,proxy
from .benchmark1_binary_block_normalizer import numeric_normalization_check
from .benchmark1_zero_residual_pad import full_pad_certificate
from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts

BASE_SCORE=268172.78036771005
EXPECTED=(86,2418,76)
REFERENCE_QC=2076880
CUTS=(10,13,12)

def construct():
    assert numeric_normalization_check()["passes"]
    assert full_pad_certificate()["passes"]
    report=cc.full_screen(CUTS)
    assert report["passes"] and not report["failures"]
    mat,table=optimized_coefficients("affine")
    expected=cc.coefficients(CUTS)
    assert (mat[128:]==expected[128:]).all()
    assert proxy(table)["top_interaction_monomials"]==0
    circuit,_=build("affine")
    stats=resource_counts(circuit)
    if (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])!=EXPECTED:
        raise RuntimeError("Unexpected independently compiled ROM costs: "+str(stats))
    return circuit,stats

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,default=Path("out/dontcare_anf_v7.kmb"))
    p.add_argument("--sample",type=int,default=256)
    a=p.parse_args()
    c,stats=construct()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(a.output,format="kmb")
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    report=dict(status="CANDIDATE_AWAITING_FULL_VERIFIER",
        math_reachable_coefficients_unchanged=True,
        normalized_inputs_checked=3145728,legal_inputs=3133452,
        native=stats,sha256=digest,
        cost_lower_bound=stats["lower_bound_qc"],
        prior_full_verified_qc=BASE_SCORE,
        target_20x_qc=REFERENCE_QC/20)
    if a.sample:
        from .hybrid_verify import verify_samples,boundary_inputs
        report["coherent"]=verify_samples(c,boundary_inputs(random_count=a.sample))
        if report["coherent"]["failures"]:
            raise RuntimeError("Quantum coherent sample failed: "+str(report))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
