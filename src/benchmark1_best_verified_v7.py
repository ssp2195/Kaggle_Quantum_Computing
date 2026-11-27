"""Pinned independently verified 8.097731x Benchmark I circuit (v7).

An original multioutput XOR/AND-network decoder uses unreachable normalized
mantissa addresses j=0..127 as Boolean don't-cares. The GF(2) extension
f00=f01 XOR f10 XOR f11 eliminates all high-pair interaction monomials
and 128 costly shared prefix ANDs. All reachable 384 coefficient rows are
EXACTLY the same as v6. Existing physically verified arithmetic, phase
padding, and binary exponent normalization remain unchanged.

Full competition-verifier public mirror PASS on 3,133,452 x 64 =
200,540,928 coherent quantum trajectories, independent of native builder.
Not an official Kaggle leaderboard submission.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from .benchmark1_dontcare_verify_candidate import construct
from .benchmark1_dontcare_boolean_rom import optimized_coefficients,proxy
from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256="bd691cb587bb51d690d86b4798b1a674fbb61313fa24d02f63ba6422ff4d413a"
VERIFIED_QC=256476.78036771005
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46014043594"
REFERENCE_QC=2076880.0
PREVIOUS_V6_QC=268172.78036771005
EXPECTED_RESOURCES=(86,2418,76)
EXPECTED_WORST_ERROR=0.9687625254666724
VERIFIED_SHOTS=200540928

def build_verified():
    c,stats=construct()
    if (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])!=EXPECTED_RESOURCES:
        raise RuntimeError("Pinned native resource mismatch: "+str(stats))
    mat,table=optimized_coefficients("affine")
    if not np.array_equal(mat[128:],cc.coefficients((10,13,12))[128:]):
        raise RuntimeError("Original reachable coefficients changed")
    if proxy(table)["top_interaction_monomials"]!=0:
        raise RuntimeError("Don't-care high-order Boolean elimination failed")
    return c

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_verified_v7.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    args=p.parse_args()
    c=build_verified()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(args.output,format="kmb")
    digest=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError("SHA-pinned v7 circuit mismatch: "+digest)
    report=dict(status="FULL_PINNED_MIRROR_PASS_NOT_OFFICIAL_KAGGLE",
       file=str(args.output),sha256=digest,
       verified_qc=VERIFIED_QC,reference_improvement=REFERENCE_QC/VERIFIED_QC,
       improvement_from_v6=(PREVIOUS_V6_QC-VERIFIED_QC)/PREVIOUS_V6_QC,
       native=resource_counts(c),
       worst_error_ratio=EXPECTED_WORST_ERROR,
       verified_shots=VERIFIED_SHOTS,full_verifier_run=VERIFIED_RUN)
    if args.verify_sample:
        from .hybrid_verify import verify_samples,boundary_inputs
        report["sample"]=verify_samples(c,boundary_inputs(random_count=args.verify_sample))
        if report["sample"]["failures"]:
            raise RuntimeError("Quantum coherent regression failed: "+str(report))
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
