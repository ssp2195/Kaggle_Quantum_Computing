"""Pinned independent Benchmark-I v5: exact reduced phase padding.

Verified with the pinned public mirror of the full competition verifier
across 3,133,452 inputs x 64 trajectories = 200,540,928. This is not a
Kaggle leaderboard submission. The original v4 remains reproducible.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_zero_residual_pad import build,full_pad_certificate,PADS,CUTS
from .benchmark1_cost_cooptim import full_screen
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256="4a232127f01b2a901209427e011ba512dd9fde991b627fbc2761399b17203ee4"
VERIFIED_QC=272300.78036771005
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45995733095"
REFERENCE_QC=2076880.0
EXPECTED_RESOURCES=(86,2602,76)
EXPECTED_WORST_ERROR=0.9687625254666724


def build_verified():
    cert=full_pad_certificate()
    assert cert["passes"] and (cert["min_linear_pad"],cert["min_quadratic_pad"])==(5,5),cert
    assert PADS==(10,5,5)
    approx=full_screen(CUTS)
    assert approx["passes"] and not approx["failures"],approx
    circuit,calls=build()
    stats=resource_counts(circuit)
    if (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])!=EXPECTED_RESOURCES:
        raise RuntimeError("Native resource mismatch: "+str(stats))
    return circuit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_verified_v5.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    a=p.parse_args()
    circ=build_verified()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    circ.to_file(a.output,format="kmb")
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError("Pinned v5 KMB hash mismatch: "+digest)
    result=dict(status="FULL_PINNED_MIRROR_PASS_NOT_KAGGLE_SUBMISSION",
      file=str(a.output),sha256=digest,verified_qc=VERIFIED_QC,
      improvement_vs_reference=REFERENCE_QC/VERIFIED_QC,
      resources=resource_counts(circ),worst_error_ratio=EXPECTED_WORST_ERROR,
      full_verifier_run=VERIFIED_RUN)
    if a.verify_sample:
        from .hybrid_verify import verify_samples,boundary_inputs
        result["coherent_sample"]=verify_samples(
          circ,boundary_inputs(random_count=a.verify_sample))
        if result["coherent_sample"]["failures"]:
            raise RuntimeError("Native coherence sample failed: "+str(result))
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
