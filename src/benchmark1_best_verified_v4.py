"""Exact pinned v4: cost-aware p7/q26 phase oracle; public full verifier PASS.

Reproduction checks byte-identical KMB SHA-256 and native gate counts.
No assertion of an official Kaggle leaderboard submission.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_cost_cooptim import build,full_screen
from .hybrid_kickmix import resource_counts

DROPS=(10,13,12)
VERIFIED_SHA256="0ab11da9ee2ff0b996cdd9bbc89eb1716d2c88bd1f1d934dab3695f7d566ed30"
VERIFIED_QC=283167.97477128694
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45857957278"

def build_verified():
    cert=full_screen(DROPS)
    if not cert["passes"] or cert["worst_error_over_tolerance"]>=1:
        raise ValueError("All-normalized-input numerical gate rejected: "+str(cert))
    return build(DROPS,swap_quadratic=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_verified_v4.kmb"))
    a=p.parse_args()
    circuit=build_verified()
    stats=resource_counts(circuit)
    if (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])!=(86,2652,86):
        raise RuntimeError("Native resource mismatch against full verifier")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(a.output,format="kmb")
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError("Pinned KMB fingerprint mismatch: "+digest)
    print(json.dumps(dict(file=str(a.output),sha256=digest,resources=stats,
        verified_qc=VERIFIED_QC,full_verifier_run=VERIFIED_RUN,
        status="FULL_PINNED_MIRROR_PASS_NOT_KAGGLE_SUBMISSION"),indent=2))

if __name__=="__main__":main()
