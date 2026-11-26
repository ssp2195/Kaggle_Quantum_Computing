"""Pinned best verified Benchmark-I circuit, Nov 26 2026.

This reconstructs exactly the native KMB accepted by the independent pinned
full 200,540,928-trajectory verifier. No Kaggle leaderboard claim.
"""
from __future__ import annotations
from pathlib import Path
import argparse,hashlib,json
from .benchmark1_single_joint_lookup import build
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256 = "0589741be529205d88ebe67d4d9d6a82faebb42debfaf1dcf38c0e39383d09b6"
VERIFIED_COST = 302641.04424652964
VERIFIED_RUN = "https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45833905951"

def build_verified():
    return build(controlled=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_joint_v3.kmb"))
    a=p.parse_args()
    circuit=build_verified()
    stats=resource_counts(circuit)
    if (stats["qubits"],stats["toffoli_like"],stats["fine_rotations"])!=(89,2764,87):
        raise RuntimeError("Native resources differ from verified full scorer")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(a.output,format="kmb")
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError("Generated KMB bytes differ from pinned full-verifier PASS: "+digest)
    print(json.dumps({"file":str(a.output),"sha256":digest,"resources":stats,
          "full_verified_qc":VERIFIED_COST,"full_verifier_run":VERIFIED_RUN,
          "status":"PINNED_MIRROR_PASS_NOT_KAGGLE_SUBMISSION"},indent=2))

if __name__=="__main__":main()
