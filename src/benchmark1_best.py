"""Rebuild the lowest-cost fully verified Benchmark I kickmix circuit.

The pinned public mirror of the competition's full inverse-square-root verifier
accepted 3,133,452 inputs × 64 shots on 2026-11-09 (200,540,928 trajectories).
The score is for Benchmark I alone, not an official Kaggle leaderboard score.

Architecture: Chebyshev quadratic 7/26, cut (10, 12, 11), grouped A then BC,
compact measurement-corrected normalization.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from .priority1_grouped import build_grouped
from .priority1_model import Approximation
from .hybrid_kickmix import resource_counts

BEST_CONFIG = Approximation("chebyshev", 2, 7, 26)
DROPS = (10, 12, 11)
VERIFIED_QC = 430_259.61862786254
VERIFIED_TOFFOLI_EQUIVALENT = 3_841.603737748773
VERIFIED_SHA256 = "938d053f2218caa1859837ffb466ded01033c2d8d6445699276e5b4e21ebbc42"
VERIFICATION_RUN = (
    "https://github.com/ssp2195/Kaggle_Quantum_Computing/"
    "actions/runs/42657329801"
)

def build():
    """Construct the complete 22-input-qubit native circuit."""
    return build_grouped(BEST_CONFIG, DROPS,
                         grouping="a_bc", normalization="compact")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("out/best_invsqrt.kmb"))
    a = p.parse_args()
    circuit = build()
    stats = resource_counts(circuit)
    expected = (112, 3218, 88)
    observed = (stats["qubits"], stats["toffoli_like"], stats["fine_rotations"])
    if observed != expected:
        raise RuntimeError(
            f"Native kickmix resource counts changed from {expected} to {observed}. "
            "The rebuilt circuit requires full reverification."
        )
    a.output.parent.mkdir(parents=True, exist_ok=True)
    circuit.to_file(a.output, format="kmb")
    digest = hashlib.sha256(a.output.read_bytes()).hexdigest()
    if digest != VERIFIED_SHA256:
        raise RuntimeError(
            f"Uncertified KMB output: SHA256 {digest} differs from full-verified "
            f"{VERIFIED_SHA256}. The generated file was NOT certified as equivalent; "
            "use the pinned build or run complete verification again."
        )
    print(json.dumps({
        "file": str(a.output), "sha256": digest, "resource_counts": stats,
        "verified_benchmark1_q_times_c": VERIFIED_QC,
        "verified_mean_toffoli_equivalent": VERIFIED_TOFFOLI_EQUIVALENT,
        "verification": VERIFICATION_RUN,
        "verified_shots": 200_540_928,
        "note": "Full independent mirror-verifier PASS; not a new Kaggle submission."
    }, indent=2))

if __name__ == "__main__":
    main()
