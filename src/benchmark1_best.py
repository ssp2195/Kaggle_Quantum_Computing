"""Build the best independently verified Benchmark I circuit.

Pinned public competition-verifier mirror:
Q=122, mean Toffoli-equivalent=3813.603737748773,
Q*C=465259.6560053503, PASS on 200540928 trajectories.
This is not a Kaggle combined leaderboard score.
"""
import argparse
import json
from pathlib import Path

from .priority1_grouped import build_grouped
from .priority1_model import Approximation
from .hybrid_kickmix import resource_counts

BEST_CONFIG = Approximation("chebyshev", 2, 7, 26)
DROPS = (10, 12, 11)
VERIFIED_QC = 465_259.6560053503
VERIFICATION_RUN = (
    "https://github.com/ssp2195/Kaggle_Quantum_Computing/"
    "actions/runs/42632214980"
)

def build():
    return build_grouped(BEST_CONFIG, DROPS,
                         grouping="a_bc", normalization="barrel")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, default=Path("out/best_invsqrt.kmb"))
    a = p.parse_args()
    circuit = build()
    stats = resource_counts(circuit)
    if (stats["qubits"], stats["toffoli_like"]) != (122, 3190):
        raise RuntimeError(
            "Kickmix build changed resource counts; repeat full verification: "
            + str(stats))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    circuit.to_file(a.output, format="kmb")
    print(json.dumps(dict(
        file=str(a.output), resources=stats,
        previously_verified_qc=VERIFIED_QC,
        evidence=VERIFICATION_RUN,
        warning="Repeat full verifier when kickmix version or output bytes change."),
        indent=2))

if __name__ == "__main__":
    main()
