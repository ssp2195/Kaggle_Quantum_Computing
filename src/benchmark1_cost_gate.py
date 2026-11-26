"""Conservative architecture admission rules for Benchmark I.

Q * native Toffoli-like is a certified lower bound on the competition Q*C:
rotation synthesis, T gates, and stochastic correction costs are nonnegative.
Passing this *lower bound* is necessary, not sufficient. Candidates must also
pass native coherence, clean-scratch, and exact full-domain verification.
"""
from dataclasses import dataclass, asdict
import argparse
import json

GOOGLE_REFERENCE_COST = 2_076_880.0
TARGET_IMPROVEMENT = 30.0
TARGET_QC = GOOGLE_REFERENCE_COST / TARGET_IMPROVEMENT
CURRENT_VERIFIED_QC = 283_167.97477128694


@dataclass(frozen=True)
class Candidate:
    name: str
    qubits: int
    native_toffolis: int
    fine_rotations: int
    passed_native_sample: bool

    @property
    def hard_cost_floor(self):
        return self.qubits * self.native_toffolis

    def report(self):
        floor = self.hard_cost_floor
        return {
            **asdict(self),
            "q_times_native_toffoli_floor": floor,
            "can_reach_30x_under_this_gate_inventory": floor <= TARGET_QC,
            "can_beat_current_verified_under_this_gate_inventory":
                floor < CURRENT_VERIFIED_QC,
            "eligible_for_full_verifier": (
                floor < CURRENT_VERIFIED_QC and self.passed_native_sample),
        }


MEASURED = [
    Candidate("best-existing-grouped-barrel", 122, 3190, 88, True),
    Candidate("symmetric-square-ab-c", 138, 3892, 88, True),
    Candidate("row-square-ab-c", 138, 3190, 88, True),
    Candidate("piecewise-linear-joint", 135, 9308, 65, False),
    Candidate("piecewise-linear-streamed", 107, 17496, 65, False),
    Candidate("compact-mbuc-verified", 112, 3218, 88, True),
    Candidate("single-lookup-cut13-full-verified", 112, 2658, 88, True),
    Candidate("sparse-one-lookup-early-A-controlled-full-verified", 89, 2764, 87, True),
    Candidate("cost-aware-refitted-p7-sq13-swapped-full-verified", 86, 2652, 86, True),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = dict(target_qc=TARGET_QC,
                  best_current_verified_qc=CURRENT_VERIFIED_QC,
                  measured_candidates=[c.report() for c in MEASURED])
    print(json.dumps(result, indent=2))
    if args.output:
        from pathlib import Path
        p = Path(args.output)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
