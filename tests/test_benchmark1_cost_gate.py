"""No architecture claims may exceed independently supported cost floors."""
from src.benchmark1_cost_gate import MEASURED, TARGET_QC, CURRENT_VERIFIED_QC

def test_target_and_verified_thresholds():
    assert 69_000 < TARGET_QC < 70_000
    assert 465_000 < CURRENT_VERIFIED_QC < 466_000

def test_all_new_variants_rejected_by_30x_floor():
    for candidate in MEASURED:
        assert candidate.hard_cost_floor > TARGET_QC

def test_no_invalid_linear_candidate_reaches_full_verifier():
    assert not MEASURED[3].report()["eligible_for_full_verifier"]
    assert not MEASURED[4].report()["eligible_for_full_verifier"]

def test_worse_native_variants_cannot_beat_current_score():
    assert MEASURED[1].hard_cost_floor > CURRENT_VERIFIED_QC
    assert MEASURED[3].hard_cost_floor > CURRENT_VERIFIED_QC
    assert MEASURED[4].hard_cost_floor > CURRENT_VERIFIED_QC
