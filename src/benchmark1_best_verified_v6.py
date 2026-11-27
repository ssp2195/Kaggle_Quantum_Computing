"""Pinned and independently verified v6 inverse-sqrt oracle.

Original four-stage base-four leading-zero normalization replaces the old
repeated binary exponent increments, retains the v5 exact five-wire phase
padding and independently certified p7/q26 truncated arithmetic. The v5
builder is untouched and kept as a reproducible fallback.

Full pinned public verifier mirror: 3,133,452 valid inputs * 64 = 200,540,928
verified trajectories. This is NOT an official Kaggle leaderboard result.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from .benchmark1_binary_block_normalizer import build,numeric_normalization_check
from .benchmark1_zero_residual_pad import full_pad_certificate
from .benchmark1_cost_cooptim import full_screen
from .hybrid_kickmix import resource_counts

VERIFIED_SHA256="1807c84b9f93008de4ebde2c2231aa8d5f08386dc03e10b63b47999c4577011d"
VERIFIED_QC=268172.78036771005
VERIFIED_RUN="https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46008864751"
REFERENCE_QC=2076880.0
EXPECTED_RESOURCES=(86,2554,76)
EXPECTED_WORST_ERROR=0.9687625254666724
EXPECTED_SHOTS=200540928

def build_verified():
    norm=numeric_normalization_check()
    if not norm["passes"] or norm["checked"]!=3133452:
        raise RuntimeError("Exact binary-normalization full-domain screen failed: "+str(norm))
    pad=full_pad_certificate()
    if not pad["passes"] or (pad["min_linear_pad"],pad["min_quadratic_pad"])!=(5,5):
        raise RuntimeError("Missing exact five-qubit padding proof: "+str(pad))
    approx=full_screen((10,13,12))
    if not approx["passes"] or approx["failures"]:
        raise RuntimeError("Exact integer polynomial has errors: "+str(approx))
    c,calls=build()
    if len(calls)!=3:raise RuntimeError("Unexpected phase stage count: "+str(calls))
    s=resource_counts(c)
    if (s["qubits"],s["toffoli_like"],s["fine_rotations"])!=EXPECTED_RESOURCES:
        raise RuntimeError("Native resource regression: "+str(s))
    return c

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,default=Path("out/best_verified_v6.kmb"))
    p.add_argument("--verify-sample",type=int,default=0)
    args=p.parse_args()
    circuit=build_verified()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(args.output,format="kmb")
    digest=hashlib.sha256(args.output.read_bytes()).hexdigest()
    if digest!=VERIFIED_SHA256:
        raise RuntimeError("Pinned v6 .kmb digest mismatch: "+digest)
    out=dict(status="FULL_PINNED_MIRROR_PASS_NOT_OFFICIAL_KAGGLE",
       file=str(args.output),sha256=digest,
       verified_qc=VERIFIED_QC,
       improvement_over_reference=REFERENCE_QC/VERIFIED_QC,
       resources=resource_counts(circuit),
       worst_error_ratio=EXPECTED_WORST_ERROR,
       verified_shots=EXPECTED_SHOTS,
       verifier_run=VERIFIED_RUN)
    if args.verify_sample:
        from .hybrid_verify import verify_samples,boundary_inputs
        sample=verify_samples(circuit,boundary_inputs(random_count=args.verify_sample))
        if sample["failures"]:raise RuntimeError("Native coherence failure: "+str(sample))
        out["sample"]=sample
    print(json.dumps(out,indent=2))

if __name__=="__main__":main()
