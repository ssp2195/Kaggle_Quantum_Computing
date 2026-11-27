"""Exact candidate for original carry-free first-row v9, full independent verifier.

Preserves the original certified narrow-target row order. Each product
or square register begins at zero; the true first row can be written
without a carry adder, and the final reverse cleanup can erase it
without a subtractor. Native product phases and all unused ROM row
extensions remain identical to the full verified 85-qubit v8.

No Kaggle competitor circuit is used. Only the pinned independent
public contest verifier may establish a promoteable score.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
from .benchmark1_first_row_specialization import build
from .benchmark1_phase_pad9_symmetric import overflow_certificate
from .benchmark1_symmetric_square_rows_v7 import exact_math_check
from .benchmark1_binary_block_normalizer import numeric_normalization_check
from . import benchmark1_cost_cooptim as cc
from .hybrid_kickmix import resource_counts

EXPECTED=(85,2232,76)
SAMPLED_SHA="cb6c18f085c0c051cfa03856414c3f43037ad7abf8703a0a1be51ba4c44ac66c"
V8_QC=241424.49222389946
REF=2076880

def construct():
    for check in (overflow_certificate(),exact_math_check(),numeric_normalization_check()):
        assert check["passes"],check
    assert cc.full_screen((10,13,12))["passes"]
    c=build(True,True)
    s=resource_counts(c)
    if (s["qubits"],s["toffoli_like"],s["fine_rotations"])!=EXPECTED:
        raise RuntimeError("Changed compiler resources: "+str(s))
    return c,s

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=Path('out/firstrow_v9.kmb'))
    p.add_argument('--sample',type=int,default=256)
    a=p.parse_args()
    c,s=construct()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(a.out,format='kmb')
    sha=hashlib.sha256(a.out.read_bytes()).hexdigest()
    if sha!=SAMPLED_SHA:raise RuntimeError("KMB hash unexpected "+sha)
    report=dict(status="CANDIDATE_AWAITING_COMPLETE_VERIFIER",
         native=s,sha256=sha,
         normalized_words_checked=3145728,
         allowed_inputs=3133452,
         previously_verified_qc=V8_QC,
         twenty_x_qc_threshold=REF/20)
    if a.sample:
        from .hybrid_verify import boundary_inputs,verify_samples
        report["sample"]=verify_samples(c,boundary_inputs(random_count=a.sample))
        if report["sample"]["failures"]:raise RuntimeError(str(report["sample"]))
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
