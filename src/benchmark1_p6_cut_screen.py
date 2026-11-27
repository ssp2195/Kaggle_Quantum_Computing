"""New p6 (192 intervals) normalized quadratic exact-cut feasibility screen.

Current v6 has 384 intervals and 9 address bits; p6 uses 192 intervals
and only 8 address bits. The table lookup could be significantly cheaper
if the required 24-bit relative phase tolerance can be maintained.

This first screens approximate coefficients with EXACT omitted-product
integer arithmetic on the entire 3,133,452 allowed input domain; every
candidate for native synthesis requires an additional full normalized
input screen and independent quantum verifier.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .priority1_model import Approximation
from .priority1_exact import full_cut_screen

CONFIG=Approximation("chebyshev",2,6,26)
DROPS=((9,11,9),(9,12,9),(9,12,10),(10,11,9),
       (10,12,9),(10,12,10),(10,12,11),
       (10,13,9),(10,13,10),(10,13,11),
       (11,11,9),(11,12,9),(11,12,10),
       (11,13,9),(11,13,10),(11,13,11),
       (12,12,9),(12,12,10),(12,13,9),
       (12,13,10),(12,14,9))

def main():
    results=[]
    for cuts in DROPS:
        try:
            r=full_cut_screen(CONFIG,cuts)
            out=dict(cuts=cuts,**r)
        except Exception as e:
            out=dict(cuts=cuts,error=f"{type(e).__name__}: {str(e)}")
        results.append(out);print(json.dumps(out),flush=True)
    p=Path("runs/p6_original_cut_accuracy.json")
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(results,indent=2)+"\n")
if __name__=="__main__":main()
