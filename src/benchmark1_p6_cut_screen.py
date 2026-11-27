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
DROPS=((10,14,13),(11,14,13),(11,15,13),(11,15,14),
       (11,16,13),(11,16,14),(12,15,13),(12,15,14),
       (12,16,13),(12,16,14),(12,17,14),(12,17,15),
       (11,14,14),(12,14,14),(13,16,14))

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
