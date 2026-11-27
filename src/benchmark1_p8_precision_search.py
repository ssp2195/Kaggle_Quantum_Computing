"""Original p8/p9 segmentation and precision exploration.

Search a distinct algorithmic regime: 768 or 1536 normalized intervals,
12/11 mantissa residual bits, reduced quadratic product widths and
potential 25-bit coefficient precision. Every candidate uses the exact
omitted-column integer model over all 3,133,452 legal inputs. This is
pure mathematical feasibility, not a quantum circuit or verified score.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .priority1_model import Approximation
from .priority1_exact import full_cut_screen

def generate():
    for p in (8,9):
        for q in (25,26):
            cuts=(
              [(7,10,9),(8,10,9),(8,10,10),(8,11,9),(8,11,10),
               (9,10,9),(9,11,9),(9,11,10),(9,12,9),
               (9,12,10),(10,11,9),(10,11,10),(10,12,10)]
              if p==8 else
              [(7,8,8),(7,9,8),(8,9,8),(8,10,8),
               (8,10,9),(9,9,8),(9,10,8),(9,10,9)])
            for v in cuts:yield p,q,v

def main():
    rows=[]
    for p,q,drops in generate():
        cfg=Approximation("chebyshev",2,p,q)
        try:
            r=full_cut_screen(cfg,drops)
            out=dict(segment_bits=p,precision=q,cuts=list(drops),
                failures=r["failed_inputs"],worst_error_ratio=r["worst_error_over_tolerance"],
                passes=r["passes"],legal_inputs=r["checked_inputs"])
        except Exception as e:
            out=dict(segment_bits=p,precision=q,cuts=list(drops),
                     error=type(e).__name__+":"+str(e))
        print(json.dumps(out),flush=True)
        rows.append(out)
    path=Path("runs/p8_p9_accuracy_screen.json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,indent=2)+'\n')
    valid=[x for x in rows if x.get("passes")]
    print("SUMMARY "+json.dumps(dict(passing=len(valid),min_errors=sorted(valid,key=lambda x:x["worst_error_ratio"])[:8])),flush=True)
if __name__=="__main__":main()
