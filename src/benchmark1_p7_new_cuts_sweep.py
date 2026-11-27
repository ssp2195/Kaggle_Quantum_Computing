"""Original exact-cut co-optimization sweep beyond the v8 cut pattern.

Each candidate jointly refits intercept, slope, and curvature
coefficients against EXACT omitted partial-product columns. A stricter
complete normalized-domain certificate tests every 3*128*8192 word
rather than sampling the physical 3,133,452 input values.

No native circuit or Q*C improvement is claimed from a numerical
screen alone. Promotable cases must subsequently pass native resource,
coherent sample and the pinned 200,540,928-trajectory verifier.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from . import benchmark1_cost_cooptim_exp as exp

CONFIGS=((11,13,12),(10,14,12),(10,13,13),(11,14,13),
         (11,14,12),(11,13,13),(10,14,13),(12,13,12))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max',type=int,default=8)
    a=ap.parse_args()
    output=[]
    for cuts in CONFIGS[:a.max]:
        try:
            report=exp.full_screen(cuts)
            record=dict(cuts=cuts,**report)
        except Exception as e:
            record=dict(cuts=cuts,error=type(e).__name__+":"+str(e)[:600])
        output.append(record)
        print(json.dumps(record),flush=True)
    p=Path('runs/p7_joint_refit_new_cuts.json')
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':main()
