"""Find minimum safe padding for each residual phase over the full input set.

The result is a rigorous integer no-wrap certificate, not a numerical phase
approximation or a sampled argument. No competing synthesis is imported.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product
from .hybrid_model import MAX_X

def analyze(chunk=200_000):
    coeff=cc.coefficients((10,13,12))
    b=coeff[:,1]>>cc.BZ
    c=coeff[:,2]>>cc.CZ
    wb=max(1,int(b.max()).bit_length())
    wc=max(1,int(c.max()).bit_length())
    lb=wb+cc.LO-(10-cc.BZ)
    lc=wc+(2*cc.LO-13)-(12-cc.CZ)
    best_b=[0]*11;best_c=[0]*11;num=[0]*11
    max_b=[0]*11;max_c=[0]*11
    for begin in range(1,MAX_X+1,chunk):
        x=np.arange(begin,min(MAX_X+1,begin+chunk),dtype=np.int64)
        h=np.floor(np.log2(x)).astype(np.int64)
        shift=(21-h)//2
        norm=x<<(2*shift)
        index=norm>>cc.LO
        low=norm&((1<<cc.LO)-1)
        linear=truncated_product(b[index],low,10-cc.BZ)
        square=truncated_product(low,low,13)
        quadratic=truncated_product(c[index],square,12-cc.CZ)
        for s in range(11):
            mask=shift==s
            if not np.any(mask):continue
            num[s]+=int(np.count_nonzero(mask))
            vb=int(np.max(linear[mask]))
            vc=int(np.max(quadratic[mask]))
            max_b[s]=max(vb,max_b[s]);max_c[s]=max(vc,max_c[s])
    for s in range(11):
        best_b[s]=max(0,s+max_b[s].bit_length()-lb)
        best_c[s]=max(0,s+max_c[s].bit_length()-lc)
    pb=max(best_b);pc=max(best_c)
    return dict(passes=pb<=10 and pc<=10,
        checked=MAX_X,width_linear=lb,width_quadratic=lc,
        min_linear_pad=pb,min_quadratic_pad=pc,
        max_bit_heights_by_shift=[dict(shift=s,inputs=num[s],
          linear_bit_length=max_b[s].bit_length(),
          quadratic_bit_length=max_c[s].bit_length(),
          required_linear_pad=best_b[s],required_quadratic_pad=best_c[s])
          for s in range(11)])

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=Path('runs/minimal_phase_padding.json'))
    args=p.parse_args()
    v=analyze()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(v,indent=2)+'\n')
    print(json.dumps(v,indent=2))
    if not v['passes']:raise SystemExit(1)

if __name__=="__main__":main()
