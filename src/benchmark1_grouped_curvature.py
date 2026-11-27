"""Independent grouped-curvature coefficient quantization.

The piecewise quadratic coefficient C(j) is shared between pairs of
adjacent normalized intervals. Their individual A(j), B(j) are jointly
refitted against the exact binary carry/omitted-product approximation.
If a full 3*128*8192 normalized-domain screen passes, compile with the
independently full-verified v6 binary normalizer and phase-padding logic.

This is a NEW approximation/lookup transformation; never claim a
quantum score before native and full independent verifier pass.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from functools import lru_cache
import numpy as np
from . import benchmark1_cost_cooptim as cc
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import build as phase_builder

CUTS=(10,13,12)
H=8192
ALLR=np.arange(H,dtype=np.int64)
MESH=np.unique(np.r_[np.arange(0,H,41),np.arange(H-150,H),
          np.linspace(0,H-1,160,dtype=int)]).astype(np.int64)
REFERENCE=2076880
INCUMBENT=268172.78036771005

def per_row(j,B,C,r,square):
    z=(j+r/H)/(1<<cc.P)
    target=(1<<cc.Q)/np.sqrt(z)
    l=truncated_product(B>>cc.BZ,r,CUTS[0]-cc.BZ)
    q=truncated_product(C>>cc.CZ,square,CUTS[2]-cc.CZ)
    wanted=target+l/(1<<(cc.LO-CUTS[0]))-q/(1<<(2*cc.LO-CUTS[1]-CUTS[2]))
    a=int(np.floor((wanted.max()+wanted.min())/2))
    weight=np.sqrt(z)/((1<<cc.Q)*cc.TOL)
    best=(float('inf'),None)
    for A in (a-1,a,a+1,a+2):
        score=float(np.max(np.abs((A-wanted)*weight)))
        if score<best[0]:best=(score,(A,B,C))
    return best

def per_pair(j,org,r,square,*,cspan=5,bspan=20):
    startC=round(float(org[j,2]+org[j+1,2])/(2*(1<<cc.CZ)))
    best=(float('inf'),None)
    for dc in range(-cspan,cspan+1):
        C=(startC+dc)*(1<<cc.CZ)
        if C<0:continue
        local=[]
        for k in (j,j+1):
            opt=(float('inf'),None)
            bcenter=int(org[k,1])>>cc.BZ
            for dB in range(-bspan,bspan+1):
                B=(bcenter+dB)*(1<<cc.BZ)
                if B<0:continue
                score=per_row(k,B,C,r,square)
                if score[0]<opt[0]:opt=score
            local.append(opt)
        loss=max(local[0][0],local[1][0])
        if loss<best[0]:best=(loss,local)
    return best

@lru_cache(maxsize=1)
def coefficients():
    org=cc.coefficients(CUTS)
    mesh_square=truncated_product(MESH,MESH,CUTS[1])
    dense_square=truncated_product(ALLR,ALLR,CUTS[1])
    data=np.array(org,copy=True)
    for j in range(1<<cc.P,4<<cc.P,2):
        sample=per_pair(j,org,MESH,mesh_square)
        if sample[1] is None:raise ValueError(f"No candidate for j={j}")
        data[j],data[j+1]=sample[1][0][1],sample[1][1][1]
        exact=max(
           per_row(j,int(data[j,1]),int(data[j,2]),ALLR,dense_square)[0],
           per_row(j+1,int(data[j+1,1]),int(data[j+1,2]),ALLR,dense_square)[0])
        if exact>=1:
            improved=per_pair(j,org,ALLR,dense_square,cspan=7,bspan=25)
            if improved[1] is not None:
                data[j],data[j+1]=improved[1][0][1],improved[1][1][1]
    data.setflags(write=False)
    return data

def full_screen():
    mat=coefficients()
    r=ALLR
    square=truncated_product(r,r,CUTS[1])
    worst=0;bad=0;at_j=None
    for j in range(1<<cc.P,4<<cc.P):
        A,B,C=map(int,mat[j])
        z=(j+r/H)/(1<<cc.P)
        l=truncated_product(B>>cc.BZ,r,CUTS[0]-cc.BZ)
        q=truncated_product(C>>cc.CZ,square,CUTS[2]-cc.CZ)
        err=np.abs((A/(1<<cc.Q)-l/(1<<(cc.Q+cc.LO-CUTS[0]))+
                     q/(1<<(cc.Q+2*cc.LO-CUTS[1]-CUTS[2])))*np.sqrt(z)-1)/cc.TOL
        nb=int(np.count_nonzero(err>=1))
        bad+=nb
        if float(err.max())>worst:worst=float(err.max());at_j=j
    return dict(passes=bad==0,failures=bad,
         max_relative_error_over_tolerance=worst,worst_j=at_j,
         normalized_words=3*(1<<cc.P)*H,
         curvature_pair_groups=(3*(1<<cc.P))//2,
         coefficient_shape=list(mat.shape),
         cuts=CUTS)

def build():
    report=full_screen()
    if not report["passes"]:raise ValueError("Grouped curvature failed: "+str(report))
    oldc,oldn,oldu=cc.coefficients,cc.normalize_compact,cc.denormalize_compact
    cc.coefficients=lambda cuts=CUTS:coefficients()
    cc.normalize_compact=normalize_binary_block
    cc.denormalize_compact=denormalize_binary_block
    try:
        circuit,calls=phase_builder()
    finally:
        cc.coefficients=oldc
        cc.normalize_compact=oldn
        cc.denormalize_compact=oldu
    return circuit

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--native',action='store_true')
    p.add_argument('--output',type=Path,default=Path('out/grouped_curvature.kmb'))
    a=p.parse_args()
    report=dict(math=full_screen())
    if a.native and report["math"]["passes"]:
        c=build()
        stats=resource_counts(c)
        report["native"]=stats
        report["could_reach_20x"]=stats["lower_bound_qc"]<REFERENCE/20
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        if stats["lower_bound_qc"]<INCUMBENT:
            from .hybrid_verify import verify_samples,boundary_inputs
            report["coherent"]=verify_samples(c,boundary_inputs(random_count=128))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
