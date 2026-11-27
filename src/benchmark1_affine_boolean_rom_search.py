"""Original Clifford-affine address recoding for a multioutput Boolean ROM.

ANF coefficients c_m of the exact reachable coefficient table may become
sparser when inputs are represented in a different reversible affine
coordinate system y=M*x XOR v. Such input transformations use ONLY
Clifford X and CNOT gates; no Toffoli, no ancilla, no fine Z rotations.

Search 9-variable GL(9,2) elementary CNOT transformations together with
input bit complements, minimizing the actual shared AND-prefix cost of
the multioutput GF(2) ANF decoder. The table's *function* is unchanged:
f_new(y)=f_old(M^-1*(y XOR v)) over all 512 states, including all
384 legal normalized intervals. We independently check this identity
exhaustively before any native compilation or full-verifier promotion.

Unlike variable relabeling, invertible affine mixing can transform
active ANF monomial structure. This experiment reads no competitor
implementation and does not alter the mathematical inverse-sqrt oracle.
"""
from __future__ import annotations
import argparse,json,random,math
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import optimized_coefficients
from .benchmark1_pruned_anf_rom import transform

N=9
MASK=(1<<N)-1

def ops():
    return [('x',i,-1) for i in range(N)]+[
        ('cx',i,j) for i in range(N) for j in range(N) if i!=j]

def apply_coeff(c,op):
    kind,i,j=op
    result=c.copy()
    bit=1<<i
    if kind=='x':
        # Substitute old x_i = new y_i XOR 1.
        for m in range(1<<N):
            if not m&bit:
                result[m]^=c[m|bit]
    elif kind=='cx':
        # Substitute old x_i = new y_i XOR new y_j.
        for m in range(1<<N):
            if m&bit:
                result[(m^bit)|(1<<j)]^=c[m]
    else:raise ValueError("Unknown affine op")
    return result

def prefix_proxy(coeff):
    prefixes=set()
    active=0
    for m,v in enumerate(coeff):
        if not v:continue
        active+=1
        t=m
        while t:
            if t.bit_count()>1:prefixes.add(t)
            t^=1<<(t.bit_length()-1)
    return len(prefixes),active

def transform_bit_value(x,op):
    kind,i,j=op
    if kind=='x':return x^(1<<i)
    return x^(((x>>j)&1)<<i)

def check_affine_identity(original,co,sequence):
    # For every new address y, undo physical CNOT/X operations and
    # compare its transformed function with ANF evaluation.
    new_values=transform(co,N)
    for y in range(1<<N):
        x=y
        for op in reversed(sequence):
            x=transform_bit_value(x,op)
        if new_values[y]!=original[x]:
            raise RuntimeError(f"Original table lost: y={y} x={x}")
    return True

def optimize(seeds=8,iterations=1600):
    src,_=optimized_coefficients('affine')
    mat,_=optimized_coefficients('affine')
    A,B,C=mat[:,0],mat[:,1]>>2,mat[:,2]>>4
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    table=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
           for a,b,c in zip(A,B,C)]
    initial=transform(table,N)
    start_cost,start_active=prefix_proxy(initial)
    best=(start_cost,initial,[])
    rng=random.Random(20261127)
    elementary=ops()
    # First optimize to a strict fixed point.
    current=initial
    greedy_ops=[]
    while True:
        candidates=[(prefix_proxy(apply_coeff(current,op))[0],op)
                    for op in elementary]
        cost,op=min(candidates,key=lambda x:x[0])
        if cost>=prefix_proxy(current)[0]:break
        current=apply_coeff(current,op)
        greedy_ops.append(op)
    if prefix_proxy(current)[0]<best[0]:
        best=(prefix_proxy(current)[0],current,greedy_ops.copy())
    print(json.dumps(dict(stage='greedy',original_prefixes=start_cost,
         original_active_monomials=start_active,
         best_prefixes=best[0],actions=greedy_ops)),flush=True)
    # Annealing explores combinations of changes that may individually
    # increase the prefix proxy while jointly lowering it.
    for seed in range(seeds):
        state=best[1] if seed%2==0 else initial
        path=best[2].copy() if seed%2==0 else []
        cost=prefix_proxy(state)[0]
        for it in range(iterations):
            temp=3.5*(1-it/iterations)+0.05
            op=elementary[rng.randrange(len(elementary))]
            neighbor=apply_coeff(state,op)
            nc=prefix_proxy(neighbor)[0]
            if nc<cost or rng.random()<math.exp(min(0,(cost-nc)/temp)):
                state=neighbor;cost=nc;path.append(op)
                # Two consecutive identical generators cancel exactly.
                if len(path)>=2 and path[-1]==path[-2]:
                    path.pop();path.pop()
                if nc<best[0]:
                    check_affine_identity(table,state,path)
                    best=(nc,state,path.copy())
                    print(json.dumps(dict(stage='anneal',seed=seed,
                         iteration=it,best_prefixes=nc,
                         active_monomials=prefix_proxy(state)[1],
                         action_count=len(path))),flush=True)
    check_affine_identity(table,best[1],best[2])
    out=dict(passes=True,
        original_prefix_ands=start_cost,
        best_prefix_ands=best[0],
        improvement_in_prefix_ands=start_cost-best[0],
        original_active_terms=start_active,
        final_active_terms=prefix_proxy(best[1])[1],
        affine_sequence=[list(op) for op in best[2]],
        native_coherence="NOT_TESTED_YET",
        original_table_rows=512,
        all_table_values_checked=512,
        legal_reachable_rows=384)
    print('SUMMARY '+json.dumps(out),flush=True)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--seeds',type=int,default=8)
    ap.add_argument('--iterations',type=int,default=1600)
    args=ap.parse_args()
    out=optimize(args.seeds,args.iterations)
    path=Path('runs/affine_rom_boolean_search.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
