"""Original folded-normalized-range direct inverse-sqrt phase oracle.

Normalize x by powers of four as the verified v6 circuit does, then
fold its mantissa z in [1,4) into w in [1,2):
   if z<2:  w=z,  phase=f(w)
   if z>=2: w=z/2, phase=f(w)/sqrt(2).
A coherent parity flag selects one of two independently fitted p7
quadratic coefficient sets (128 intervals each). The packed QROM has
256 entries/address8, versus the v6 512/address9. Unfold restores
all original input wires exactly.

Coefficients are re-fitted specifically for exact native omitted-column
partial products, B divisible by 4, C by 16, and direct phased
intercept/linear/quadratic. All 2*128*16384=4,194,304 normalized
mantissa/parity pairs are checked BEFORE native compilation.

This source and coefficients are authored independently; no competitor
submission, lookup, KMB or ROM table is used.
"""
from __future__ import annotations
import argparse,hashlib,json
from functools import lru_cache
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .hybrid_kickmix import resource_counts
from .benchmark1_binary_block_normalizer import (
    normalize_binary_block,denormalize_binary_block)
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled
from .priority1_truncated import _mul_cut

P=7;LO=14;H=1<<LO;Q=26
BZ=2;CZ=4
CUTS=(10,13,10)
TOL=2**-24
MAX_X=3133452
REFERENCE=2076880.0
INCUMBENT=268172.78036771005
MESH=np.unique(np.r_[np.arange(0,H,137),np.arange(0,256),
      np.arange(H-256,H),np.linspace(0,H-1,200,dtype=int)]).astype(np.int64)
ALLR=np.arange(H,dtype=np.int64)

def _fit_one(par,jj,seedB,seedC,r,square):
    dl,ds,dq=CUTS
    z=(jj+r/H)*(2**par)/(1<<P)
    want=(1<<Q)/np.sqrt(z)
    weight=np.sqrt(z)/((1<<Q)*TOL)
    B0=round(float(seedB)/(1<<BZ))
    C0=round(float(seedC)/(1<<CZ))
    best=(float("inf"),None)
    for dc in range(-8,9):
        C=C0+dc
        if C<0:continue
        quad=truncated_product(np.int64(C),square,dq-CZ)/2**(2*LO-ds-dq)
        for db in range(-8,9):
            B=B0+db
            if B<0:continue
            lin=truncated_product(np.int64(B),r,dl-BZ)/2**(LO-dl)
            wanted=want+lin-quad
            mid=(float(wanted.min())+float(wanted.max()))/2
            for A in (int(np.floor(mid)),int(np.ceil(mid))):
                loss=float(np.max(np.abs((A-wanted)*weight)))
                if loss<best[0]:best=(loss,(A,B<<BZ,C<<CZ))
    if best[1] is None:raise ValueError("Unable to fit folded interval "+str((par,jj)))
    return best

@lru_cache(maxsize=1)
def coefficients():
    nodes=.5+.5*np.cos(np.pi*(np.arange(3)+.5)/3)
    vand=np.polynomial.polynomial.polyvander(nodes,2)
    jj=np.arange(128,256,dtype=np.float64)
    squares=truncated_product(MESH,MESH,CUTS[1])
    out=np.zeros((256,3),dtype=np.int64)
    for par in (0,1):
        y=(2**(-par/2))/np.sqrt((jj[:,None]+nodes[None,:])/(1<<P))
        v=np.linalg.solve(vand,y.T).T
        b=np.rint((-v[:,1]*(1<<Q))/(1<<BZ)).astype(np.int64)*(1<<BZ)
        c=np.rint((v[:,2]*(1<<Q))/(1<<CZ)).astype(np.int64)*(1<<CZ)
        for idx,j in enumerate(range(128,256)):
            score,pair=_fit_one(par,j,int(b[idx]),int(c[idx]),MESH,squares)
            out[par*128+idx]=pair
    out.setflags(write=False)
    return out

def full_normalized_screen():
    table=coefficients()
    r=ALLR
    square=truncated_product(r,r,CUTS[1])
    worst=0.;bad=0;arg=None
    for par in (0,1):
        for j in range(128,256):
            A,B,C=(int(v) for v in table[par*128+j-128])
            lin=truncated_product(np.int64(B>>BZ),r,CUTS[0]-BZ)
            quad=truncated_product(np.int64(C>>CZ),square,CUTS[2]-CZ)
            z=(j+r/H)*(2**par)/(1<<P)
            phase=A/(1<<Q)-lin/(1<<(Q+LO-CUTS[0]))+quad/(1<<(Q+2*LO-CUTS[1]-CUTS[2]))
            err=np.abs(phase*np.sqrt(z)-1)/TOL
            bad+=int(np.count_nonzero(err>=1))
            k=int(np.argmax(err))
            if float(err[k])>worst:
                worst=float(err[k]);arg=(par,j,int(r[k]))
    return dict(passes=bad==0,failures=bad,
        normalized_folded_pairs=2*128*H,
        worst_error_ratio=worst,worst_parity_j_r=arg,
        segments_per_parity=128,physical_rom_entries=256,
        cuts=CUTS,B_zero_bits=BZ,C_zero_bits=CZ)

def pad_certificate(chunk=180000):
    table=coefficients()
    B=table[:,1]>>BZ
    C=table[:,2]>>CZ
    wb=max(1,int(B.max()).bit_length())
    wc=max(1,int(C.max()).bit_length())
    wl=wb+LO-(CUTS[0]-BZ)
    wq=wc+(2*LO-CUTS[1])-(CUTS[2]-CZ)
    maxL=[0]*11;maxQ=[0]*11
    for first in range(1,MAX_X+1,chunk):
        x=np.arange(first,min(MAX_X+1,first+chunk),dtype=np.int64)
        k=np.floor(np.log2(x)).astype(np.int64)//2
        s=10-k
        z=x<<(2*s)
        par=(z>>21)&1
        folded=np.where(par==0,z<<1,z)
        j=(folded>>LO)-128
        if np.any((j<0)|(j>=128)):raise RuntimeError("Folded address not in [128,256)")
        addr=j+par*128
        r=folded&(H-1)
        lin=truncated_product(B[addr],r,CUTS[0]-BZ)
        square=truncated_product(r,r,CUTS[1])
        quad=truncated_product(C[addr],square,CUTS[2]-CZ)
        for sh in range(11):
            mask=s==sh
            if mask.any():
                maxL[sh]=max(maxL[sh],int(lin[mask].max()))
                maxQ[sh]=max(maxQ[sh],int(quad[mask].max()))
    pl=max(0,max(i+x.bit_length()-wl for i,x in enumerate(maxL)))
    pq=max(0,max(i+x.bit_length()-wq for i,x in enumerate(maxQ)))
    return dict(passes=pl<=10 and pq<=10,
        legal_inputs=MAX_X,linear_word_width=wl,quadratic_word_width=wq,
        required_linear_pad=pl,required_quadratic_pad=pq,
        max_linear_by_shift=maxL,max_quadratic_by_shift=maxQ)

def build():
    import kickmix as km
    model=full_normalized_screen()
    if not model["passes"]:raise ValueError("Folded quadratic accuracy error: "+str(model))
    pads=pad_certificate()
    if not pads["passes"]:raise ValueError("Folded phase physical pad error: "+str(pads))
    mat=coefficients()
    A,B,C=mat[:,0],mat[:,1]>>BZ,mat[:,2]>>CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    tA=[int(v) for v in A]
    tBC=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    combined=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
              for a,b,c in zip(A,B,C)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,scale=normalize_binary_block(b,x)
    parity=b.alloc_qubits(1)[0]
    b.reset(parity)
    b.cx(x[21],parity)
    # w=z for parity=0; w=z/2 for parity=1. The high bit of
    # z is zero exactly when a one-bit left rotation is needed.
    b.x(parity)
    b.left_rotate(x,1,control=parity)
    b.x(parity)
    address=km.array([x[i] for i in range(LO,21)]+[parity])
    if len(address)!=8:raise RuntimeError("Unexpected folded address size")
    coeff=b.alloc_qubits(wa+wb+wc)
    b.init_lookup(table=combined,address=address,target=coeff)
    a,beta,gamma=coeff[:wa],coeff[wa:wa+wb],coeff[wa+wb:]
    with b.mark("folded_constant"):
        phase_with_pad(b,a,scale,Q,+1,10)
        b.del_lookup(table=tA,address=address,target=a)
        b.free(a)
    dl,ds,dq=CUTS
    with b.mark("folded_linear"):
        linear=b.alloc_qubits(wb+LO-(dl-BZ));b.reset(linear)
        mul_controlled(b,beta,x[:LO],linear,dl-BZ)
        phase_with_pad(b,linear,scale,Q+LO-dl,-1,pads["required_linear_pad"])
        mul_controlled(b,beta,x[:LO],linear,dl-BZ,inverse=True)
        b.free(linear)
    with b.mark("folded_quadratic"):
        sq=b.alloc_qubits(2*LO-ds);b.reset(sq)
        _mul_cut(b,x[:LO],x[:LO],sq,ds)
        prod=b.alloc_qubits(wc+len(sq)-(dq-CZ));b.reset(prod)
        mul_controlled(b,sq,gamma,prod,dq-CZ)
        phase_with_pad(b,prod,scale,Q+2*LO-ds-dq,+1,pads["required_quadratic_pad"])
        mul_controlled(b,sq,gamma,prod,dq-CZ,inverse=True)
        b.free(prod)
        _mul_cut(b,x[:LO],x[:LO],sq,ds,inverse=True)
        b.free(sq)
    b.del_lookup(table=tBC,address=address,target=coeff[wa:])
    b.free(coeff[wa:])
    b.x(parity)
    b.right_rotate(x,1,control=parity)
    b.x(parity)
    b.cx(x[21],parity)
    b.free(parity)
    denormalize_binary_block(b,x,flags,scale)
    return b.finish_circuit()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native",action="store_true")
    ap.add_argument("--sample",type=int,default=128)
    ap.add_argument("--output",type=Path,default=Path("out/original_folded_p7.kmb"))
    args=ap.parse_args()
    result=dict(math=full_normalized_screen(),padding=pad_certificate())
    if args.native and result["math"]["passes"] and result["padding"]["passes"]:
        c=build()
        costs=resource_counts(c)
        result["native"]=costs
        result["toffoli_floor_can_beat_v6"]=costs["lower_bound_qc"]<INCUMBENT
        result["toffoli_floor_can_reach_20x"]=costs["lower_bound_qc"]<REFERENCE/20
        args.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(args.output,format="kmb")
        result["sha256"]=hashlib.sha256(args.output.read_bytes()).hexdigest()
        if args.sample and costs["lower_bound_qc"]<INCUMBENT:
            from .hybrid_verify import verify_samples,boundary_inputs
            result["coherent"]=verify_samples(c,boundary_inputs(random_count=args.sample))
            if result["coherent"]["failures"]:
                raise RuntimeError("Quantum coherent folded phase failed: "+str(result))
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
