"""Original sparsity-pruned multioutput Boolean coefficient ROM.

A Möbius transform produces a shared multi-output ANF. Skip AND subtrees
whose entire output polynomial is zero, and emit only the necessary shared
prefix ANDs. The coherent cleanup is the same XOR map in reverse, rather
than a measured unlookup. Search hybrid XAG/QROM arrangements at native
cost, retaining the independently verified v6 polynomial and normalizer.

No competitor source, lookup synthesis, or circuit is imported.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_controlled_cut import mul_controlled
from .priority1_truncated import _mul_cut
from .hybrid_kickmix import resource_counts
from . import benchmark1_cost_cooptim as cc

BASE=268172.78036771005;GOAL=2076880/20
CUTS=(10,13,12)

def transform(table, n):
    coeff=list(map(int,table))
    if len(coeff)!=(1<<n):raise ValueError("Invalid ANF table shape")
    for i in range(n):
        for mask in range(1<<n):
            if mask&(1<<i):coeff[mask]^=coeff[mask^(1<<i)]
    return coeff

def _descendant_mask(n,mask,start):
    remaining=((1<<n)-1)^((1<<start)-1)
    remaining&=~mask
    ans=0
    sub=remaining
    while True:
        ans|=1<<(mask|sub)
        if sub==0:break
        sub=(sub-1)&remaining
    return ans

def binary_pruned_rom(builder,address,target,table):
    n=len(address)
    coeff=transform(table,n)
    if any(v.bit_length()>len(target) for v in coeff):raise ValueError("Target too small")
    active=[i for i,v in enumerate(coeff) if v]
    # Mark all masks needed as a prefix of at least one nonzero monomial
    present=set()
    for mask in active:
        subset=mask
        while True:
            present.add(subset)
            if subset==0:break
            highest=subset.bit_length()-1
            subset^=1<<highest
    def xor_mask(ctl,w):
        i=0
        while w:
            if w&1:
                if ctl is None:builder.x(target[i])
                else:builder.cx(ctl,target[i])
            i+=1;w>>=1
    xor_mask(None,coeff[0])
    counter=[0]
    def visit(start,parent,mask):
        for j in range(start,n):
            child=mask|(1<<j)
            if child not in present:continue
            if parent is None:
                v=address[j];scratch=False
            else:
                v=builder.alloc_qubits(1)[0]
                builder.init_and(parent,address[j],v)
                counter[0]+=1;scratch=True
            xor_mask(v,coeff[child])
            visit(j+1,v,child)
            if scratch:
                builder.del_and(parent,address[j],v)
                builder.free(v)
    visit(0,None,0)
    return dict(used_anf_terms=len(active),used_ancilla_ands=counter[0],
      active_prefix_terms=len(present),total_monomials=1<<n)

def build(*,init="qrom",erase_a="qrom",erase_bc="qrom"):
    if not all(x in ('qrom','anf') for x in (init,erase_a,erase_bc)):
        raise ValueError("mode unknown")
    if not cc.full_screen(CUTS)["passes"]:raise ValueError("Math failed")
    import kickmix as km
    mat=cc.coefficients(CUTS)
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    tA=[int(v) for v in A]
    tBC=[int(b)|(int(c)<<wb) for b,c in zip(B,C)]
    tALL=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
          for a,b,c in zip(A,B,C)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,shift=normalize_binary_block(b,x)
    reg=b.alloc_qubits(wa+wb+wc)
    addr=x[cc.LO:]
    def load(tgt,table,mode,unload=False):
        if mode=="qrom":
            if unload:b.del_lookup(table=table,address=addr,target=tgt)
            else:b.init_lookup(table=table,address=addr,target=tgt)
        else:
            binary_pruned_rom(b,addr,tgt,table)
    load(reg,tALL,init)
    a=reg[:wa];beta=reg[wa:wa+wb];gamma=reg[wa+wb:]
    phase_with_pad(b,a,shift,cc.Q,+1,10)
    load(a,tA,erase_a,True)
    b.free(a)
    dl,ds,dq=CUTS
    product=b.alloc_qubits(wb+cc.LO-(dl-cc.BZ));b.reset(product)
    mul_controlled(b,beta,x[:cc.LO],product,dl-cc.BZ)
    phase_with_pad(b,product,shift,cc.Q+cc.LO-dl,-1,5)
    mul_controlled(b,beta,x[:cc.LO],product,dl-cc.BZ,inverse=True)
    b.free(product)
    square=b.alloc_qubits(2*cc.LO-ds);b.reset(square)
    _mul_cut(b,x[:cc.LO],x[:cc.LO],square,ds)
    product=b.alloc_qubits(wc+len(square)-(dq-cc.CZ));b.reset(product)
    mul_controlled(b,square,gamma,product,dq-cc.CZ)
    phase_with_pad(b,product,shift,cc.Q+2*cc.LO-ds-dq,+1,5)
    mul_controlled(b,square,gamma,product,dq-cc.CZ,inverse=True)
    b.free(product)
    _mul_cut(b,x[:cc.LO],x[:cc.LO],square,ds,inverse=True)
    b.free(square)
    load(reg[wa:],tBC,erase_bc,True)
    b.free(reg[wa:])
    denormalize_binary_block(b,x,flags,shift)
    return b.finish_circuit()

def main():
    rows=[]
    for init,era,erbc in [
        ("qrom","qrom","qrom"),
        ("qrom","anf","qrom"),
        ("qrom","qrom","anf"),
        ("qrom","anf","anf"),
        ("anf","qrom","qrom"),
        ("anf","anf","anf")
    ]:
        try:
            c=build(init=init,erase_a=era,erase_bc=erbc)
            costs=resource_counts(c)
            entry=dict(init=init,erase_a=era,erase_bc=erbc,
              native=costs,lower_bound_below_v6=costs["lower_bound_qc"]<BASE,
              lower_bound_can_reach_20x=costs["lower_bound_qc"]<GOAL)
            if entry["lower_bound_below_v6"] and (init,era,erbc)!=("qrom","qrom","qrom"):
                from .hybrid_verify import boundary_inputs,verify_samples
                entry["coherent"]=verify_samples(c,boundary_inputs(random_count=128))
        except Exception as e:
            entry=dict(init=init,erase_a=era,erase_bc=erbc,
                       error=type(e).__name__+":"+str(e)[:350])
        rows.append(entry);print(json.dumps(entry),flush=True)
    path=Path('runs/pruned_anf_rom_sweep.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rows,indent=2)+'\n')
if __name__=="__main__":main()
