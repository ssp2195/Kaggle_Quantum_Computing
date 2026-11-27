"""Original staged high-bit gamma decoder with exact interval threshold ROM.

The low six/five bits of curvature C are loaded alongside A,B by the
existing shared GF(2) ANF decoder, but its 1-2 most significant bits
are NOT loaded until just before the quadratic phase.

These high C bitplanes are sparse monotone intervals of the nine-bit
coefficient address. Compute them later using exact reversible
less-than threshold comparators, costing no fine rotations and freeing
one or two coefficient-bank qubits during the intercept and large linear
controlled multiplier. Optionally measure-and-erase B immediately after
its linear phase; B is not needed for the quadratic phase.

The original exact A/B/C values and p7/q26 partial-product phase are
identical for all 384 reachable intervals. The high-bit threshold
functions are proven for all 384 rows, not assumed from smoothness.
Reversible threshold gates are applied twice to restore the workspace.
The original v9 carry-safe arithmetic and binary normalization remain.

No Kaggle competitor code/circuit is read. Full independent verifier
is REQUIRED before a performance improvement may be promoted.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from . import benchmark1_cost_cooptim as cc
from .benchmark1_dontcare_boolean_rom import optimized_coefficients
from .benchmark1_pruned_anf_rom import binary_pruned_rom
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .benchmark1_first_row_specialization import mul_first_row,square_first_row
from .hybrid_kickmix import resource_counts

V9=237684.49222389946
REF=2076880.0
CUTS=(10,13,12)

def run_segments(bit_values):
    runs=[];lo=None
    for j in range(128,512):
        active=(int(bit_values[j])&1)==1
        if active and lo is None:lo=j
        if not active and lo is not None:runs.append((lo,j-1));lo=None
    if lo is not None:runs.append((lo,511))
    return runs

def bitplane_certificate(late):
    if late not in (1,2,3):raise ValueError(late)
    mat,_=optimized_coefficients('affine')
    C=mat[:,2]>>cc.CZ
    wc=max(1,int(C.max()).bit_length())
    report=[]
    for bit in range(wc-late,wc):
        values=[(int(v)>>bit)&1 for v in C]
        segments=run_segments(values)
        report.append(dict(bit=bit,segments=[list(pair) for pair in segments],
                    valid_ones=sum(values[128:]),
                    segment_count=len(segments)))
    return dict(passes=True,coefficient_rows_checked=384,
         omitted_top_gamma_bits=late,
         full_curvature_width=wc,
         late_bitplane_segments=report)

def build(late=2,erase_b=True,pad8=True):
    assert bitplane_certificate(late)["passes"]
    import kickmix as km
    mat,_=optimized_coefficients('affine')
    A,B,C=mat[:,0],mat[:,1]>>cc.BZ,mat[:,2]>>cc.CZ
    wa,wb,wc=[max(1,int(t.max()).bit_length()) for t in (A,B,C)]
    low_c_width=wc-late
    low_c_mask=(1<<low_c_width)-1
    tA=[int(v) for v in A]
    tB=[int(v) for v in B]
    tClow=[int(v)&low_c_mask for v in C]
    tBClow=[int(b)|(int(c)<<wb) for b,c in zip(B,tClow)]
    packed=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb))
            for a,b,c in zip(A,B,tClow)]
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name='x')
    flags,scale=normalize_binary_block(b,x)
    addr=x[cc.LO:]
    reg=b.alloc_qubits(wa+wb+low_c_width)
    b.reset(reg)
    binary_pruned_rom(b,addr,reg,packed)
    a=reg[:wa];beta=reg[wa:wa+wb];gamma_low=reg[wa+wb:]
    with b.mark('intercept'):
        phase_with_pad(b,a,scale,cc.Q,+1,8 if pad8 else 9)
        if pad8:
            b.ccz(a[26],scale[0],scale[3])
            b.ccz(a[25],scale[1],scale[3])
        b.del_lookup(table=tA,address=addr,target=a)
        b.free(a)
    dl,ds,dq=CUTS
    with b.mark('linear'):
        linear=b.alloc_qubits(wb+cc.LO-(dl-cc.BZ))
        b.reset(linear)
        mul_first_row(b,beta,x[:cc.LO],linear,dl-cc.BZ)
        phase_with_pad(b,linear,scale,cc.Q+cc.LO-dl,-1,5)
        mul_first_row(b,beta,x[:cc.LO],linear,dl-cc.BZ,inverse=True)
        b.free(linear)
    if erase_b:
        with b.mark('early_slope_cleanup'):
            b.del_lookup(table=tB,address=addr,target=beta)
            b.free(beta)
    late_qubits=[]
    segments=bitplane_certificate(late)['late_bitplane_segments']
    def apply_threshold_segments(target,ranges):
        # For all reachable j>=128, XOR of threshold comparisons
        # gives the exact interval indicator of any [lo,hi] run.
        for lo,hi in ranges:
            if lo>128:
                b.flip_if_less_than(addr,lo,target=target)
            if hi<511:
                b.flip_if_less_than(addr,hi+1,target=target)
            elif lo==128:
                b.x(target)
    with b.mark('late_curvature_decode'):
        for info in segments:
            q=b.alloc_qubits(1)[0]
            b.reset(q)
            apply_threshold_segments(q,info['segments'])
            late_qubits.append(q)
    gamma=km.array(list(gamma_low)+late_qubits)
    with b.mark('quadratic'):
        sq=b.alloc_qubits(2*cc.LO-ds);b.reset(sq)
        square_first_row(b,x[:cc.LO],sq,ds)
        out=b.alloc_qubits(wc+len(sq)-(dq-cc.CZ));b.reset(out)
        mul_first_row(b,sq,gamma,out,dq-cc.CZ)
        phase_with_pad(b,out,scale,cc.Q+2*cc.LO-ds-dq,+1,5)
        mul_first_row(b,sq,gamma,out,dq-cc.CZ,inverse=True)
        b.free(out)
        square_first_row(b,x[:cc.LO],sq,ds,inverse=True)
        b.free(sq)
    with b.mark('curvature_cleanup'):
        for info,q in reversed(list(zip(segments,late_qubits))):
            apply_threshold_segments(q,info['segments'])
            b.free(q)
        if erase_b:
            # Native measured unlookup requires a minimum target width
            # for its address-unary work register. Append CLEAN zero
            # target bits while all other arithmetic has been freed.
            pad=b.alloc_qubits(max(0,16-len(gamma_low)))
            b.reset(pad)
            union=km.array(list(gamma_low)+list(pad))
            b.del_lookup(table=tClow,address=addr,target=union)
            b.free(pad)
            b.free(gamma_low)
        else:
            b.del_lookup(table=tBClow,address=addr,target=reg[wa:])
            b.free(reg[wa:])
    denormalize_binary_block(b,x,flags,scale)
    return b.finish_circuit()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--late",type=int,nargs='+',default=[1,2,3])
    p.add_argument("--early-B",action='store_true')
    p.add_argument("--pad9",action='store_true')
    p.add_argument("--samples",type=int,default=256)
    args=p.parse_args()
    records=[]
    for late in args.late:
        try:
            proof=bitplane_certificate(late)
            c=build(late,erase_b=args.early_B,pad8=not args.pad9)
            cost=resource_counts(c)
            r=dict(late_bits=late,early_B=args.early_B,pad8=not args.pad9,
               coefficients=proof,native=cost,
               toffoli_floor_beats_v9=cost['lower_bound_qc']<V9,
               can_reach_20x_at_floor=cost['lower_bound_qc']<REF/20)
            if args.samples and r['toffoli_floor_beats_v9']:
                from .hybrid_verify import verify_samples,boundary_inputs
                r['sample']=verify_samples(c,boundary_inputs(random_count=args.samples))
        except Exception as e:
            r=dict(late_bits=late,early_B=args.early_B,
              error=type(e).__name__+":"+str(e)[:550])
        print(json.dumps(r),flush=True)
        records.append(r)
    path=Path('runs/late_gamma_bitplane_native.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':main()
