"""Original 48-interval cubic Horner inverse-square-root quantum phase oracle.

This is a qualitatively different algorithm from the v9 quadratic
oracle: normalize x to [1,4), 6-bit j in [16,63], 16-bit residual r.
A third-degree Chebyshev interpolation has only 48 reachable table rows.
Use signed-alternating nonnegative coefficients
  A - B*r/2^16 + C*(r/2^16)^2 - D*(r/2^16)^3
at Q=29, evaluated with integer nested Horner, namely
  t = C - T_16(D,r);
  u = B - T_16(t,r);
  v = T_16(u,r);
  phase = (A-v) * 2^s / 2^29,
where T_d is EXACTLY the sum of partial-product columns i+j>=d, not
floor(ab/2^d). All target arithmetic is unsigned, with t,u>=0.

A is phased and erased independently. B,C,D are jointly fetched but D
is erased early, using measurement-aware cleanup. B,C are restored by
reversing Horner and independently erased. Exact circuit behavior is
checked against the actual T_16 model for all 48*65536 = 3,145,728
normalized mantissas before any full verifier claim.

Independent synthesis: no Kaggle competitor's circuit, ROM, or source.
"""
from __future__ import annotations
import argparse,hashlib,json
from functools import lru_cache
from pathlib import Path
import numpy as np
from .priority1_exact import truncated_product
from .benchmark1_first_row_specialization import mul_first_row
from .benchmark1_binary_block_normalizer import normalize_binary_block,denormalize_binary_block
from .benchmark1_zero_residual_pad import phase_with_pad
from .hybrid_kickmix import resource_counts

P=4;LO=20-P;H=1<<LO;Q=29
CUT=LO
REF=2076880.0
V9_QC=237684.49222389946
TOL=2**-24

@lru_cache(maxsize=1)
def coefficient_table():
    nodes=.5+.5*np.cos(np.pi*(np.arange(4)+.5)/4)
    indices=np.arange(1<<P,4<<P,dtype=np.float64)
    v=np.linalg.solve(np.polynomial.polynomial.polyvander(nodes,3),
       (((indices[:,None]+nodes[None,:])/(1<<P))**(-.5)).T).T
    # For a decreasing, convex inverse square root, signs are (+,-,+,-).
    raw=np.rint(v*(1<<Q)).astype(np.int64)
    out=np.zeros((1<<(P+2),4),dtype=np.int64)
    out[1<<P:]=raw*np.array([1,-1,1,-1],dtype=np.int64)
    if np.any(out<0):raise ValueError("Negative Horner coefficient")
    out.setflags(write=False)
    return out

def full_normalized_screen():
    coeff=coefficient_table()
    r=np.arange(H,dtype=np.int64)
    worst=0.;arg=None;bad=0;min_t=2**40;min_u=2**40
    for j in range(1<<P,4<<P):
        A,B,C,D=map(int,coeff[j])
        inner=C-truncated_product(np.int64(D),r,CUT)
        outer=B-truncated_product(inner,r,CUT)
        phase_num=A-truncated_product(outer,r,CUT)
        min_t=min(min_t,int(inner.min()))
        min_u=min(min_u,int(outer.min()))
        z=(j+r/H)/(1<<P)
        err=np.abs((phase_num/(1<<Q))*np.sqrt(z)-1)/TOL
        bad+=int(np.count_nonzero(err>=1))
        max_ix=int(np.argmax(err))
        if float(err[max_ix])>worst:
            worst=float(err[max_ix]);arg=(j,int(r[max_ix]))
    widths=[max(1,int(v.max()).bit_length()) for v in coeff.T]
    return dict(passes=bad==0 and min_t>=0 and min_u>=0,
       failures=bad,normalized_values=(3<<P)*H,
       worst_error_over_tolerance=worst,worst_j_r=arg,
       minimum_first_horner_result=min_t,
       minimum_second_horner_result=min_u,
       coefficient_widths=widths,coefficient_rows=48,
       physical_table_rows=64,
       partial_product_drops=[CUT]*3,
       coefficient_precision=Q)

def _packed(mat,idx,offsets):
    return [sum(int(mat[j,k])<<offsets[i] for i,k in enumerate(idx))
            for j in range(len(mat))]

def build(*, rom_mode="qrom",free_d_early=True):
    if not full_normalized_screen()["passes"]:
        raise ValueError("Full 3.145M normalized cubic Horner screen failed")
    import kickmix as km
    coeff=coefficient_table()
    wa,wb,wc,wd=[max(1,int(v.max()).bit_length()) for v in coeff.T]
    ta=list(map(int,coeff[:,0]))
    # B,C,D combined low to high, to retain stable native cleanup indices.
    tbcd=_packed(coeff,[1,2,3],[0,wb,wb+wc])
    td=list(map(int,coeff[:,3]))
    tb=list(map(int,coeff[:,1]))
    tc=list(map(int,coeff[:,2]))
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,name="x")
    flags,s=normalize_binary_block(b,x)
    address=x[LO:]
    residual=x[:LO]
    with b.mark("cubic_A_intercept_lookup_phase_erase"):
        a=b.alloc_qubits(wa)
        b.init_lookup(table=ta,address=address,target=a)
        phase_with_pad(b,a,s,Q,+1,10)
        b.del_lookup(table=ta,address=address,target=a)
        b.free(a)
    with b.mark("cubic_Horner_coefficient_load"):
        bank=b.alloc_qubits(wb+wc+wd)
        if rom_mode=="qrom":
            b.init_lookup(table=tbcd,address=address,target=bank)
        elif rom_mode=="anf":
            from .benchmark1_pruned_anf_rom import binary_pruned_rom
            b.reset(bank)
            binary_pruned_rom(b,address,bank,tbcd)
        else:raise ValueError("Invalid ROM mode")
        beta=bank[:wb];gamma=bank[wb:wb+wc];delta=bank[wb+wc:]
    # Compute gamma <- gamma - T_16(delta,r).
    def product(word):
        out=b.alloc_qubits(len(word))
        b.reset(out)
        mul_first_row(b,word,residual,out,CUT)
        return out
    def unproduct(word,out):
        mul_first_row(b,word,residual,out,CUT,inverse=True)
        b.free(out)
    with b.mark("cubic_inner_D_times_r"):
        p=product(delta)
        b.isub(p,target=gamma)
        unproduct(delta,p)
        if free_d_early:
            # Discard the unchanged D coefficient before allocating a
            # 20/24-bit multiply target; release its workspace.
            b.del_lookup(table=td,address=address,target=delta)
            b.free(delta)
    with b.mark("cubic_middle_CminusDr_times_r"):
        p=product(gamma)
        b.isub(p,target=beta)
        unproduct(gamma,p)
    with b.mark("cubic_outer_phase"):
        p=product(beta)
        phase_with_pad(b,p,s,Q,-1,10)
        unproduct(beta,p)
    with b.mark("cubic_reverse_middle"):
        p=product(gamma)
        b.iadd(p,target=beta)
        unproduct(gamma,p)
        b.del_lookup(table=tb,address=address,target=beta)
        b.free(beta)
    with b.mark("cubic_reverse_inner"):
        if free_d_early:
            delta=b.alloc_qubits(wd)
            b.init_lookup(table=td,address=address,target=delta)
        p=product(delta)
        b.iadd(p,target=gamma)
        unproduct(delta,p)
        b.del_lookup(table=td,address=address,target=delta)
        b.free(delta)
        b.del_lookup(table=tc,address=address,target=gamma)
        b.free(gamma)
    denormalize_binary_block(b,x,flags,s)
    return b.finish_circuit()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native',action='store_true')
    ap.add_argument('--sample',type=int,default=128)
    ap.add_argument('--rom-mode',choices=('qrom','anf'),default='qrom')
    ap.add_argument('--retain-d',action='store_true')
    ap.add_argument('--output',type=Path,default=Path('out/cubic48_horner.kmb'))
    args=ap.parse_args()
    report=dict(math=full_normalized_screen())
    if args.native:
        circ=build(rom_mode=args.rom_mode,free_d_early=not args.retain_d)
        stats=resource_counts(circ)
        report.update(native=stats,rom_mode=args.rom_mode,
            reuse_d=not args.retain_d,
            can_beat_v9_at_floor=stats['lower_bound_qc']<V9_QC,
            can_reach_20x_at_floor=stats['lower_bound_qc']<REF/20)
        args.output.parent.mkdir(parents=True,exist_ok=True)
        circ.to_file(args.output,format='kmb')
        report['sha256']=hashlib.sha256(args.output.read_bytes()).hexdigest()
        if args.sample and report['can_beat_v9_at_floor']:
            from .hybrid_verify import boundary_inputs,verify_samples
            report['coherent']=verify_samples(circ,boundary_inputs(random_count=args.sample))
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
