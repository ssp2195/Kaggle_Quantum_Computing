"""Independent measurement-aware coefficient-bit plane reuse feasibility.

Find coefficient output bits that can be synthesized cheaply *after*
the main multioutput ANF decode, allowing their original zero-valued
register qubits to act as clean Boolean-prefix scratch during decoding.
This may lower native peak qubits without recomputing the 374-entry
nonlinear decoder.

For each valid j=128..511, count 1-runs, address parity complexity,
and individual bit ANF prefix complexity. Unreachable j<128 are
unconstrained and may be chosen for minimum additional implementation
cost. All reconstructed valid values must be exact.
"""
from __future__ import annotations
import json
from pathlib import Path
from .benchmark1_dontcare_boolean_rom import optimized_coefficients
from .benchmark1_pruned_anf_rom import transform

def main():
    mat,_=optimized_coefficients('affine')
    A,B,C=mat[:,0],mat[:,1]>>2,mat[:,2]>>4
    widths=[max(1,int(v.max()).bit_length()) for v in (A,B,C)]
    wa,wb,wc=widths
    rows=[int(a)|(int(b)<<wa)|(int(c)<<(wa+wb)) for a,b,c in zip(A,B,C)]
    physical=rows[128:]
    reports=[]
    for bit in range(sum(widths)):
        values=[(row>>bit)&1 for row in physical]
        transitions=sum(a!=b for a,b in zip(values,values[1:]))
        one_runs=sum(v==1 and (i==0 or values[i-1]==0)
                     for i,v in enumerate(values))
        # Complexity of individual output under the optimized full-table ANF.
        bit_tab=[(row>>bit)&1 for row in rows]
        c=transform(bit_tab,9)
        active=[i for i,v in enumerate(c) if v]
        pref=set()
        for m in active:
            while m:
                if m.bit_count()>1:pref.add(m)
                m^=1<<(m.bit_length()-1)
        # Find the best parity-affine hypothesis over reachable rows.
        parity_best=(385,None,None)
        for mask in range(512):
            base=[((j&mask).bit_count()&1) for j in range(128,512)]
            errors=sum(a!=b for a,b in zip(values,base))
            n=min(errors,384-errors)
            if n<parity_best[0]:
                parity_best=(n,mask,1 if errors>192 else 0)
        category=('A' if bit<wa else ('B' if bit<wa+wb else 'C'))
        entry=dict(output_bit=bit,word=category,
             local_bit=(bit if category=='A' else (bit-wa if category=='B' else bit-wa-wb)),
             valid_one_runs=one_runs,valid_transitions=transitions,
             individual_anf_prefixes=len(pref),
             individual_anf_terms=len(active),
             best_parity_errors=parity_best[0],
             best_parity_mask=parity_best[1],
             best_parity_offset=parity_best[2],
             minimum_valid_ones=sum(values),maximum_valid_ones=sum(values))
        reports.append(entry)
    by_runs=sorted(reports,key=lambda x:(x['valid_transitions'],x['individual_anf_prefixes']))
    by_anf=sorted(reports,key=lambda x:(x['individual_anf_prefixes'],x['valid_transitions']))
    result=dict(widths=widths,total_coeff_bits=sum(widths),
                best_by_transitions=by_runs[:20],
                best_by_anf_prefixes=by_anf[:20],
                all_bitplanes=reports)
    path=Path('runs/coeff_bitplane_stats.json')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='all_bitplanes'},indent=2))

if __name__=='__main__':main()
