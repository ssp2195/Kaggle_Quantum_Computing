"""Independent native v5 stage cost decomposition, per-operation and live wires.

Diagnostic only. Writes original circuit to .kmx in the CI runner and scans
per-stage INFO_TAGs to locate Toffoli and rotation bottlenecks.
"""
from collections import Counter,defaultdict
from pathlib import Path
import json,re
from .benchmark1_best_verified_v5 import build_verified
from .hybrid_kickmix import resource_counts

def profile(out=Path("runs/v5_stage_profile.json")):
    c=build_verified()
    tmp=out.parent/"profile-v5.kmx"
    tmp.parent.mkdir(parents=True,exist_ok=True)
    c.to_file(tmp,format="kmx")
    stage=defaultdict(Counter); stack=[]; notes=[]
    for line in tmp.read_text().splitlines():
        token=line.split(" ",1)[0] if line.strip() else ""
        if token=="INFO_TAG":
            notes.append(line)
            m=re.search(r'@(?:enter|exit):[^"]+',line)
            if m:
                label=m.group(0)
                if label.startswith("@enter:"):
                    stack.append(label[len("@enter:"):])
                elif label.startswith("@exit:"):
                    wanted=label[len("@exit:"):]
                    if stack and stack[-1]==wanted:
                        stack.pop()
                    else:
                        raise ValueError("Unbalanced native tag: "+label)
            continue
        # Attribute to the top-level algorithmic mark, rather than internal
        # kickmix lookup, adder or rotation helper marks.
        roots=("constant_A_phase","linear_phase","quadratic_phase")
        key=next((x for x in stack if x in roots),"unmarked")
        stage[key][token]+=1
    result=dict(resources=resource_counts(c),info_tags=notes[:25],
                num_info_tags=len(notes),
                opcode_counts=dict(sum((Counter({k:v}) for z in stage.values() for k,v in z.items()),Counter())),
                stage_opcode_counts={k:dict(v) for k,v in stage.items()})
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2)[:16000])
    tmp.unlink()
if __name__=="__main__":profile()
