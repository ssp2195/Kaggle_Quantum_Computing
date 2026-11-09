"""Native centered signed-digit Horner comparisons and full-domain screens."""
import argparse
import json
from pathlib import Path
from .benchmark1_centered_horner import CenteredConfig,screen,build
from .hybrid_kickmix import resource_counts
from .hybrid_verify import verify_samples,boundary_inputs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=["screen","build","verify"])
    p.add_argument("--segment-bits",type=int,default=6)
    p.add_argument("--guard",type=int,default=4)
    p.add_argument("--lookup",choices=["qrom","xag"],default="qrom")
    p.add_argument("--cleanup",choices=["mbuc","coherent"],default="mbuc")
    p.add_argument("--output",type=Path)
    p.add_argument("--circuit",type=Path)
    p.add_argument("--random-count",type=int,default=8192)
    a=p.parse_args()
    cfg=CenteredConfig(a.segment_bits,a.guard,a.lookup,a.cleanup)
    if a.command=="screen":
        result=screen(cfg)
    elif a.command=="build":
        if a.output is None:p.error("--output required")
        c=build(cfg)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        c.to_file(a.output,format="kmb")
        result=dict(config=vars(cfg),resources=resource_counts(c),
                    mathematical=screen(cfg),file=str(a.output))
    else:
        if a.circuit is None:p.error("--circuit required")
        import kickmix as km
        c=km.Circuit.from_file(a.circuit,format="kmb")
        result=verify_samples(c,boundary_inputs(random_count=a.random_count))
    print(json.dumps(result,indent=2,sort_keys=True))
    if a.output and a.command!="build":
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+"\n")
    if result.get("failures",0):raise SystemExit(1)


if __name__=="__main__":
    main()
