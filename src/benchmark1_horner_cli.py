"""Generate and evaluate Horner architecture with native kickmix."""
import argparse
import json
from pathlib import Path
from .benchmark1_horner import HornerConfig,screen,build
from .hybrid_kickmix import resource_counts
from .hybrid_verify import boundary_inputs,verify_samples


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=["screen","build","verify"])
    p.add_argument("--guard",type=int,default=4)
    p.add_argument("--cleanup",choices=["mbuc","coherent"],default="mbuc")
    p.add_argument("--normalization",choices=["barrel","serial"],default="barrel")
    p.add_argument("--output",type=Path)
    p.add_argument("--circuit",type=Path)
    p.add_argument("--random-count",type=int,default=4096)
    a=p.parse_args()
    cfg=HornerConfig(a.guard,a.cleanup,a.normalization)
    if a.command=="screen":
        result=screen(cfg)
    elif a.command=="build":
        if a.output is None:p.error("--output required")
        circ=build(cfg)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        circ.to_file(a.output,format="kmb")
        result=dict(config=vars(cfg),resources=resource_counts(circ),
                    mathematical=screen(cfg),file=str(a.output))
    else:
        import kickmix as km
        if a.circuit is None:p.error("--circuit required")
        circ=km.Circuit.from_file(a.circuit,format="kmb")
        result=verify_samples(circ,boundary_inputs(random_count=a.random_count))
    print(json.dumps(result,indent=2,sort_keys=True))
    if a.output is not None and a.command!="build":
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+"\n")
    if result.get("failures",0):raise SystemExit(1)


if __name__=="__main__":
    main()
