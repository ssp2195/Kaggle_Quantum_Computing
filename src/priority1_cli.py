"""Priority 1: compile and measure competing direct phase approximations."""
import argparse
import json
from pathlib import Path
from .priority1_model import Approximation, full_domain_screen
from .priority1_kickmix import build_priority1
from .priority1_truncated import build_truncated, plan
from .priority1_exact import full_cut_screen
from .priority1_streamed import build_streamed
from .priority1_grouped import build_grouped
from .hybrid_kickmix import resource_counts
from .hybrid_verify import verify_samples, boundary_inputs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",choices=["screen","build","verify"])
    p.add_argument("--family",choices=["taylor","chebyshev","minimax","chebyshev-tuned"],default="chebyshev")
    p.add_argument("--degree",type=int,default=2)
    p.add_argument("--segment-bits",type=int,default=6)
    p.add_argument("--precision",type=int,default=26)
    p.add_argument("--output",type=Path)
    p.add_argument("--circuit",type=Path)
    p.add_argument("--full",action="store_true")
    p.add_argument("--truncated",action="store_true")
    p.add_argument("--normalization",choices=["serial","barrel"],default="serial")
    p.add_argument("--grouping",choices=["ab_c","a_bc"],
                   help="two-stage coefficient lookup grouping")
    p.add_argument("--streamed",action="store_true",
                   help="load and erase each polynomial coefficient independently")
    p.add_argument("--drops",type=str,
                   help="exhaustively screen comma-separated linear,square,quadratic drop counts")
    p.add_argument("--truncation-budget",type=float,default=.025,
                   help="maximum relative phase-error allowance per discarded product")
    p.add_argument("--random-count",type=int,default=2048)
    a=p.parse_args()
    if a.command=="screen":
        result=full_domain_screen(
            Approximation(a.family,a.degree,a.segment_bits,a.precision))
    else:
        import kickmix as km
        if a.command=="build":
            if a.output is None:
                p.error("--output required")
            cfg=Approximation(a.family,a.degree,a.segment_bits,a.precision)
            drops=tuple(map(int,a.drops.split(","))) if a.drops else None
            if a.grouping:
                if drops is None:
                    p.error("--grouping requires --drops")
                circuit=build_grouped(cfg,drops,grouping=a.grouping,
                                      normalization=a.normalization)
            elif a.streamed:
                if drops is None:
                    p.error("--streamed requires --drops")
                circuit=build_streamed(cfg,drops,normalization=a.normalization)
            elif a.truncated:
                circuit=build_truncated(cfg,tolerance_fraction=a.truncation_budget,
                                        explicit_drops=drops)
            else:
                circuit=build_priority1(cfg)
            a.output.parent.mkdir(parents=True,exist_ok=True)
            circuit.to_file(a.output,format="kmb")
            result=dict(**resource_counts(circuit),file=str(a.output),
                        architecture=("grouped_"+a.grouping if a.grouping else
                                      "streamed" if a.streamed else
                                      "truncated" if a.truncated else "full_product"))
            if a.truncated or a.streamed or a.grouping:
                result["normalization"]=a.normalization
                result["truncation_plan"]=(full_cut_screen(cfg,drops) if drops is not None
                                            else plan(cfg,a.truncation_budget))
        else:
            if a.circuit is None:
                p.error("--circuit required")
            circuit=km.Circuit.from_file(a.circuit,format="kmb")
            xs=range(1,3_133_453) if a.full else boundary_inputs(random_count=a.random_count)
            result=verify_samples(circuit,xs,print_progress=a.full)
    print(json.dumps(result,indent=2,sort_keys=True))
    if a.output and a.command!="build":
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+"\n")
    if result.get("failures",0):
        raise SystemExit(1)

if __name__=="__main__":
    main()
