"""Resource gates for genuine p6 quadratic in-place-Horner architectures.

Only test native resource bounds for candidates passing full classical math;
do not attribute winner status without independent quantum verifier.
"""
import json
from .benchmark1_horner_xag import XAGHornerConfig,screen,build
from .benchmark1_centered_horner import CenteredConfig,screen as centered_screen,build as centered_build
from .hybrid_kickmix import resource_counts

REFERENCE_QC=2076880.
CURRENT_QC=272300.78036771005
TARGET20=REFERENCE_QC/20

def main():
 rows=[]
 configs=[]
 for guard in (4,5,6):
  for lookup in ("qrom","xag"):
   configs.append(("horner",XAGHornerConfig(guard=guard,lookup=lookup,cleanup="mbuc")))
 for p in (6,7):
  for guard in (3,4,5):
   for lookup in ("qrom","xag"):
    configs.append(("centered",CenteredConfig(segment_bits=p,guard=guard,lookup=lookup,cleanup="mbuc")))
 for fam,conf in configs:
  try:
   math=screen(conf) if fam=="horner" else centered_screen(conf)
   report=dict(family=fam,parameters=str(conf),math=math)
   if math["passes"]:
    c=build(conf) if fam=="horner" else centered_build(conf)
    s=resource_counts(c)
    report["native"]=s
    report["can_reach_20x"]=s["lower_bound_qc"]<TARGET20
    report["can_beat_incumbent"]=s["lower_bound_qc"]<CURRENT_QC
  except Exception as e:
   report=dict(family=fam,parameters=str(conf),exception=type(e).__name__+":"+str(e)[:400])
  print(json.dumps(report),flush=True)
  rows.append(report)
 import pathlib
 p=pathlib.Path("runs/architecture_horner_sweep.json");p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(rows,indent=2))
if __name__=="__main__":main()
