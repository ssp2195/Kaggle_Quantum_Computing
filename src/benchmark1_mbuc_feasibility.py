"""Feasibility of staged measurement-based AND uncomputation in kickmix.

A cubic Boolean temporary w=(a AND b) AND c can be cleaned for ZERO
additional Toffolis *provided* the intermediate t=a AND b remains live:
  measure w in X basis and conditional CZ(t,c)
  measure t in X basis and conditional CZ(a,b)
Both corrections are Clifford.

This does NOT compute AND for free: the forward two CCX remain necessary.
If t is gone, the correction for measuring w directly can involve the
cubic phase (-1)^(mu*a*b*c), which a CZ on original input pairs cannot
generally cancel without recomputing a helper or paying a non-Clifford cost.
"""
import argparse
import json


def chain_circuit(cleanup="mbuc"):
    if cleanup not in ("mbuc","coherent"):
        raise ValueError("cleanup must be mbuc/coherent")
    import kickmix as km
    b=km.CircuitBuilder()
    x=b.create_quantum_register(3,name="x")
    t,w=b.alloc_qubits(2)
    b.init_and(x[0],x[1],t)
    b.init_and(t,x[2],w)
    if cleanup=="mbuc":
        b.del_and(t,x[2],w)
        b.del_and(x[0],x[1],t)
    else:
        b.ccx(t,x[2],w)
        b.ccx(x[0],x[1],t)
    b.free(km.array([t,w]))
    return b.finish_circuit()


def compare():
    import kickmix as km
    from .hybrid_kickmix import resource_counts
    result={}
    for mode in ("mbuc","coherent"):
        c=chain_circuit(mode)
        sim=km.Simulator(batch_size=64)
        sim.use_same_registers_as(c)
        for x in range(8):
            sim.clear_for_shot()
            for s in range(64):
                sim.write_within_shot("x",s,x)
            sim.do(c)
            for s in range(64):
                assert sim.read_within_shot("x",s,out=int)==x
                assert abs(float(sim.read_shot_phase(s)))<1e-12
        result[mode]=dict(**resource_counts(c),
                          shots_verified=8*64)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output")
    a=p.parse_args()
    result=compare()
    print(json.dumps(result,indent=2))
    if a.output:
        from pathlib import Path
        x=Path(a.output)
        x.parent.mkdir(parents=True,exist_ok=True)
        x.write_text(json.dumps(result,indent=2)+"\n")
