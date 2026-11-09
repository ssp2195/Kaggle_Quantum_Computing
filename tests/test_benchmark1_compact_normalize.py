"""Measurement-based prefix erasure must preserve all superposition phases."""
import pytest


def test_native_compact_normalization_roundtrip():
    km=pytest.importorskip("kickmix")
    from src.benchmark1_compact_normalize import (
        normalize_compact,denormalize_compact)
    from src.priority1_normalize import normalize_barrel,denormalize_barrel
    from src.hybrid_kickmix import resource_counts
    for mode in ("original","compact"):
        b=km.CircuitBuilder()
        x=b.create_quantum_register(22,name="x")
        flags,shift=((normalize_compact(b,x) if mode=="compact"
                      else normalize_barrel(b,x)))
        if mode=="compact":
            denormalize_compact(b,x,flags,shift)
        else:
            denormalize_barrel(b,x,flags,shift)
        c=b.finish_circuit()
        print("COMPACT_NORM",mode,resource_counts(c))
        sim=km.Simulator(batch_size=64)
        sim.use_same_registers_as(c)
        for xv in (1,2,3,4,5,16,31,64,255,256,1024,65535,1048576,
                   2097151,3133452):
            sim.clear_for_shot()
            for shot in range(64):
                sim.write_within_shot("x",shot,xv)
            sim.do(c)
            for shot in range(64):
                assert sim.read_within_shot("x",shot,out=int)==xv
                assert abs(float(sim.read_shot_phase(shot)))<1e-10
                assert sim.read_shot_phase(shot) is not None


def test_native_full_phase_compact_normalizer():
    km=pytest.importorskip("kickmix")
    from src.priority1_grouped import build_grouped
    from src.priority1_model import Approximation
    from src.hybrid_kickmix import resource_counts
    from src.hybrid_verify import verify_samples,boundary_inputs
    c=build_grouped(Approximation("chebyshev",2,7,26),
                    (10,12,11),grouping="a_bc",normalization="compact")
    res=resource_counts(c)
    report=verify_samples(c,boundary_inputs(random_count=8192))
    print("COMPACT_ORACLE_RESOURCES",res)
    print("COMPACT_ORACLE_SAMPLE",report)
    assert report["failures"]==0,report
