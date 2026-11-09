"""Normalization round-trip is the identity, not an inverse-sqrt oracle."""
import pytest
km=pytest.importorskip("kickmix")
from src.priority1_normalize import normalize_barrel,denormalize_barrel
from src.hybrid_verify import boundary_inputs


def test_normalizer_roundtrip(tmp_path):
    builder=km.CircuitBuilder()
    x=builder.create_quantum_register(22,name="x")
    preds,scale=normalize_barrel(builder,x)
    denormalize_barrel(builder,x,preds,scale)
    circuit=builder.finish_circuit()
    scratch=km.array(sorted(set(km.q(i) for i in range(circuit.num_qubits))
                            -set(circuit.register_data["x"])))
    samples=boundary_inputs(random_count=4096)
    sim=km.Simulator(batch_size=256)
    sim.use_same_registers_as(circuit)
    for offset in range(0,len(samples),256):
        batch=samples[offset:offset+256]
        sim.clear_for_shot()
        for slot,val in enumerate(batch):
            sim.write_within_shot("x",slot,val)
        sim.do(circuit)
        for slot,val in enumerate(batch):
            assert sim.read_within_shot("x",slot,out=int)==val
            assert sim.read_within_shot(scratch,slot,out=int)==0
            # Phase is reported in half-turns. A round-trip is identity.
            phase=(float(sim.read_shot_phase(slot))/2)%1
            assert min(phase,1-phase)<1e-12, (val,phase)
