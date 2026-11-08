"""CLI: numerical screen, native kickmix build, and development verifier."""
import argparse
import json
from pathlib import Path
from .direct_phase import MAX_X, INPUT_BITS, dyadic_partition, verify_partition
from .kickmix_phase import build_circuit, estimated_rotations
from .verify import verify_kmb


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("screen")
    s.add_argument("--degrees", type=int, nargs="+", default=[1, 2, 3])
    s.add_argument("--max-x", type=int, default=MAX_X)
    s.add_argument("--safety", type=float, default=.75)
    s.add_argument("--output", type=Path)
    b = sub.add_parser("build")
    b.add_argument("--degree", type=int, choices=[1, 2, 3], default=2)
    b.add_argument("--max-x", type=int, default=MAX_X)
    b.add_argument("--input-bits", type=int, default=INPUT_BITS)
    b.add_argument("--output", type=Path, required=True)
    v = sub.add_parser("verify")
    v.add_argument("circuit")
    v.add_argument("--max-x", type=int, default=MAX_X)
    v.add_argument("--shots-per-x", type=int, default=64)
    v.add_argument("--development", action="store_true")
    a = ap.parse_args()
    if a.cmd == "screen":
        rows = []
        for degree in a.degrees:
            blocks = dyadic_partition(degree, max_x=a.max_x, safety=a.safety)
            report = dict(degree=degree, **verify_partition(blocks, max_x=a.max_x),
                          estimated_parity_rotations=estimated_rotations(blocks))
            rows.append(report)
            print(json.dumps(report))
        if a.output:
            a.output.parent.mkdir(parents=True, exist_ok=True)
            a.output.write_text(json.dumps(rows, indent=2)+"\n")
    elif a.cmd == "build":
        if a.input_bits < a.max_x.bit_length():
            ap.error("Input register too narrow")
        blocks = dyadic_partition(a.degree, max_x=a.max_x)
        circuit = build_circuit(blocks, input_bits=a.input_bits)
        a.output.parent.mkdir(parents=True, exist_ok=True)
        circuit.to_file(a.output, format="kmb")
        stats = circuit.operation_statistics()
        print(json.dumps(dict(output=str(a.output), blocks=len(blocks),
                              qubits=circuit.num_qubits,
                              toffoli_like=stats.num_toffoli_like,
                              t_like=stats.num_t_like,
                              z_pow=stats.num_z_pow())))
    else:
        print(json.dumps(verify_kmb(a.circuit, max_x=a.max_x,
                                   shots_per_x=a.shots_per_x,
                                   competition_contract=not a.development),
                         indent=2))


if __name__ == "__main__":
    main()
