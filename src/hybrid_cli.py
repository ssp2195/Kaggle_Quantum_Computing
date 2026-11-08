"""Normalized quadratic hybrid-phase research CLI."""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
from .hybrid_model import MAX_X, Config, full_domain_screen
from .hybrid_kickmix import build_hybrid, resource_counts
from .hybrid_verify import verify_samples, boundary_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("screen")
    p.add_argument("--pairs", nargs="*", default=["7:30", "8:24", "8:26", "8:28", "9:26", "9:28", "9:30", "10:28"])
    p.add_argument("--output", type=Path)
    p = sub.add_parser("build")
    p.add_argument("--segment-bits", type=int, default=8)
    p.add_argument("--precision", type=int, default=28)
    p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("verify")
    p.add_argument("circuit", type=Path)
    p.add_argument("--full", action="store_true")
    p.add_argument("--random-count", type=int, default=256)
    p.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.command == "screen":
        rows = []
        for pair in args.pairs:
            p, q = map(int, pair.split(":"))
            result = full_domain_screen(Config(p, q))
            rows.append(result)
            print(json.dumps(result), flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(rows, indent=2) + "\n")
        return
    import kickmix as km
    if args.command == "build":
        circuit = build_hybrid(Config(args.segment_bits, args.precision))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        circuit.to_file(args.output, format="kmb")
        print(json.dumps(dict(output=str(args.output), **resource_counts(circuit))))
        return
    circuit = km.Circuit.from_file(args.circuit, format="kmb")
    inputs = range(1, MAX_X + 1) if args.full else boundary_inputs(random_count=args.random_count)
    tic = time.monotonic()
    report = verify_samples(circuit, inputs, print_progress=args.full)
    report["elapsed_seconds"] = time.monotonic() - tic
    print(json.dumps(report, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    if report["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
