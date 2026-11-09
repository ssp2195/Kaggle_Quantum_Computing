"""Build and test architecture prototypes against kickmix."""
import argparse
import json
from pathlib import Path

from .benchmark1_piecewise_linear import LinearConfig, build as build_linear, screen
from .priority1_grouped import build_grouped
from .priority1_model import Approximation
from .hybrid_kickmix import resource_counts
from .hybrid_verify import verify_samples, boundary_inputs


def build_candidate(name):
    if name in ("linear-joint", "linear-streamed"):
        return build_linear(LinearConfig(10, 28, 7), grouping=name.split("-")[1])
    if name in ("square-rows", "square-symmetric"):
        return build_grouped(Approximation("chebyshev", 2, 7, 26),
                             drops=(10, 12, 11), grouping="ab_c",
                             normalization="barrel", square_method=name.split("-")[1])
    raise ValueError(f"Unknown architecture: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["screen", "build", "verify"])
    parser.add_argument("--architecture",
                        choices=["linear-joint", "linear-streamed",
                                 "square-rows", "square-symmetric"],
                        default="linear-joint")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--circuit", type=Path)
    parser.add_argument("--random-count", type=int, default=2048)
    args = parser.parse_args()
    if args.command == "screen":
        result = screen(LinearConfig(10, 28, 7))
    elif args.command == "build":
        if not args.output:
            parser.error("--output required to build")
        circuit = build_candidate(args.architecture)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        circuit.to_file(args.output, format="kmb")
        result = dict(architecture=args.architecture, file=str(args.output),
                      **resource_counts(circuit))
    else:
        import kickmix as km
        if args.circuit is None:
            parser.error("--circuit required to verify")
        circuit = km.Circuit.from_file(args.circuit, format="kmb")
        result = verify_samples(circuit, boundary_inputs(random_count=args.random_count))
    print(json.dumps(result, indent=2))
    if args.output and args.command != "build":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    if result.get("failures", 0):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
