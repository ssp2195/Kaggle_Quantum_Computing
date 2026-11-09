#!/usr/bin/env python3
"""Package an optimized Benchmark I circuit plus official reference files.

The competition requires all three tasks in one ZIP. Only invsqrt is replaced.
Never put authentication tokens into this file or a GitHub repository.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

def create_submission(optimized: Path, reference: Path, destination: Path) -> list[str]:
    if not optimized.is_file() or optimized.suffix not in {".kmb", ".kmx"}:
        raise ValueError("Optimized inverse-square-root circuit must be a .kmb or .kmx file")
    root = reference / "sample_submission" if (reference / "sample_submission").is_dir() else reference
    dqi = next((root / f"dqi.{extension}" for extension in ("kmb", "kmx")
                if (root / f"dqi.{extension}").is_file()), None)
    if dqi is None:
        raise FileNotFoundError("Missing official sample_submission/dqi.kmb or .kmx")
    df = root / "dfthcsos"
    if not (df / "submission.npz").is_file():
        raise FileNotFoundError("Missing official sample_submission/dfthcsos/submission.npz")
    angles = []
    i = 0
    while True:
        options = [df / f"angle_k_{i}.{ext}" for ext in ("kmb", "kmx")]
        found = next((p for p in options if p.is_file()), None)
        if found is None:
            break
        angles.append(found)
        i += 1
    if not angles:
        raise FileNotFoundError("Missing dfthcsos/angle_k_0.kmb or .kmx")
    unexpected = [p.name for p in df.glob("angle_k_*.*")
                  if p.suffix in {".kmb", ".kmx"}
                  and p.stem not in {f"angle_k_{j}" for j in range(len(angles))}]
    if unexpected:
        raise ValueError(f"Noncontiguous or duplicate ANGLE indices: {unexpected}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    files = [(optimized, f"invsqrt{optimized.suffix}"),
             (dqi, f"dqi{dqi.suffix}"),
             (df / "submission.npz", "dfthcsos/submission.npz")]
    files.extend((p, f"dfthcsos/{p.name}") for p in angles)
    with ZipFile(destination, "w", ZIP_DEFLATED) as zipf:
        for src, arcname in files:
            zipf.write(src, arcname)
    with ZipFile(destination) as zipf:
        assert zipf.testzip() is None
        assert sorted(zipf.namelist()) == sorted(x[1] for x in files)
    return [x[1] for x in files]

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--invsqrt", type=Path, required=True, help="Optimized 8/26 .kmb")
    p.add_argument("--reference", type=Path, required=True,
                   help="Official competition-data directory or sample_submission directory")
    p.add_argument("--output", type=Path, default=Path("out/submission.zip"))
    args = p.parse_args()
    members = create_submission(args.invsqrt, args.reference, args.output)
    print(f"Created {args.output} with {len(members)} files:")
    for m in members:
        print(" ", m)

if __name__ == "__main__":
    main()
