# Kaggle Quantum Computing — Benchmark I direct-phase pilot

Experimental branch for the Google Quantum AI Verifiable Quantum Circuit Compilation competition. The original main branch is unchanged.

**Target:** apply exp(2πi/sqrt(x)) to each basis state of the named 22-qubit input register x, for 1 <= x <= 3,133,452. Preserve x, restore all scratch to zero, and meet tolerance 2^-24/sqrt(x) turns. Official verification uses 64 shots per input (200,540,928 trajectories). Score is Q times *official* Toffoli-equivalent cost including non-Clifford rotation synthesis. The authoritative requirements are in Section III.A.3 of arXiv:2610.23832 and the Kaggle competition.

This branch compares piecewise linear, quadratic and cubic Chebyshev direct-phase approximations on dyadically aligned ranges. A separate optional native kickmix builder implements the phase as conditional parity gadgets, without constructing an inverse-square-root output register.

## Setup

Requires Python 3.11+, numpy and pytest. To generate or simulate native circuits install **kickmix** from the kickmix/ subdirectory of https://github.com/quantumlib/Qualtran following that project's build instructions. No stub package is substituted. The mathematical experiments run without kickmix.

Commands (from the repository root):

    python -m pip install numpy pytest
    python -m pytest -q
    python -m src.cli screen --output results/screen.json

Native experimental circuit on a **small development domain**:

    python -m src.cli build --degree 2 --max-x 31 --input-bits 5 --output out/toy.kmb
    python -m src.cli verify out/toy.kmb --max-x 31 --shots-per-x 64 --development

Full candidate (expensive and **not yet validated**):

    python -m src.cli build --degree 2 --output out/invsqrt_experiment.kmb
    python -m src.cli verify out/invsqrt_experiment.kmb

The full development verifier checks exactly one named register x, its width, all 3,133,452 allowed inputs, input restoration, all unused workspace, phases modulo one turn, and 64 shots per input. **Run the separate official Kaggle verifier before submission**: this development verifier does not calculate official resource cost, reference fixed-point phase values, salted SHA-256 sampling or the official rotation fallback synthesis. Do not submit using this verifier alone.

## Findings

See results/findings.md and results/numerical_screen.json. All three approximations pass the numerical tolerance on the complete integer domain, but straightforward parity decompositions appear much more expensive in fine Z_POW rotations than the reference. Numerical success is **not** a full kickmix or Kaggle verification result.

Code consists of a numerical model, native circuit builder, development verifier, CLI and tests. No RL or extraneous framework has been added.


## Authoritative Kaggle verification (competition files not redistributed)

Once you have accepted the Kaggle rules and configured the Kaggle CLI, download
the **official** data and use the verifier shipped with the competition:

    mkdir -p downloads official
    kaggle competitions download google-verifiable-quantum-circuit-compilation -p downloads
    python -m zipfile -e downloads/google-verifiable-quantum-circuit-compilation.zip official
    bash scripts/verify_official.sh out/invsqrt_experiment.kmb official

For a fast *non-authoritative* 64-input smoke run, use
`MAX_X=64 MAX_WORKERS=1 bash scripts/verify_official.sh ...`.
Use the default full range of 3,133,452 inputs and 64 trajectories each
before declaring the candidate valid. Inspect the verifier's JSON `status`,
`verified_shots`, `num_qubits`, `toffoli_equivalent` and `score`.
Do not commit downloaded competition data, credentials, or generated .kmb files.
The verifier and baseline belong to the competition, not this repository.

Native smoke tests are wired to GitHub Actions and install a current upstream
kickmix build; if the upstream API changes, inspect the workflow logs.
The competition's private test salt cannot be reproduced locally.
