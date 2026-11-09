# Joint coefficient / reversible-truncation optimization (Benchmark I)

## Objective

Target a **30x** reduction from the published Google Benchmark I reference
cost approximately 2,076,880, i.e. **Q*C <= 69,229.3**. If the other two
benchmarks remain at their reference scores, this corresponds to the
**estimated** combined Kaggle score `30**(1/3) = 3.1072`.

The user's actual reported combined score for the original 8/26 hybrid is
**1.33**, and no experiment in this file is an actual Kaggle leaderboard
submission.

## Why this experiment is different

Previously we fitted continuous Chebyshev polynomials, rounded coefficients,
and *then* truncated the arithmetic. That can be suboptimal: omission of
low partial products changes the implemented polynomial and its worst-case
error. Here we search the **integer coefficients** against the **exact
omitted-partial-product arithmetic implemented by `_mul_cut`**.

- Mathematical form: `2**(-k) * (A - truncated(B*u) + truncated(C*truncated(u*u)))`
  with each term scaled by its exact fixed-point denominator.
- Piecewise quadratic, `p=7` segment bits, `q=25` fractional coefficient
  bits; `384` active normalized intervals.
- Drop settings (linear, square, quadratic): `(9,10,9)`.
- For each normalized interval, inspect all `8192` discrete mantissa
  offsets and minimize the maximum weighted relative phase error via
  deterministic integer-coordinate search (not a claim of global minimax).
- Coefficient selection is reproducible, with no trained model, random seed,
  private dataset, or new dependency beyond NumPy.

The preliminary exhaustive **normalized-mantissa** sweep found a maximum
error ratio approximately **0.761882**, with no intervals exceeding the
tolerance. This is a **numerical screening result**, not a native circuit
verification. The explicit `full_cut_screen` test also evaluates all
3,133,452 admissible Benchmark I integers; GitHub Actions runs the test.

## Commands

    python -m pytest -q tests/test_priority1_integer_tuned.py

    # Requires actual native kickmix C++20 extension:
    python -m src.priority1_cli build --truncated \
        --family chebyshev-cut-tuned --degree 2 \
        --segment-bits 7 --precision 25 \
        --drops 9,10,9 --output out/cut_tuned_7_25.kmb

    python -m src.priority1_cli verify \
        --circuit out/cut_tuned_7_25.kmb --random-count 8192

The native unit test checks serializing and reloading .kmb and validating
input restoration, zero scratch, and phase over 64 trajectories per
representative input.

[Native candidate workflow](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/workflows/benchmark1-integer-tuned.yml)

**No claim of a 30x gain:** even a 149-qubit circuit with 2,961 native
Toffoli-like gates already has a `Q*Toffoli` cost floor of 441,189,
before rotation synthesis. It is >6x above the target cost; coefficient
tuning alone cannot close this gap. A fundamentally different direct
phase representation with vastly fewer nonlinear gate operations or
workspace qubits would be required.

## Native kickmix outcome (completed)

[Native CI run 42632099404](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42632099404)
**passed**. It constructed and reloaded the real C++20 kickmix 22-bit
circuit, verified 64 measurement trajectories per selected input, and
recorded:

- Qubits: **143**
- Native Toffoli-like gates: **3,155**
- Static fine Z_POW rotations: **90**
- Exhaustive 3,133,452-input mathematical phase screen: **PASS**
- Worst mathematical error/tolerance: **0.7617151634552422**

The computed lower bound Q * Toffoli-like = **451,165** is
already **6.52 times the 30x target cost** of 69,229, even
pretending rotations are free. This variant therefore cannot
meet the 30x objective.

[Full pinned competition verifier](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42632573053)
was launched separately for 200,540,928 trajectories and the rotation-
inclusive resource score; do not report a verified Q*C until its PASS
and measured score are available.
