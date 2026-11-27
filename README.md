# Benchmark I — current best fully verified 8.098× circuit (2026-11-27)

The v7 compiler is an **independently authored, don’t-care-aware multioutput Boolean ROM synthesizer**. The base-four-normalized mantissa interval address has 512 encodings, but only 384 are reachable. The 128 unreachable encodings can therefore be assigned arbitrary table values without altering any valid input’s phase. Our reversible ANF decoder assigns these unused entries using the bitwise **XOR extension** `f00=f01 XOR f10 XOR f11`, cancelling the complete highest-order interaction of the two top address bits. It avoids **128 shared nonlinear AND-prefix computations** while preserving every reachable A/B/C coefficient, the original measurement-corrected cleanup, the exact five-qubit residual phase padding and the independently verified four-stage binary normalizer.

- **86 qubits, 2,418 native Toffoli-like operations, 76 fine Z_POW rotations**.
- **Full pinned independent verifier mirror: PASS** across **3,133,452 inputs × 64 = 200,540,928 coherent trajectories**, with all scratch restored.
- Rotation-inclusive **Q×C = 256,476.78036771005**, a **8.097731× reference improvement** and **4.3614% cost reduction versus verified v6**.
- Worst phase error/tolerance remains **0.9687625254666724**, at x=244524.
- Exact `.kmb` SHA-256: `bd691cb587bb51d690d86b4798b1a674fbb61313fa24d02f63ba6422ff4d413a`.
- Source: `src/benchmark1_best_verified_v7.py`, `src/benchmark1_dontcare_boolean_rom.py`, `src/benchmark1_pruned_anf_rom.py`.
- **Rebuild:** `python -m src.benchmark1_best_verified_v7 --output out/best_verified_v7.kmb`.
- Independent result receipt: `results/benchmark1_best_verified_v7.json`.
- [Full pinned mirror PASS](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46014043594) · [Native decoder comparisons](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46013515647).

The previous v6 and v5 builders are preserved and independently SHA-pinned. **20× remains open** (requires Q×C ≤103,844); this is **not an official Kaggle leaderboard submission**.

---

# Historical — previous fully verified 7.745× circuit (2026-11-27)

New **original four-stage binary zero-block normalization** removes the old ten repeated binary exponent increments. The four reversible tests shift by 8, 4, 2, and 1 base-four digit pairs. The four resulting coherent flags directly form the phase exponent; uncomputation reverses the block rotations and uses measurement-conditioned Clifford correction while the original controls are again live. The previously verified five-qubit residual-phase padding, full-domain integer coefficient refit, and p7/q26 arithmetic are preserved.

- **86 qubits, 2,554 Toffoli-like operations, 76 fine rotations.**
- Full pinned independent competition verifier mirror: **PASS, 200,540,928 trajectories** (3,133,452 valid inputs × 64).
- Rotation-inclusive **Q×C = 268,172.78036771005**, approximately **7.744559×** improvement over Google's approximate reference Q×C = 2,076,880; **1.516% lower** Q×C than v5 (272,300.78036771005).
- Same phase-error/tolerance maximum **0.9687625254666724** at x=244524.
- Native KMB SHA-256 `1807c84b9f93008de4ebde2c2231aa8d5f08386dc03e10b63b47999c4577011d`.
- Rebuild with `python -m src.benchmark1_best_verified_v6 --output out/best_verified_v6.kmb`.
- Proof and reproducibility receipt: `results/benchmark1_best_verified_v6.json`.
- [Full pinned verifier PASS](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/46008864751).

The previous full-verified v5 source remains untouched and exact-SHA reproducible. **20× is not yet achieved**; the necessary total cost is at most 103,844. This is *not* an official Kaggle leaderboard submission.

---

# Historical — previous fully verified 7.627× circuit (2026-11-27)

This **independently developed exact phase-padding reduction** uses the original hardware-aware, integer-refitted piecewise-quadratic oracle but eliminates unused padding and corresponding gates from its two residual-dependent variable rotations. Across every legal input, the linear and quadratic phase word each require exactly **five** high-zero padding qubits to avoid wraparound; the constant-phase word retains ten. The original phase scale and arithmetic approximation are unchanged.

- **86 qubits, 2,602 native Toffoli-like gates, 76 fine Z_POW rotations.**
- **Full independently pinned public competition-verifier mirror: PASS**, on **3,133,452 inputs × 64 = 200,540,928 trajectories**.
- Rotation-inclusive effective cost: **3,166.288143810582 Toffoli-equivalent**, Q*C **272,300.78036771005**, a **7.627154× reference improvement** and **3.83772% Q*C reduction versus the previous 7.33× winner**.
- Unchanged worst phase error/tolerance: **0.9687625254666724** at x=244524.
- Exact KMB SHA-256: `4a232127f01b2a901209427e011ba512dd9fde991b627fbc2761399b17203ee4`.
- Independently checked minimum residual padding: **5/5**, all 3,133,452 legal inputs; exact normalized approximation also covers **3,145,728 values**.
- **Rebuild:** `python -m src.benchmark1_best_verified_v5 --output out/best_verified_v5.kmb`
- **Machine-readable receipt:** `results/benchmark1_best_verified_v5.json`
- [Full official-style mirror PASS](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45995733095) · [Exhaustive minimum-padding proof](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45995550913).

The prior `src.benchmark1_best_verified_v4` builder remains retained and SHA-pinned. This is **not an official Kaggle leaderboard submission**. The **20× milestone is not yet achieved**.

---

# Historical — previous fully verified 7.33× circuit (2026-11-26)

This design optimizes **integer coefficients and omitted product columns jointly**, rather than fitting continuous Chebyshev coefficients and truncating afterward. It retains one joint A/B/C lookup, early measurement-corrected A erasure, compact normalization, and uses the lower-cost quadratic multiplier operand orientation.

- **86 qubits, 2,652 native Toffoli-like gates, 86 fine Z_POW rotations.**
- Full pinned independent competition-verifier mirror: **PASS**, **200,540,928** trajectories across **3,133,452** inputs.
- Actual effective non-Clifford gate cost **3,292.650869433569**.
- Verified **Q*C = 283167.97477128694**, **7.334445×** improvement over the approximate Google reference Q*C=2,076,880.
- Worst phase error/tolerance: **0.9687625254666724**, at x=244524.
- Native KMB SHA-256: `0ab11da9ee2ff0b996cdd9bbc89eb1716d2c88bd1f1d934dab3695f7d566ed30`.
- Full verifier evidence: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45857957278
- Exact source + hash-checked reproducibility builder: `python -m src.benchmark1_best_verified_v4 --output out/best_verified_v4.kmb`
- Machine-readable result: `results/benchmark1_best_verified_v4.json`
- **More than 10× (Q*C <207688) has not been achieved**; 30× also remains open.

The result is from a pinned **public mirror verifier**, not a new official Kaggle leaderboard submission. Previous full-PASS 6.86× and 5.65× builders remain in the repository for regressions.

---

# Historical: previous fully verified 6.86× circuit (2026-11-26)

The current best full pinned-mirror verified 22-bit inverse-square-root phase oracle uses hardware-aware sparse p7/q26 coefficient tables (B divisible by 4, C by 16), one joint A/B/C lookup with early A erasure, compact normalization and controlled additions.

- **89 qubits, 2,764 native Toffoli-like operations, 87 fine rotations.**
- Pinned mirror full-domain **PASS** for **200,540,928 trajectories**.
- Effective Toffoli-equivalent cost **3400.461171309322**.
- **Q*C = 302641.04424652964**, improvement **6.8625×** relative to the approximate Google reference 2,076,880.
- Worst phase error / tolerance **0.9172047394125554** (PASS).
- Exact native KMB SHA-256: `0589741be529205d88ebe67d4d9d6a82faebb42debfaf1dcf38c0e39383d09b6`.
- **10× goal Q*C < 207,688 remains unachieved.** 30× goal remains unachieved.
- Full independent verifier: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45833905951
- Rebuild pin: `python -m src.benchmark1_best_joint_v3 --output out/best_joint_v3.kmb`
- Machine-readable result: `results/benchmark1_best_verified_v3.json`

This is a pinned public competition-verifier **mirror**, not an official Kaggle submission. The previous verified 5.65× implementation remains below for exact regression.

---

# Historical: previous fully verified 5.65× result (2026-11-26)

The current best reproducible circuit uses **one joint A/B/C coefficient lookup**
and early measurement-corrected erasure of the consumed A bits. It retains
p7/q26 Chebyshev quadratic phasing with exact cut settings **(10,13,10)**
and compact normalization.

- **112 qubits; 2,658 Toffoli-like operations; 88 fine rotations**
- Full pinned mirror-verifier **PASS** on **200,540,928 trajectories**
- **Q*C = 367,309.6152203282**, approximately **5.65×** better than the
  Google reference (~2,076,880).
- Exact KMB SHA-256:
  `1726e2a438e4296590210c2422004ed1feda9153c09f28a9edc31dce930e0be6`
- **10× objective (Q*C <= 207,688) not yet achieved.**
- **30× objective (Q*C <= 69,229.33) not yet achieved.**

Build the pinned circuit with:
```bash
python -m src.benchmark1_single_lookup_verified --output out/best_invsqrt_v2.kmb
```

The builder rejects any byte-level mismatch against the verified KMB.
[Full verifier PASS](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45828484630) ·
[Independent exact-hash rebuild](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45831536475) ·
[Machine-readable results](results/benchmark1_best_verified_v2.json).

This is independent mirror verification, **not a new official Kaggle
leaderboard submission**. The previously verified 4.83× builder remains
available at `src/benchmark1_best.py` for regression and comparison.

---

## Historical reference — previous fully verified 4.83× circuit


**Active branch:** `benchmark-1`.

**Previously best fully verified circuit:** grouped Chebyshev quadratic 7/26,
cut settings `(10,12,11)`, coefficient grouping **A then BC**,
with **compact measurement-corrected barrel normalization**.
Uses **112 qubits**, **3,218 Toffoli-like operations**, and
**88 fine rotations**. Full pinned competition-verifier mirror:
**PASS** on **200,540,928 trajectories**, `Q*C = 430,259.61862786254`
(~**4.83×** improvement versus Google's ~2,076,880 reference).

**30× target:** `Q*C ≤ 69,229.33` — **not achieved**.
The best verified implementation still needs a further **6.21×**
reduction in cost.

Rebuild the exact, SHA-256-checked candidate using pinned native kickmix:

```bash
python -m src.benchmark1_best --output out/best_invsqrt.kmb
```

The builder explicitly checks native qubit and gate counts **and** the
serialized circuit's SHA-256. If the hash changes under a different kickmix
release, the rebuilt binary must be reverified before being submitted.

Evidence: [full verification PASS](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42657329801).
See [measurement during computation: native tests and cost limits](results/measurement_during_computation.md),
[XAG / Horner / measurement cleanup research](results/xag_horner_mbuc_20261109.md)
and [architecture comparison](results/architecture_comparison_20261109.md).

The real Kaggle leaderboard still needs an actual new submission.
Do **not** submit exploratory XAG, Horner, or linear candidates; their
current end-to-end native resource costs are higher or incorrect.

---

# Benchmark I — Priority 1 optimization experiments

**Active research branch:** `benchmark-1`. For current results, algorithms,
validation status, and the explicit 30× reference-cost stretch objective, see
[Priority 1 report](results/priority1_v2_report.md).

**Reproduce newer experiments** (requires `numpy`, `scipy`, and the
native C++20 `kickmix` extension from Qualtran):

```bash
python -m pytest -q tests/test_priority1_model.py tests/test_priority1_exact.py
python -m src.priority1_cli build --truncated --family chebyshev \
  --degree 2 --segment-bits 7 --precision 26 --drops 10,11,12 \
  --output out/chebyshev_7_26.kmb
python -m src.priority1_cli verify --circuit out/chebyshev_7_26.kmb --full
```

The 30× target is **not yet achieved**. Native sampled verification,
complete-domain mathematical screening, and full official-style resource
verification have different evidentiary strength; consult the report.

---


## Current focus: normalized hybrid direct phase

The original dyadic parity-gadget pilot is retained as a comparison. The
stronger `src/hybrid_*.py` experiment implements a **full 22-bit** reversible
kickmix circuit: normalize, lookup piecewise quadratic phase coefficients,
apply constant/linear/quadratic phase terms, then restore every register.
There is no Newton iteration or assembled numeric inverse-square-root output.

Numerical/physical results: [hybrid_report.md](results/hybrid_report.md).
The numerical sweep covers every one of the 3,133,452 integer inputs.
Native kickmix tests construct full-width .kmb circuits and check randomized
measurement trajectories. GitHub Actions additionally contains a three-way
full-domain run of 200,540,928 trajectories per candidate.

    python -m pip install numpy pytest
    python -m pytest -q tests/test_hybrid_model.py
    python -m src.hybrid_cli screen --output runs/hybrid_numerical.json

Install the **native kickmix C++ extension** from the upstream Qualtran
`kickmix/` project, then run:

    python -m src.hybrid_cli build --segment-bits 8 --precision 26 --output out/hybrid_8_26.kmb
    python -m src.hybrid_cli verify out/hybrid_8_26.kmb --random-count 2048
    python -m src.hybrid_cli verify out/hybrid_8_26.kmb --full --output runs/hybrid_8_26_full.json

Use `bash scripts/verify_official.sh out/hybrid_8_26.kmb official` with the
Kaggle-provided verifier and dataset before submission. Our development
verifier is **not** the official scorer, and its reported Toffoli count
**excludes** analytic fine-rotation synthesis charges.

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
