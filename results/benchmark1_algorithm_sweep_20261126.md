# Benchmark I architecture and algorithm experiments — 26 November 2026

## Promoted winner

Cost-aware integer phase synthesis simultaneously optimizes approximation
coefficients against the *exact omitted partial products*, then selects a
hardware-efficient multiplication orientation. p7/q26 Chebyshev-based
coefficient seeds are only initialization, not the actual optimization
objective. The exact 384-interval x 8192-offset normalized model is
verified before circuit compilation. QROM loading/unloading remains
coherent, with measurement-corrected cleanup only for valid temporary
predicates.

- **Full pinned independent mirror PASS** on 3,133,452 integers x 64
  trajectories = 200,540,928, run:
  https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/45857957278
- **86 qubits**, **2,652 native Toffoli-like gates**, **86 fine rotations**.
- **Q*C 283167.97477128694**, reference improvement **7.334445x**.
- Worst normalized error/tolerance 0.9687625254666724.
- Exact KMB SHA-256 `0ab11da9ee2ff0b996cdd9bbc89eb1716d2c88bd1f1d934dab3695f7d566ed30`.
- This was **not** an official Kaggle leaderboard submission.

Previous full PASS: 6.8625x Q*C302641.04, run 45833905951;
reoriented multiplication: 6.9318x Q*C299615.04,
run 45857288125.

## Negative architectural results from isolated research branches

| Candidate | Representative native Q/T/fine rotations | Decision |
| --- | --- | --- |
| p6 sparse interpolation, cuts (11,12,13) | 94 / 2866 / 91 | full numerical and native sampled PASS, Q*T floor too high |
| p8 sparse interpolation, padded coherent-corrected unlookup | 84 / 3030 / 86 | sampled PASS, resource bound uncompetitive |
| p9 sparse interpolation | 82 / 3804 / 81 | sampled PASS, resource bound uncompetitive |
| Fused linear/quadratic carry accumulator | 104 / 4196 / 64 | much higher full cost floor |
| Factored high/low quadratic residual | 89 / 3124 / 88 | higher cost floor |
| Direct quadratic monomial accumulator | 85 / 5214 / 91 | too many adders |
| In-place two-stage Horner | 115 / 4293 / 37 | excess arithmetic gates |
| Upper-triangle reversible square | 86 / 3026 / 87 | qubit reduction offset by extra gates |
| Piecewise-linear p10/q27 approximation | numerical error >2.3x tolerance on sample | rejected before native build |
| Stronger p7 coefficient-bit zeroing B divisible by 8 | numerical error >1.25x tolerance on sample | rejected before native build |

Lower cost bounds use Q times native Toffolis and **do not** include
fine-rotation synthesis. Passing a cost floor is insufficient for promotion;
the highest verified full Q*C determines the incumbent.

The best winner's native floor **86×2652 = 228072**, already greater
than the **10x target of 207688**. Therefore neither 10x nor 30x
is achieved by this implemented circuit; a more fundamental
Boolean-phase representation/space scheduling change is required.

The experimental branches are intentionally disposable after preservation
of this report and the exact pinned v4 source on `benchmark-1`.
