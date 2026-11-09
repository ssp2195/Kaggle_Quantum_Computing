# Benchmark I — Priority 1: approximation + direct phase architecture

Research branch: `benchmark-1`. All results below were obtained using actual
native kickmix C++ circuit construction and operation statistics, except
where expressly labeled mathematical-only.

## Objective and honest 30x target

- Current Google reference cost: approximately **2,076,880**.
- First verified hybrid 8/26 cost: **876,513.3341652812**, PASS on the
  pinned public competition-verifier mirror for all **200,540,928** trajectories.
- Desired **30x improvement vs reference**: benchmark cost <=
  **69,229.33**, corresponding to an estimated *combined* score of
  `30**(1/3) ~= 3.107` if Benchmarks II and III are unchanged.
- No result documented below achieves 30x. Do not report sampled tests
  or static operation counts as an official cost.

## Full-domain mathematical approximation screens

For the normalized inverse-square-root function, quantized polynomial
coefficients use `2**-q` fractional units. The simulator uses actual
bitwise integer arithmetic and `Z_POW` phase operations.

| Family | Degree | Segment bits | Coefficient precision | Intervals | Worst error / tolerance | Fails |
|---|---:|---:|---:|---:|---:|---:|
| Taylor (baseline) | 2 | 8 | 26 | 768 | 0.617764 | 0 |
| Chebyshev interpolation | 2 | 6 | 26 | 192 | 0.735365 | 0 |
| Chebyshev interpolation | 2 | 7 | 26 | 384 | 0.600155 | 0 |
| Chebyshev interpolation | 3 | 4 | 26 | 48 | 0.734057 | 0 |
| Quantization-tuned Chebyshev | 2 | 6 | 25 | 192 | 0.859530 | 0 |
| Weighted relative minimax (rounded coefficients) | 2 | 7 | 26 | 384 | 0.651389 | 0 |
| Weighted relative minimax (rounded coefficients) | 2 | 6 | 26 | 192 | 1.651256 | 84,361 |

**Important:** Chebyshev interpolation and weighted minimax are distinct.
Integer rounding can make a continuous minimax solution worse. Seven
one-unit changes to Chebyshev coefficients allow the 6/25 candidate to
pass every integer in the finite domain.

## Actual native kickmix resource comparison

These are full 22-bit-input .kmb circuits. Fine rotations are the static
count of class >=4 rotations. Native sample tests reported zero failures
for the candidates below.

| Candidate | Q | Toffoli-like | Fine Z_POW |
|---|---:|---:|---:|
| Original Taylor 8/26 | 165 | 4490 | 116 executed |
| Chebyshev quadratic 6/26 full product | 189 | 4737 | 127 |
| Chebyshev quadratic 7/26 full product | 177 | 4441 | 121 |
| Chebyshev cubic 4/26 full product | 280 | 10720 | 210 |
| Quantization-tuned Chebyshev quadratic 6/25 full product | 183 | 4551 | 124 |
| Weighted minimax quadratic 7/26 full product | 178 | 4495 | 122 |
| Chebyshev 6/26 conservatively truncated | 161 | 3623 | 100 |
| Quantization-tuned 6/25 conservatively truncated | 156 | 3561 | 99 |
| Chebyshev 7/26 conservatively truncated | 153 | 3529 | 97 |
| Chebyshev 7/26 medium truncation | 152 | 3361 | 94 |
| Chebyshev 7/26 aggressive bounded truncation | **151** | **3287** | **93** |

The key observation is that **the smaller polynomial table alone does
not produce a cheaper circuit**. Wider fractional operands and extra
polynomial products often dominate. Omitting low partial products
under an explicit phase-error budget is more effective.

## Reproducibility

Numerical polynomial and quantization check:

    python -m pytest -q tests/test_priority1_model.py
    python -m src.priority1_cli screen --family chebyshev --degree 2 --segment-bits 7 --precision 26

Build the conservative truncation circuit with native kickmix:

    python -m src.priority1_cli build --truncated --family chebyshev \
      --degree 2 --segment-bits 7 --precision 26 \
      --truncation-budget 0.12 --output out/cheb7_cut.kmb

Check all 3,133,452 inputs, each with 64 simulated measurement outcomes:

    python -m src.priority1_cli verify --circuit out/cheb7_cut.kmb --full

The exact finite-domain truncation experiment adds:

    python -m pytest -q tests/test_priority1_exact.py
    python -m src.priority1_cli build --truncated --family chebyshev \
      --degree 2 --segment-bits 7 --precision 26 --drops 10,8,15 \
      --output out/cheb7_exact.kmb

For an **official-style Q*C** measurement, run the competition verifier,
or the integrity-checked public mirror used by the GitHub workflow.
Full verified results for new variants should be inserted in this report
only after the corresponding workflow has returned PASS and score.

Native experiment logs:
- [Polynomial comparisons](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42625864096)
- [Truncated arithmetic](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42626861981)
- [Larger bounded truncation](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42627819141)
- [Exhaustive-cut candidates](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42628635243)
- [Full Q*C scoring of tuned 6/25](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42627447581)

## Branch note

The connector can create a new ref but cannot invoke GitHub's branch
rename endpoint or delete old refs. The branch `benchmark-1` was created
at the former branch's latest commit, preserving all history; the old
branch remains as an unmodified legacy ref. A true rename can be done
through GitHub repository Settings > Branches (or the branch UI) if desired.

## Strongest exact finite-domain truncation candidates (native kickmix)

To improve on the conservative uniform bound, the integer model
`src/priority1_exact.py` explicitly computes omitted low partial
products, and checks **every x in 1..3,133,452**. The model uses exact
integer multiplication-and-shift structure and double-precision relative
phase comparison; the .kmb is also tested with native kickmix on thousands
of boundary and pseudorandom inputs, 64 shots each.

| Cut settings (linear, square, quadratic) | Q | Toffoli-like | Static fine rotations | Full-domain worst error/tolerance | Native sampled failures |
|---|---:|---:|---:|---:|---:|
| (10, 8, 15) | 149 | 3,125 | 88 | 0.739087 | 0 |
| (10, 11, 12) | 149 | 3,009 | 88 | 0.717347 | 0 |
| **(10, 12, 11)** | **149** | **2,961** | **88** | **0.717347** | **0** |

This is the strongest generated direct-phase circuit in the current
experiment. [Native kickmix logs](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42628635243).
The full 200,540,928-trajectory scoring job is separately tracked at
[best candidate verifier](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42629669231).
Until that full run returns PASS, the table gives *native resource counts*,
not a final verified score.

### Why 30x is still a different research problem

The desired Benchmark I cost is 2,076,880 / 30 ≈ **69,229**.
The best current circuit has a **hard cost floor**
`Q * Toffoli_like = 149 * 2961 = 441,189`, even if all
arbitrary rotations were free. This is over six times above the
30×-improvement target. At Q=149, a 30× improvement requires fewer than
465 Toffoli-equivalent operations, including every synthesized rotation;
that cannot arise by coefficient tuning or modest product truncation.

The next high-risk research direction must therefore fundamentally
replace the multiplier/QROM work with a much cheaper diagonal-phase
representation. The present phase-arithmetic architecture has been
optimized materially, but does **not** meet the 30× objective.

## First completed full-scale score in this round

The conservative, tuned Chebyshev 6/25 truncated direct-phase circuit passed
the independent pinned public mirror of the official Benchmark I verifier for
all **200,540,928** trajectories. Its rotation-inclusive cost is
**667,314.430894474** (= 156 qubits × 4277.656608297911 Toffoli-equivalent
operations), an improvement over our earlier 876,513.3341652812 candidate.
Worst error/tolerance is 0.8682547320913627.

- [Full verified output](priority1_full_verified_score_6_25.json)
- [CI evidence](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42627447581)

The newer Chebyshev 7/26 with exact drops (10,12,11) has lower native
resource counts but its full-scale rotation-inclusive score is to be
reported separately, only after its own full verifier completes.


## Coefficient streaming to reduce peak ancilla count

A separate native kickmix candidate was implemented in
`src/priority1_streamed.py`. The normalization, phase model, and exact
truncation (10,12,11) are unchanged, but coefficients are loaded and
uncomputed separately for the constant, linear and quadratic phases.
This reduces the live workspace during multiplication.

| Coefficient strategy | Q | Toffoli-like | Static fine rotations | Sampled failures |
|---|---:|---:|---:|---:|
| Joint QROM, all coefficients live | 149 | 2961 | 88 | 0 |
| One coefficient at a time | **111** | 4063 | 88 | 0 |

[Full native kickmix test evidence](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42630388609).

The Q gain is large (38 qubits); the Toffoli penalty is also large
(1102 extra operations). Neither implementation is declared winner
until the rotation-inclusive resource scoring comparison is complete.

## Grouped coefficient lookup and barrel normalization

Four real native kickmix architectures were compiled and tested over
deterministic boundary/pseudorandom inputs with 64 shots each. All passed
the sampled oracle/scratch/phase verifier. All implement the same p7/q26
Chebyshev quadratic with exact dropped partial products (10,12,11).

| QROM coefficient grouping | Normalization | Q | Toffoli-like | Static fine rotations | Sample failures |
|---|---|---:|---:|---:|---:|
| A,B then C | Serial | 138 | 3512 | 88 | 0 |
| A then B,C | Serial | 122 | 3512 | 88 | 0 |
| A,B then C | Barrel | 138 | 3190 | 88 | 0 |
| **A then B,C** | **Barrel** | **122** | **3190** | **88** | **0** |

[Native sampled CI evidence](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42631502737).
[Full competition-style verification](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/workflows/benchmark1-grouped-barrel-full.yml)
has been launched for A then B,C / barrel. This is the strongest
work-space tradeoff so far, but does NOT yet have a verified Q*C score.
Its Toffoli-only lower bound is 122×3190=389,180, still >5.6× the 30×
benchmark target of 69,229. Solving 30× requires a radically different
diagonal-phase circuit, beyond optimizing this arithmetic architecture.

We kept the legacy experiment branch unchanged; the connector cannot
actually rename/delete a remote GitHub ref. `benchmark-1` is the new
active branch with the complete inherited history.
