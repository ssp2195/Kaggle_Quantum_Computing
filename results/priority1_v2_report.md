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
