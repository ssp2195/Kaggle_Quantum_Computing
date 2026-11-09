# Benchmark I architectural experiments — 2026-11-09

## Decision

**Retain the grouped Chebyshev quadratic phase circuit with barrel normalization.**
Among the independently compiled and verified architectures in this repository,
it remains the lowest-cost full-domain PASS. The 30x target is NOT achieved.

### Best independently verified result

- Architecture: quantized Chebyshev quadratic, segment bits 7, coefficient
  precision 26, exact drops (10,12,11); coefficient grouping A then BC;
  barrel normalization.
- Real kickmix circuit: 122 qubits, 3,190 Toffoli-like operations, 88 fine
  Z_POW rotations.
- Mean Toffoli-equivalent (including rotation synthesis): 3,813.603737748773.
- Benchmark I Q*C: **465,259.6560053503**, improvement **about 4.46x**
  relative to published baseline approximately 2,076,880.
- Full pinned public competition-verifier mirror: **PASS** on **200,540,928**
  trajectories (3,133,452 inputs x 64 shots).
- Evidence: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42632214980
- Rebuild: python -m src.benchmark1_best --output out/best_invsqrt.kmb
- This is NOT a new Kaggle leaderboard submission score.

## Full-domain mathematical feasibility

A normalized Chebyshev piecewise-linear phase with p=10 segment bits,
q=28 fractional coefficient bits, and 7 omitted low columns in the
linear product passed the exhaustive model screen:
**3,133,452 integer inputs; 0 failures; maximum relative phase error
0.8008876778 of the allowed tolerance at x=1,056,255.**
There are 3072 nontrivial intervals and 4096 table rows.
This mathematical result DOES NOT certify the generated quantum circuit.

## Actual native kickmix circuit measurements

| Architecture | Peak qubits | Native Toffoli-like | Fine rotations | Sample correctness | Q*Toff lower bound | Decision |
|---|---:|---:|---:|---|---:|---|
| Best grouped quadratic A/BC barrel | 122 | 3190 | 88 | FULL PASS | 389180 | Keep |
| Quadratic AB/C barrel, original squaring | 138 | 3190 | 88 | PASS | 440220 | Worse Q*C than best |
| Quadratic AB/C barrel, symmetric squaring | 138 | 3892 | 88 | PASS | 537096 | Reject |
| Linear p10/q28, joint QROM, XOR cleanup | 135 | 9308 | 65 | FAIL | 1256580 | Reject |
| Linear p10/q28, streamed QROM, XOR cleanup | 107 | 17496 | 65 | FAIL | 1872072 | Reject |

References:
- Native architecture CI: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42648568168
- Native stage-cost profiler: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42647183393
- Stage and numerical reports: results/architecture_stage_profile.json
  and results/architecture_linear_screen.json

The naive symmetric-square arithmetic identity is exact; implementing
individual cross-products using controlled constant adders doubled the
quadratic-stage operation cost from 814 to 1516 native Toffoli-like gates.
Therefore do not use this circuit in a submission.

For linear circuits, the first native compilation failed during the
measurement-based QROM erasure because extra clean workspace was needed.
Replacing the unlookup with an XOR self-inverse QROM pass made compilation
possible but yielded INVALID sampled circuits. The native phase and
scratch-error ratios were far outside tolerance. Independently of this
bug, Q*Toffoli is already far above the current verified Q*C. These are
rejects, not claims of valid alternative implementations.

## Exact measured bottleneck

| Component | Toffoli-like operations |
|---|---:|
| Normalize plus restore | 228 |
| AB QROM plus cleanup | 551 |
| C QROM plus cleanup | 551 |
| Constant phase network | 180 |
| Linear reversible product and phase, including cleanup | 866 |
| Quadratic reversible product and phase, including cleanup | 814 |
| TOTAL | 3190 |

Isolated-component peak qubit counts are not additive. Native stats
were obtained by building kickmix circuits, not inferred from arithmetic
gate formulas.

The coefficient QROMs together cost 1102 Toffoli-like gates, and
the linear/quadratic product stages another 1680. At 122 qubits, a
30x target cost approximately 69229 permits fewer than 568 total
Toffoli-equivalent operations including synthesized fine rotations.
Therefore coefficient tuning inside this architecture cannot achieve 30x.

## Independent direct-diagonal numerical screens

On every one of the 2^20 normalized mantissa inputs:

- Walsh parity truncation required around 1800 retained independent
  phase terms to meet the target relative numerical accuracy.
- Boolean algebraic-normal-form truncation required around 2000
  retained monomial terms to pass.

These are *mathematical full-grid screens*, not native kickmix circuits,
not formal lower bounds, and not rejection of possible SHARED
nonlinear predicate networks. Independent phase gadgets at these
sizes are poor candidates for the 30x resource target.

Source: src/priority1_walsh_screen.py, experiments/anf_phase_screen.py.
Results: results/priority1_walsh_feasibility.md,
results/architecture_anf_screen.json.

## Architecture selected for next major construction: XAG + in-place Horner

Public independent analyses of approximately 63-qubit, 20x-class
inverse-square-root circuits show a fundamentally more compact approach:

1. An 8-bit XOR/AND graph (XAG) coefficient lookup instead of
   generic kickmix binary-address lookup.
2. Measured temporary ANDs, corrected coherently, to erase lookup
   dependencies without repeating all Toffoli computation.
3. In-place signed-digit centered quadratic Horner updates rather than
   storing wide full product registers.
4. Exponent/shift metadata parked in zeros of the original input,
   narrow barrel shifting, and two-pass phase taps.

A reported 63-qubit, ~102664-cost PUBLIC implementation uses around
172 ANDs for the primary 8-bit lookup, ~41 fine rotations, and a
small workspace peak. This is **external research**, not an implementation,
replication, or verified result of this repository.

Sources, methodology, and reproducibility cautions:
- https://github.com/moodyastra/quantum/blob/0458f0c929f13a28d83b46ecafa2b31ebd4bfded/docs/analysis/invsqrt.md
- https://github.com/moodyastra/quantum/blob/0458f0c929f13a28d83b46ecafa2b31ebd4bfded/docs/analysis/invsqrt_critique.md
- https://github.com/clduab11/quantum-hackathon/blob/b92db3adf01a5b102aa85c2975f16e402e50cbd3/docs/research/PRA-320-invsqrt-source-audit.md

These analyses warn that some public binaries depend on frozen
coefficient/synthesis data and are not yet reproducible from the
distributed generators. Do not use another competitor's unpublished
search products or assert independent construction without attribution.
We have NOT imported any external .kmb or its vendor data.

This architecture is the best RESEARCH DIRECTION identified from the
cost anatomy and public comparable examples, but it has **not yet been
built in this repository**. The actual best VERIFIED implementation
here remains grouped Chebyshev 7/26 at cost 465259.66.

A reproducible original XAG/Horner implementation must be evaluated
with the real kickmix C++ backend and full official-style verifier
before the 30x goal can be claimed.

## Cost admission gate

Run: python -m src.benchmark1_cost_gate

The gate rejects a candidate from full-cost submission experiments
if Q times its native Toffoli-like count is above the current best
verified Q*C, or it failed sampled native coherence. Passing this
necessary condition does not imply actual eligibility for the
30x goal or mathematical correctness.

## Open target

Reference Q*C approximately 2076880 / 30 = **69229.33**.
Current verified own best 465259.66: requires a further
**6.72x improvement**. No tested native architecture in this round
meets or approaches 69229. The user's last actual combined
Kaggle leaderboard score was 1.33 from the previous submission;
this experiment has not submitted new circuits to Kaggle.
