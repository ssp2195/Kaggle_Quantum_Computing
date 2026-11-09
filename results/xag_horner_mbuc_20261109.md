# XAG / in-place Horner and measurement-based cleanup — 2026-11-09

## Scope and evidence labels

All candidate resource numbers below are from actual **native C++20 kickmix**
`.kmb` circuits; mathematical results are exhaustive over the full
**3,133,452 integer input** domain; quantum circuit correctness is from
native sampled coherent measurements with 64 trajectories per input.
The best prior fully verified circuit has
**Q=122, C=3813.603737748773, Q*C=465259.6560053503**, passing
**200,540,928** trajectories. The requested 30x target is **Q*C <= 69229**.

Full competition-style verification of the new compact-normalization circuit
was launched at
https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42656330676.
Do not report its final Q*C until that verifier completes.

## A. In-place Horner phase

Normalized polynomial:
`P(u)=A-B*u+C*u*u=A-(B-C*u)*u`.

After quantizing integer coefficients, use a small guard `g`:
```
R = (B << g) - CUT(C * mantissa, low_bits - g)
Y = (A << g) - CUT(R * mantissa, low_bits)
Z_POW from Y with k-dependent binary phase shift
undo Y, undo R, erase coefficients, restore x
```

Each CUT excludes particular low bit-products **before** addition, preserving
the exact binary arithmetic implemented by the circuit. In-place subtraction
must propagate borrows into **all higher destination bits**; truncating
a borrow at a row boundary caused our first native smoke tests to fail.
The corrected code uses `dest[start:]`. This bug and fix were independently
caught by real kickmix phase simulation.

### Measured Horner results

| Candidate | Q | Native Toffoli-like | Fine rotations | Exact math worst/tolerance | Native sample |
|---|---:|---:|---:|---:|---|
| p7 q26, g3, AND MBUC | 149 | 3351 | 36 | 0.818315 | PASS |
| p7 q26, g3, fully coherent AND erase | 149 | 3891 | 36 | 0.818315 | PASS |
| p7 q26, g4, AND MBUC | 153 | 3451 | 37 | 0.686444 | PASS |
| p7 q26, g4, coherent erase | 153 | 4037 | 37 | 0.686444 | PASS |
| p6 q26, g5, 8-bit ANF-XAG + MBUC | 159 | 3682 | 38 | 0.784008 | PASS |
| p6 q26, g5, built-in QROM + MBUC | 159 | 3467 | 38 | 0.784008 | PASS |

Evidence:
- https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42653214903
- https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42653946003

MBUC saves **540** Toffoli-like gates at g3 and **586** at g4, while
maintaining coherent phase and restoring scratch. Despite this substantial
savings, the best g3 variant has a hard native gate-cost floor
`149 * 3351 = 499299`, exceeding the prior fully verified Q*C even
**before paying for synthesized rotations**.

## B. Centered signed-digit Horner

Use `v = u - 1/2`, represented with no extra input register by flipping
the top mantissa bit to create a two's-complement signed offset.
Fit `P(v) = A-B*v+C*v*v` directly. For its sign bit, use a single
controlled reverse-direction add rather than computing a full AND row.

| Candidate | Q | Native Toffoli-like | Fine rotations | Math worst/tolerance | Sample |
|---|---:|---:|---:|---:|---|
| p6 q26, g4, QROM | 154 | 3409 | 37 | 0.775636 | PASS |
| p6 q26, g4, XAG | 154 | 3624 | 37 | 0.775636 | PASS |
| p7 q26, g3, QROM | 146 | 3347 | 36 | 0.789843 | PASS |

Evidence:
https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42655007647

The centered form needs fewer guard bits on some interval choices, but
does **not** yet reduce total Q*C below the best grouped quadratic circuit.

## C. Early measurement-corrected erasure of normalization flags

The old barrel normalizer keeps 10 prefix flags alive through the complete
phase calculation. The compact normalizer:

1. Computes the prefix ANDs on the original input;
2. Encodes exponent k into a four-qubit register;
3. Applies native `del_and` with conditional CZ correction *before* the
   controlled barrel rotation, freeing all 10 prefix qubits;
4. Rotates the input and evaluates the diagonal phase;
5. Unrotates the input, **recomputes** prefix flags, and subtracts k;
6. Erases the recomputed flags, leaving all scratch clean.

This is a causal correction scheme: the original AND controls are
available while each temporary is measured. Using a measurement
before the input has been restored would not allow the same
unconditional correctness argument.

Measured complete 22-bit grouped quadratic, `p7/q26` with exact partial
product drops (10,12,11):

| Normalization | Q | Toffoli-like | Fine rotations | Native sample |
|---|---:|---:|---:|---|
| Existing barrel | 122 | 3190 | 88 | Full verifier PASS |
| Compact measurement-corrected barrel | **112** | **3218** | **88** | PASS on **553024 trajectories** |

The compact scheme changes the isolated normalization from 228 to
256 Toffoli-like operations, while reducing the full-oracle peak Q by
10. The resource floor is 112*3218=360416.
The native smoke test passed; the full 200540928-trajectory verifier
was launched separately.

Evidence:
https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42655591929

## D. What measurement helps and what it cannot do

For a temporary Boolean AND t=a*b, X-basis measurement gives outcome
mu and a phase `(-1)**(mu*a*b)`. The conditional CZ(a,b) exactly cancels
this for every outcome. The `del_and` primitive implements this correction.

For a cubic temporary w=(a*b)*c, if t=a*b remains live, measure w and
correct CZ(t,c), then measure t and correct CZ(a,b). These are Clifford
fixups. The staged toy circuit is in
`src/benchmark1_mbuc_feasibility.py`.
It is *not* evidence that arbitrary cubic byproducts can be removed
without a helper. Without t live, the w byproduct requires a conditional
cubic phase on a,b,c, typically CCZ, which is non-Clifford.

Crucially, measuring the **Horner result accumulator** is not a
drop-in replacement for reversing arithmetic. The measured output bits
are generally high-degree Boolean functions of x. Their byproducts
cannot in general be corrected by only Z/CZ gates acting on the original
inputs. A correction decomposition must be derived, or a helper retained.
Also, a coherent AND cannot be computed for free using only stabilizer
operations: forward CCX or equivalent non-Clifford resource remains.

This is consistent with the measured >500-Toffoli savings for MBUC
compared with coherent AND erasure, **but does not imply zero-cost
arithmetic or an automatic 30x solution**.

## E. Decision and next architectural work

Among all fully end-to-end generated candidate circuits, the 8-bit
naive ANF-XAG and in-place Horner approach has not beaten the original
grouped quadratic design. The best structural improvement in this
round is the **compact early-MBUC normalizer**, pending full scoring.

For a credible 30x route, remaining high-impact changes are:

- Exact **multi-output** XAG synthesis for coefficient bit planes,
  rather than all monomial ANDs. Goal below about 170–200 ANDs for
  an 8-bit coefficient lookup, with preserved attribution/reproducibility.
- A narrow bit-sliced phase accumulator that emits phase taps without
  loading and retaining separate wide A/B/C coefficients.
- Reuse vacated input bits and prove peak live-qubit schedules, rather
  than independently allocating 60+ coefficient and temporary qubits.
- Stage-local measurement cleanup only when the resulting Boolean
  byproduct is correctable by available Clifford phases.

An independent publicly reported 63-qubit/~20x architecture with
optimized multi-output XAG and in-place rows is analyzed at:
https://github.com/moodyastra/quantum/blob/0458f0c929f13a28d83b46ecafa2b31ebd4bfded/docs/analysis/invsqrt.md
It is *external evidence*, **not** an original output of this repository,
and may depend on frozen synthesis data. No external circuit binaries
were imported.

**No native candidate built in this phase has met the 30x target.**

## F. Further native compact-Horner and causal measurement results

The compact normalizer also integrates with in-place Horner:

| Guard bits | Peak qubits | Toffoli-like | Fine rotations | Status |
|---|---:|---:|---:|---|
| 3 | 139 | 3379 | 36 | Native sampled PASS |
| 4 | 143 | 3479 | 37 | Native sampled PASS |

Evidence: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42655937489

Both are lower-workspace than the earlier Horner variants, but their
native Toffoli cost floors of 469681 and 497497, respectively, already
exceed our prior 465260 full cost when rotations are assumed free.
Therefore neither warrants a new full verifier run.

The independent native staged cubic-AND toy passed every 3-bit input
under 64 trajectories, preserving phase and cleaning scratch. At
**five qubits** it costs **2 Toffoli-like gates with MBUC**, compared
with **4 with coherent CCX erasure**; all forward CCX gates remain.

Evidence: https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42656187092

**Consequence:** measurement erasure can halve the gate count of an AND
compute/erase pair and can release high-degree temporary computations
when intermediate factors remain live, but it cannot by itself make
reversible Horner arithmetic or QROM free. The compact-normalization
schedule is the strongest new measured application to this task.

## G. Strong feasibility bound for 30x

Even the **mandatory 22 input qubits** multiplied by the current
compact-normalizer's **3218 native Toffoli-like operations** gives
`22 * 3218 = 70796`, which is greater than the approximately
`69229` 30x target — **before any synthesized rotation cost**.
Therefore **workspace optimization alone cannot reach 30x** with
our present gate inventory, even in the physically impossible ideal
case of zero ancillas.

At the actually measured 112 qubits, the allowed total
Toffoli-equivalent cost would be only `69229/112 ≈ 618`.
Our compact-normalizer already has **3218** non-Clifford
Toffoli-like gates, and its 88 fine rotations add substantial synthesis
cost. Causal measurement cleanup has already been used on the safe
temporary ANDs, so simply repeating MBUC cannot close that gap.

A 63-qubit circuit has a total allowance about `1099` eq.
For comparison, the externally reported 63-qubit / ~20x design uses
roughly 1334 mean Toffolis *plus* ~41 fine rotations. Even starting
from that stronger public architecture, 30x requires further
reduction in nonlinear gate operations, workspace, or both.

This arithmetic bound does not rule out another architecture.
It rules out expecting the requested 30x solely from cleaning
the current polynomial/oracle workspace with measurements.

## H. Upstream kickmix adder already uses measurement cleanup

A source audit of upstream native kickmix's
[reversible quantum adder](https://github.com/quantumlib/Qualtran/blob/main/kickmix/src/kickmix/gen/adders/gen_iadd.cc)
shows that `gen_iadd` already includes X-basis measurements
(`loop.hmr`) and conditional Clifford correction (`loop.cz_if`),
and in certain carry-cleanup paths calls `builder.del_and`.
The controlled increment generator
[`gen_iadd1.cc`](https://github.com/quantumlib/Qualtran/blob/main/kickmix/src/kickmix/gen/adders/gen_iadd1.cc)
uses a similar measurement-corrected cleanup chain.

This is **source-derived evidence**, distinct from our measured native
circuit metrics. Our MB-versus-coherent Horner comparison deliberately
varies only the separate partial-product AND row cleanup; both
implementations also inherit kickmix's existing measurement-aided
adder cleanup. It would therefore be incorrect to claim that replacing
`iadd` with a generic MBUC adder could automatically halve all its
Toffoli-like costs again. Only a validated different carry/dataflow
architecture could save further operations.
