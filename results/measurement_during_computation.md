# Benchmark I — Measurement *during computation*: feasibility and cost audit

**Status: native phase-instrument tests PASS; no demonstrated improvement beyond
existing MBUC for fixed Boolean predicates.** Resource estimates for a
full streaming phase oracle remain mathematical **proxies only**.

**Reproduction:** [native 2- and 3-bit tests](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42671400111).
Implementation: `src/benchmark1_measure_compute.py`;
tests: `tests/test_benchmark1_measure_compute.py`.

## 1. The useful measurement-based phase-computation identity

Suppose a temporary holds a Boolean function of the *still-coherent* input:

```
|psi>|0> -> sum_x alpha_x |x>|g(x)>
```

Apply Z_POW(t,e), where Z_POW on |1> contributes `exp(i*pi*e)`,
and measure t in the X basis, outcome m in {0,1}. The unnormalized
logical output is

```
sum_x alpha_x exp(i*pi*e*g(x)) (-1)**(m*g(x)) |x> / sqrt(2)
```

Both outcomes have probability 1/2 for **every** logical input state.
The factor `(-1)**(m*g(x))` is the only byproduct.

For g=a*b with live a,b, use conditional CZ(a,b) to cancel it at
Clifford cost. In native kickmix this is `hmr(t,m)` followed by
`ccz(a,b,m)`. For g=(a*b)*c, keep the helper t=a*b until the measured
outer result is corrected by conditional CZ(t,c); then erase t using
its own MBUC. This is measurement **integrated with phase production**,
rather than a request to preserve an output value.

**Important:** A Z-basis measurement of a data-dependent output
instead **dephases** the input across the two g=0/1 sectors. No
classical feedback can restore the lost off-diagonal coherences
without the quantum information that was disclosed by the measurement.

If g is cubic in original x bits and the intermediate AND is gone,
the random phase correction generally requires a conditional CCZ or
some newly computed helper: no free Clifford-only feedback is available.
This is a concrete special case of the causal correction criterion.

## 2. Actual native kickmix circuit comparison

A complete 2-input quadratic and a staged 3-input cubic
phase instrument were built with true native C++ kickmix, exercised
for **every computational input** and **64 independent measurement
trajectories per input**, with input restoration and coherent phase
checks.

| Boolean phase predicate | Mode | Q | Toffoli-like | Fine rotations | Checked trajectories | Failures |
|---|---|---:|---:|---:|---:|---:|
| a*b | Coherent reverse | 3 | 2 | 1 | 256 | 0 |
| a*b | Native compute-phase-del_and | 3 | 1 | 1 | 256 | 0 |
| a*b | Inline X measurement and conditional CZ | 3 | 1 | 1 | 256 | 0 |
| (a*b)*c | Coherent reverse | 5 | 4 | 1 | 512 | 0 |
| (a*b)*c | Native compute-phase-del_and | 5 | 2 | 1 | 512 | 0 |
| (a*b)*c | Inline X measurement and conditional CZ | 5 | 2 | 1 | 512 | 0 |

**Measured finding:** Moving the X-basis measurement into the phase
application preserves the correct oracle, **but costs exactly the
same** as ordinary kickmix MBUC of the computed phase temporary.
The 50% saving relative to coherent inverse is already captured by
the existing MBUC architecture. It cannot be claimed a second time.

This is consistent with Gidney's temporary AND: 4 T gates to compute,
0 T gates for its measurement-based erasure (Quantum 2, 74 (2018):
https://quantum-journal.org/papers/q-2018-06-18-74/).
More generally, Clifford operations, stabilizer states and Pauli
measurements do not generate the required nonlinear/non-Clifford
resource for free. Teleportation moves this resource into ancilla
state preparation, which must still be charged.

## 3. Could measured phase predicates avoid the arithmetic entirely?

The existing grouped p7/q26 quadratic approximation contains
coefficients of widths (27,18,11), a 13-bit residual mantissa, and
truncation parameters (linear,square,quadratic)=(10,12,11).

If we distribute the polynomial into independent controlled Boolean
phase terms instead of explicitly assembling the integer products,
the *raw unoptimized source monomial* count is:

| Term family | Count |
|---|---:|
| Linear coefficient-bit × mantissa-bit terms | **179** |
| Unique quadratic coefficient-bit × mantissa-bit-pair terms | **202** |
| Of the latter: cubic cross-bit terms | 161 |
| Of the latter: squared-bit degree-2 terms | 41 |
| Total raw nonlinear phase terms | **381** |

These counts are for a *directly distributed source polynomial*;
they **do not** prove exact equality to the original nested
partial-product truncation (carry and column-cut behavior differs).
They also do not prove a lower bound: terms can be shared,
combined, omitted, or represented by a different phase network,
subject to complete-domain error verification.

**Cost proxy.** The actual best full-verified circuit has Q=112,
3218 Toffoli-like, 88 fine rotations and
C=3841.603737748773. After accounting for the one T-like gate,
the observed rotation component is approximately **7.083565**
Toffoli-equivalent units per fine rotation. Using that
**candidate-dependent observed rate only as a comparison model**,
381 separately synthesized angles would cost the following
*before* predicate computation, coefficient lookup, normalization,
workspace, or conditional shift selection:

| Hypothetical peak Q | Rotation-only Q*C at 381 angles |
|---:|---:|
| 22 (absolute input-qubit minimum) | 59,374 |
| 40 | 107,954 |
| 50 | 134,942 |
| 63 | 170,027 |
| 112 | 302,270 |

30x requires **Q*C ≤69,229**. At Q≥40 the naive one-rotation-per-term
approach already exceeds the target from synthesis alone under
the current effective rate. At the impossible zero-workspace
Q=22 limit, even the naive unshared forward AND operations and
coefficient selection would consume the remaining budget.
Moreover, the 4-bit exponent k makes each direct angle
**variable**; implementing each with the existing `_phase`
barrel-shift routine can require many additional rotations and gates.

Code for these numbers: `src/benchmark1_measure_streaming_cost.py`.
A separate native one-term experiment is included in the CI; consult
its measured resource output before claiming a variable-exponent score.

## 4. Feasibility of other measurement-during-computation methods

**(a) Pauli/Clifford measurement-only AND computation**:
Not a viable deterministic coherent replacement for non-Clifford
Boolean arithmetic without additional nonstabilizer resource states.

**(b) Magic-state-injected / gate-teleported Toffoli**:
Potentially valid but resource preparation and adaptive corrections
must be included. In kickmix, a forward Toffoli-like gate already costs
one Toffoli equivalent; replacing it with about four T gates offers
no inherent score saving, and can increase peak qubits.

**(c) Digit-by-digit measured carries, used as classical controls**:
If the measured carries depend on the coherent data, a Z-basis
measurement leaks input information. A teleportation or masked
measurement scheme is required; there is no free transition to a
classical carry in the coherent phase oracle.

**(d) Phase-only X-basis measurement of low-degree predicates**:
Correct and supported; it removes the need to retain the *terminal*
function bit. It is most promising if a new compact representation
uses a small set of **shared, degree-2-in-live-values** predicates
without generic full-width addition or QROM. This is a different
architecture from monomial-by-monomial streaming and remains
unimplemented at full width.

**(e) Repeat-until-success or adaptive angle injection**:
Potentially useful for reducing *rotation synthesis* in another
resource model, but requires deterministic correction, appropriate
bounded branching and accounting for consumed T/magic states and
ancillas. We have not built a verifier-accepted RUS architecture here;
it must not be counted as an improvement.

## 5. Bound from the current best fully verified circuit

With current native nonlinear gate inventory unchanged,
**Q × 3218 ≥22 × 3218 = 70,796**, already greater than the
**~69,229** 30x target *even if rotations were free and there were
no ancillas*. At the measured Q=112, 3218 forward Toffoli-like gates
alone cost at least **360,416** before synthesized fine rotations.

So any strategy that **only** measures more workspace while preserving
the present arithmetic network cannot achieve 30x.

## Decision

**Do not replace the verified 112-qubit best with inline-measured phase
gadgets:** the native cost is identical to existing measurement-based
uncomputation for the same predicate, and raw streaming adds too
many rotations. Retain the MBUC early-normalization improvement.

A meaningful next experiment requires **joint Boolean DAG and phase-tap
optimization**: identify shared live factors so only a small number
of X-measured quadratic phase predicates remain, with causal Clifford
corrections, plus tight angle-synthesis and qubit-lifetime scheduling.
Do not confuse the design hypothesis with a validated score.

No new full-domain circuit or Kaggle submission is claimed here.
