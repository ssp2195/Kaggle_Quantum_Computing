# Measurement *during computation*: feasibility study

## Aim

Determine whether measuring temporary Boolean function values **as part of
generating the phase** can reduce Toffoli-equivalent cost below existing
compute--phase--measurement-erase circuits. The 30x Benchmark I target is
approximately 69,229 for Q*C, compared with our verified 112-qubit best
at Q*C 430,259.62.

## Exact direct phase instrument

Prepare `|x>|g(x)>`. Apply `Z_POW(t,e)` to the temporary, then
X-basis measure it with outcome `m`. The induced logical phase, before
feedback, is `e^{i*pi*e*g(x)}*(-1)^{m*g(x)}`. A conditional Clifford
phase can cancel the second factor when `g` is quadratic in live bits.

For `g(x)=a*b`, applying `CZ(a,b)` conditioned on the measurement
result deterministically implements the exact desired phase without
leaving `g` computed. The forward AND remains necessary.

For `g(x)=(a*b)*c`, keep `t=a*b` live, measure the outer AND and
correct `CZ(t,c)`, then clean `t`. This is a legal causal correction.
If `t` has been discarded, the byproduct is a cubic phase requiring
non-Clifford correction in general.

The experiment compares three native kickmix circuits for both
quadratic and cubic predicates:

- `coherent`: compute AND, apply phase, coherent CCX reversal.
- `mbuc`: compute AND, apply phase, native `del_and` with correction.
- `inline`: compute AND, apply phase, explicit `hmr` plus
  classical-outcome-controlled CZ. This measurement is part of the
  phase application.

All variants preserve input and workspace exactly. The CI checks every
input and 64 measurement trajectories per input, as well as the native
Toffoli/rotation/qubit costs.

## Fundamental limitations

**Measuring a data-dependent output bit in Z does not implement a phase
oracle**: it reveals information about x and decoheres superpositions
across the two measurement outcomes. Correcting branch phases cannot
restore the lost coherence.

**Measuring arbitrary Horner output bits in X is not free**, either. The
byproduct `(-1)^{m*f(x)}` can be higher-degree or otherwise outside
the span of currently available quadratic Clifford corrections. Cheap
measurement of those outputs requires retained live helpers or an
additional nonlinear correction circuit.

**Gates with measurement/classical feedforward are not automatically
non-Clifford for free**. A coherent nonlinear AND needs non-stabilizer
resources: e.g. Gidney's temporary logical AND costs four T gates
forward, although its measurement-assisted cleanup costs no T gates.
Native kickmix scores forward Toffoli-like computation as a nonzero
cost. Magic-state injection moves this cost into the resource state;
it cannot be silently excluded from the competition resource budget.

**Measurement-augmented coefficient or phase lookup** might be useful
only if replacing whole arithmetic blocks with a much smaller set of
quadratic-in-live-bits phase predicates. The previous full-grid Walsh
and ANF screens needed ~1800 and ~2000 independent terms, respectively.
That is not compatible with a 30x cost under straightforward native
rotation synthesis. Shared nonlinear features and compact phase taps
remain a different, unverified architecture.

### Literature

Craig Gidney, *Halving the cost of quantum addition*,
Quantum 2, 74 (2018),
https://quantum-journal.org/papers/q-2018-06-18-74/

Source of exact kickmix model: pinned Qualtran native extension, SHA
`9a0c11a4ee7697c082fbf8f2e1e8bf4a71fcb453`.

## Status

The above is theoretical analysis and a native test specification.
The run logs must be checked before claiming native measurement-versus
compute gains. Results will be appended after CI completes.
