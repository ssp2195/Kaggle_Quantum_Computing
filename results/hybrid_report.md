# Benchmark I: normalized hybrid quadratic phase (2026-11-08)

**Architecture:** In-place reversible 22-bit normalization by leading-zero-pair detection, QROM lookup of three nonnegative coefficients, and direct phase application of the independent constant, linear, quadratic contributions. The inverse-square-root *result register is never constructed*. No Newton–Raphson iteration is used. Independent terms are uncomputed and the input x restored.

Let x=4^k z, 1<=z<4. Divide z into dyadic intervals of width 2^-p, write z=(j+u)/2^p with 0<=u<1, and apply

    2^-k [ A_j - B_j u + C_j u^2 ],

where A_j, B_j, C_j are rounded fixed-point integer coefficients of precision q, divided by 2^q. The actual kicked phase is this quantity in *turns*. We store the scale s=10-k in four reversible qubits and apply a barrel shift on each term register, instead of 11 different sets of controlled fine Z rotations.

## Full-domain mathematical screening (all 3,133,452 integers)

| p | q | intervals | max error / tolerance | failures |
|---:|---:|---:|---:|---:|
| 7 | 30 | 384 | 2.491798 | 66,663 |
| 8 | 24 | 768 | 2.475865 | 145,624 |
| 8 | 26 | 768 | 0.617764 | 0 |
| 8 | 28 | 768 | 0.329472 | 0 |
| 8 | 30 | 768 | 0.311815 | 0 |
| 9 | 26 | 1536 | 0.632347 | 0 |
| 9 | 28 | 1536 | 0.155363 | 0 |
| 9 | 30 | 1536 | 0.061488 | 0 |
| 10 | 28 | 3072 | 0.165788 | 0 |

This is exhaustive only over the *quantized mathematical expression*. It is not quantum simulation.

## Native C++20 kickmix measurements

Real circuits with exactly one named 22-qubit input register x were constructed, serialized to .kmb, loaded back, and run on 2,509 sampled inputs times 64 trajectories per input. The selected inputs include boundaries, powers-of-four and deterministic pseudorandom values.

| p | q | Q | Toffoli-like | T-like | fine Z_POW | instructions | failed trajectories | worst sampled error/tolerance |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 26 | 165 | 4490 | 1 | 119 | 72593 | 0 / 160,576 | 0.526944 |
| 8 | 28 | 177 | 4808 | 1 | 125 | 79585 | 0 / 160,576 | 0.320188 |
| 9 | 28 | 164 | 5356 | 1 | 119 | 126735 | 0 / 160,576 | 0.146373 |

Baseline reference: Q=260, Toffoli count=7,735, fine Z_POW=37 and official task score about 2.077 million. **Do not interpret a lower native Toffoli count or lower-bound Q*C as an official score**. The official scorer charges for all fine Z_POW rotations, and the hybrid circuits use more rotations.

**Native smoke CI:** https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42509696864

**Full trajectory matrix:** https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42509940993

Full-scale jobs are configured to test **3,133,452 x values × 64 shots = 200,540,928 trajectories each** using real kickmix. This document distinguishes pending full-scale results from already completed experiments. The official verifier requires authenticated Kaggle competition files and remains authoritative for a submitted PASS and Q*C score.

## Reproduce

    python -m pytest -q tests/test_hybrid_model.py
    python -m src.hybrid_cli screen --output runs/hybrid_numerical.json

    # Requires native kickmix installed from quantumlib/Qualtran/kickmix
    python -m src.hybrid_cli build --segment-bits 8 --precision 26 --output out/hybrid_8_26.kmb
    python -m src.hybrid_cli verify out/hybrid_8_26.kmb --random-count 2048
    python -m src.hybrid_cli verify out/hybrid_8_26.kmb --full --output runs/full.json

    # Requires separately downloaded official Kaggle competition files
    bash scripts/verify_official.sh out/hybrid_8_26.kmb official

The simulation is not a substitute for the official metric notebook or Kaggle private-salt verification. The nominal numerical margin must cover circuit implementation and fine rotation synthesis; resource scoring must include all rotation costs. Prefer p=8/q=26 for the lowest native cost among these three, pending full authoritative verification.


## Analytical error certificate independent of pointwise sampling

The exact Taylor remainder of `f(z)=z^-1/2` after quadratic order satisfies

    |R3| <= (5/16) * z0^-7/2 * h^3,    h=2^-p.

Dividing by `f(z)`, using `z0>=1` and `z<=z0+h`, gives the
**relative** error bound

    |R3|/f(z) <= (5/16) h^3 sqrt(1+h).

Quantizing each of the three Taylor coefficients to nearest multiples of
`2^-q` changes their summed value by at most `(3/2)2^-q` because `0<=u<1`.
Also `f(z)>=1/2`. Hence, normalized by the allowed relative tolerance
`2^-24`, a *uniform bound over the continuum of normalized inputs* is

    error/tolerance <= (5/16) * 2^(24-3p) * sqrt(1+2^-p)
                       + 3 * 2^(24-q).

For p=8, q=28 this is approximately 0.501 and therefore is strictly
below 1. The bound also holds for p=9, q=26 and p=9, q=28. For p=8,
q=26 it is conservative (above 1); the exhaustive finite-domain screen
instead shows an observed maximum ratio 0.618 on every admissible integer.

The certificate is for exact arithmetic on the quantized coefficients.
The actual kickmix circuit is designed to realize those exact integers
with reversible adders, carry-safe multiplication and exact fractional
Z_POW exponents. An official verifier run is still required to establish
that the generated .kmb meets the competition's complete implementation
and resource-accounting rules.
