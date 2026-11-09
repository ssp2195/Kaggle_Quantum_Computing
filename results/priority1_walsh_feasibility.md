# Sparse Walsh direct-phase feasibility study

This is a **full 2^20 normalized-grid** experiment, not a random sampling or
tiny circuit toy. The normalized input ranges from `1` to `2` with a
20-bit fractional mantissa. We expand `f(z)=1/sqrt(z)` in the Walsh parity
basis, keep the largest-magnitude coefficients, and inverse-transform to
evaluate the **worst weighted relative phase error** on every grid input.

Walsh functions correspond to input parities, which are implementable using
Clifford CNOT networks and `Z_POW` rotations in kickmix. A full Benchmark I
circuit would also need normalization, exponent scaling and conditional
selection; those costs are **not** included here. Thus the screen is a
favorable test of naive sparse parity-phase synthesis, not a proven no-go
theorem for all diagonal representations.

| Top retained Walsh terms | Worst error / tolerance | Pass? |
|---:|---:|---|
| 500 | ~57.285 | No |
| 1000 | ~6.782 | No |
| 1500 | ~1.730 | No |
| 1700 | ~1.118 | No |
| 1800 | ~0.905 | Yes |
| 2000 | ~0.622 | Yes |
| 3000 | ~0.135 | Yes |

The 1,800-rotation candidate is *mathematically* valid on the represented
mantissa interval but remains far too costly to justify native full-scale
circuit generation as the next step. These rotations are not free under
the competition's synthesis metric. Even ignoring normalization and
control costs, 1,800 fine angles and the mandatory 22-qubit register
strongly disfavor the required ~69,229 resource-cost target.

Reproduce:

    python -m pytest -q tests/test_priority1_walsh_screen.py
    python -m src.priority1_walsh_screen --bits 20 --output runs/walsh_screen.json

**Decision:** do not proceed with naive largest-coefficient Walsh truncation.
A cheaper phase architecture would need shared nonlinear structure,
multi-resolution representations, or block-specific phase synthesis, rather
than one independently synthesized rotation per retained parity.
