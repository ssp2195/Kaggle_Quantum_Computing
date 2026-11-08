# Direct-phase prescreen, 2026-11-08

Greedy dyadic intervals cover **every integer** 1 through 3,133,452. A separate Chebyshev interpolant is fitted on each block and each integer point checked against a 75%-of-allowed phase error margin. This is a numerical screen, not a kickmix execution.

| Degree | Blocks | Worst error / tolerance | Naive parity rotations | Rotations after conservative pruning |
|---:|---:|---:|---:|---:|
| 1 | 11,773 | 0.749390 | 197,476 | 197,476 |
| 2 | 941 | 0.746976 | 239,573 | 198,287 |
| 3 | 269 | 0.727630 | 729,919 | 358,006 |

The angles are in turns. Each Boolean monomial is implemented as XOR parity phases. The coefficient pruning scheme drops coefficients only when their absolute sum is below 5% of the strictest phase-error tolerance within the block. The numbers are *unfused gadget estimates*, not measured instructions or official scores. Google's reference has 37 fine Z_POW rotations but also uses arithmetic; full circuit-resource comparisons require native kickmix and the Kaggle scorer.

**Finding:** the naive direct-phase architecture is mathematically feasible but likely uncompetitive due to high rotation counts. Native kickmix was unavailable in the local execution environment, so no .kmb was generated or fully verified. Research next: compress controlled-angle representation or hybridize with arithmetic rather than expanding high-degree phases into individual parity gadgets.
