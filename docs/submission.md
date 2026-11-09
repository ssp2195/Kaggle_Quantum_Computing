# Preparing the complete Kaggle submission

**Only Benchmark I (invsqrt) has been optimized.** Kaggle requires a **ZIP
containing all three benchmark tasks**, so keep the competition-provided
sample/reference DQI and DFTHCSOS files unchanged.

## Verified Benchmark I result

The publicly mirrored inverse-square-root verifier passed the complete
200,540,928 shots for the 8/26 circuit on 2026-11-09.

- Cost Q*C: **876,513.3341652812**
- Qubits: **165**
- Toffoli-equivalent cost: **5312.202025244128**
- Toffoli-like: **4490**
- Fine rotations actually executed per shot: **116**
- Worst phase-error/tolerance: **0.6177643626145476**
- Circuit SHA-256: `5ba3ac2adbf6d49def4b1ef2d00efcdbba73e82d94d178fa2613d1352e9ce26c`

[GitHub Actions full verifier](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42554578982).
This **is not** a Kaggle leaderboard result.

## Assemble the competition package

Download the complete competition data via the Kaggle website or an
authenticated Kaggle client after accepting the competition rules. You need
`sample_submission/dqi.kmb` and the entire
`sample_submission/dfthcsos/` directory, including `submission.npz`.

Get the validated 8/26 circuit from the
[mirrored-metric-8-26 workflow artifact](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42554578982),
or rebuild it with native kickmix:

    python -m src.hybrid_cli build --segment-bits 8 --precision 26 --output out/hybrid_8_26.kmb

Create the ZIP:

    python -m scripts.package_submission \
        --invsqrt out/hybrid_8_26.kmb \
        --reference official \
        --output out/submission.zip

Validate the ZIP contents:

    python -m zipfile -l out/submission.zip

Submit with an authenticated Kaggle CLI (configure authentication *locally*;
do not commit or echo tokens):

    kaggle competitions submit \
        -c google-verifiable-quantum-circuit-compilation \
        -f out/submission.zip \
        -m "Normalized quadratic direct-phase invsqrt 8/26; unchanged reference DQI/DFTHCSOS"

Check the actual Kaggle leaderboard result:

    kaggle competitions submissions -c google-verifiable-quantum-circuit-compilation

Remember the leaderboard score is a **geometric mean of the three task
improvement factors**, not the Benchmark I cost. Assuming the two reference
tasks are unchanged, the expected combined score is approximately
`(2076880 / 876513.3341652812) ** (1/3)`, which is roughly 1.33,
**not an official leaderboard score**.

No credential or token is stored anywhere in this repository.
