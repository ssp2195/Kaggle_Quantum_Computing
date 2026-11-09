# Preparing the full Kaggle submission

Only Benchmark I (`invsqrt`) has been optimized by this branch. Keep the
competition-provided sample/reference `dqi` and `dfthcsos` files unchanged
unless separate independently verified improvements are available.

## New best fully verified Benchmark I result

**Compact measurement-corrected normalization + grouped quadratic phase:**

- Full public competition-verifier mirror: **PASS** on **200,540,928** shots.
- Benchmark I cost `Q*C`: **430,259.61862786254**.
- Peak qubits: **112**.
- Native Toffoli-like operations: **3,218**.
- Mean Toffoli-equivalent cost including rotations: **3,841.603737748773**.
- Fine rotations per shot: **88**.
- Worst phase-error / tolerance: **0.717347376787342**.
- Certified circuit SHA-256:
  `938d053f2218caa1859837ffb466ded01033c2d8d6445699276e5b4e21ebbc42`.
- [Full successful verifier and .kmb artifact](https://github.com/ssp2195/Kaggle_Quantum_Computing/actions/runs/42657329801).

This is not a new Kaggle leaderboard submission. The user previously
reported an actual combined Kaggle score of 1.33 from an earlier circuit.

## Reproduce or retrieve the certified KMB

Download the certified `.kmb` from the GitHub Actions run above, or
rebuild using the same pinned native kickmix version:

```bash
git clone https://github.com/quantumlib/Qualtran.git /tmp/Qualtran
git -C /tmp/Qualtran checkout 9a0c11a4ee7697c082fbf8f2e1e8bf4a71fcb453
python -m pip install /tmp/Qualtran/kickmix

python -m src.benchmark1_best --output out/best_invsqrt.kmb
sha256sum out/best_invsqrt.kmb
```

`src.benchmark1_best` validates resource counts and rejects any serialized
circuit whose SHA-256 differs from the full-verifier-certified artifact.
Different versions of kickmix may produce a different .kmb and require
independent complete verification.

## Assemble a three-task submission ZIP

Download official competition data after accepting the rules. You need
`sample_submission/dqi.kmb` and all
`sample_submission/dfthcsos/` files, including `submission.npz`.

```bash
python -m scripts.package_submission \
    --invsqrt out/best_invsqrt.kmb \
    --reference official \
    --output out/submission.zip

python -m zipfile -l out/submission.zip
```

Submit using a locally authenticated Kaggle CLI. Do not put tokens in
source code, shell history, repository secrets or public job logs.

```bash
kaggle competitions submit \
    -c google-verifiable-quantum-circuit-compilation \
    -f out/submission.zip \
    -m "112q MBUC-normalized invsqrt; reference DQI/DFTHCSOS"

kaggle competitions submissions \
    -c google-verifiable-quantum-circuit-compilation
```

**Expected leaderboard effect**, assuming the other two benchmark costs
stay at their references: `(2076880 / 430259.61862786254)**(1/3)`
≈ **1.69**, not an official measured score. The **30× Benchmark I**
stretch goal corresponds to `Q*C ≤ 69229.33`, still **6.215×** lower
than the present verified circuit.

No Kaggle authentication credential is committed to the repository.
