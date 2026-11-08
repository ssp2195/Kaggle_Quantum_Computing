#!/usr/bin/env bash
# Run Google's authoritative Benchmark-I verifier on a native kickmix .kmb.
set -euo pipefail
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 CANDIDATE.kmb [OFFICIAL_DATA_DIR]" >&2
  exit 2
fi
candidate="$(realpath "$1")"
official="${2:-official}"
verifier="$official/verifiers/invsqrt/verify_submission_invsqrt.py"
if [[ ! -f "$verifier" ]]; then
  echo "Missing official verifier: $verifier" >&2
  echo "Accept the Kaggle rules, download competition data, and extract it first." >&2
  exit 2
fi
if [[ ! -f "$candidate" ]]; then
  echo "Missing candidate .kmb: $candidate" >&2
  exit 2
fi
# The verifier checks the output, scratch cleanup, phase error and
# kickmix-resource score including fine Z_POW synthesis.
python -I "$verifier" "$candidate" \
  --n-p 10 --beta 1 --precision-bits 24 --eta 218 \
  --max-x "${MAX_X:-3133452}" --shots-per-x 64 \
  --batch-size 256 --max-workers "${MAX_WORKERS:-4}" \
  --salt "${VERIFIER_SALT:-quantum-circuit-local-v1}"
