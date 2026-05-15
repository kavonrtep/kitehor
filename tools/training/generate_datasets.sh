#!/usr/bin/env bash
# Generate 7 simulated datasets (5 train + 2 holdout) using the filtered
# params.tsv. Each dataset gets its own master seed; per-case PRNG seeds
# are deterministically derived via FNV-1a of "master:case_id".
#
# Output: eval/training_data/sim_seed${seed}/{sequences.fasta,truth.tsv,...}

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

PARAMS="eval/training_data/params_filtered.tsv"
BIN="./hordetect/target/release/hordetect"
OUTDIR_BASE="eval/training_data"

if [[ ! -x "$BIN" ]]; then
    echo "error: hordetect release binary not found at $BIN" >&2
    exit 1
fi
if [[ ! -f "$PARAMS" ]]; then
    echo "error: $PARAMS missing — run build_filtered_params.R first" >&2
    exit 1
fi

TRAIN_SEEDS=(101 102 103 104 105)
TEST_SEEDS=(901 902)

for seed in "${TRAIN_SEEDS[@]}" "${TEST_SEEDS[@]}"; do
    outdir="${OUTDIR_BASE}/sim_seed${seed}"
    if [[ -d "$outdir" && -f "$outdir/sequences.fasta" ]]; then
        echo "[skip] $outdir already populated"
        continue
    fi
    echo "[seed=$seed] generating → $outdir"
    "$BIN" simulate-grid --params "$PARAMS" --outdir "$outdir" --seed "$seed"
done

echo "done. 5 train + 2 holdout datasets in $OUTDIR_BASE/sim_seed*"
