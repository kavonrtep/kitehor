#!/usr/bin/env bash
# Run `kitehor synth-batch` over every category in ground_truth_v2/configs/
# and write outputs to ${SYNTH_OUT:-./out}, one subdir per category.
#
# Output bundles per array:  {case_id}.fa, .truth.tsv, .periods.tsv
# Add --diagnostics manually if you need per-case .diagnostics.json.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG_ROOT="${HERE}/configs"
OUT_ROOT="${SYNTH_OUT:-${HERE}/out}"
BIN="${KITEHOR:-${HERE}/../target/release/kitehor}"

if [[ ! -x "${BIN}" ]]; then
    echo "error: ${BIN} not found or not executable" >&2
    echo "build first:  cargo build --release  (from the kitehor crate root)" >&2
    exit 1
fi

mkdir -p "${OUT_ROOT}"
echo "binary  : ${BIN}"
echo "configs : ${CFG_ROOT}"
echo "outputs : ${OUT_ROOT}"
echo

total_t0=$(date +%s)
for cat_dir in "${CFG_ROOT}"/*/; do
    cat_name="$(basename "${cat_dir}")"
    out_dir="${OUT_ROOT}/${cat_name}"
    mkdir -p "${out_dir}"
    n_in="$(find "${cat_dir}" -maxdepth 1 -name '*.yaml' | wc -l)"
    t0=$(date +%s)
    "${BIN}" synth-batch \
        --config-dir "${cat_dir}" \
        --out-dir    "${out_dir}" \
        "$@"
    t1=$(date +%s)
    n_fa="$(find "${out_dir}" -maxdepth 1 -name '*.fa' | wc -l)"
    printf "  %-20s configs=%-5d fa=%-5d wall=%ds\n" \
        "${cat_name}" "${n_in}" "${n_fa}" "$((t1 - t0))"
done
total_t1=$(date +%s)

n_total=$(find "${OUT_ROOT}" -name '*.fa' | wc -l)
echo
echo "TOTAL: ${n_total} bundles in $((total_t1 - total_t0))s"
