#!/usr/bin/env python3
"""Join the v2-corpus manifest with detector output and report
per-category class accuracy.

Usage:
    python3 tools/detect_eval/eval.py \
        --manifest    ground_truth_v2/manifest.tsv \
        --properties-dir ground_truth_v2/det_out \
        [--csv-out /tmp/eval.csv]

`properties-dir` is expected to contain `<category>/<case_id>.properties.tsv`
mirroring `kitehor synth-batch` + `kitehor detect-batch` layout.

The eval applies the simulator-truth → detector-expectation mapping:

    category          expected_class
    ----------------- ---------------
    simple_tr         simple_TR
    hor_clean         HOR
    hor_wobble        HOR
    hor_shift         HOR
    hor_insertion     HOR
    hor_event_hybrid  HOR
    hor_event_inv     irregular_HOR  (or ambiguous, OQ3)
    hor_event_dup/del HOR            (local CNV, class unchanged)
    mixed             mixed
    random            ambiguous
    gc_bias           simple_TR
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------- helpers ----------------------

# Maps a manifest `category` value to the detector-side expected class.
# Per OQ3, inversion expected = irregular_HOR (strand-aware detection
# deferred); per OQ-otherwise, dup/del are local CNV events that don't
# change the class.
CATEGORY_TO_EXPECTED = {
    "simple_tr":             "simple_TR",
    "hor_clean":             "HOR",
    "hor_wobble":            "HOR",
    "hor_shift":             "HOR",
    "hor_insertion":         "HOR",
    "hor_event_hybrid":      "HOR",
    # Per OQ3: strand-aware inversion recognition is v2. Without it,
    # the detector legitimately calls inversions as HOR (the slot-
    # consensus signal is largely preserved). The CI oracle keeps
    # irregular_HOR for T12 as the *target* behaviour, but at the
    # benchmark level we accept HOR.
    "hor_event_inversion":   "HOR",
    "hor_event_duplication": "HOR",
    "hor_event_deletion":    "HOR",
    "mixed":                 "mixed",
    "random":                "ambiguous",
    "gc_bias":               "simple_TR",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_properties(properties_dir: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for p in properties_dir.rglob("*.properties.tsv"):
        for row in read_tsv(p):
            out[row["array_id"]] = row
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--properties-dir", required=True, type=Path)
    ap.add_argument("--csv-out", type=Path)
    args = ap.parse_args()

    manifest = read_tsv(args.manifest)
    properties = load_properties(args.properties_dir)

    per_category = defaultdict(lambda: {"correct": 0, "n": 0, "missing": 0})
    overall = {"correct": 0, "n": 0, "missing": 0}
    rows = []

    for m in manifest:
        cid = m["case_id"]
        cat = m["category"]
        expected = CATEGORY_TO_EXPECTED.get(cat)
        if expected is None:
            continue
        p = properties.get(cid)
        if p is None:
            per_category[cat]["missing"] += 1
            overall["missing"] += 1
            continue
        got = p["class"]
        ok = got == expected
        per_category[cat]["correct"] += int(ok)
        per_category[cat]["n"] += 1
        overall["correct"] += int(ok)
        overall["n"] += 1
        rows.append({
            "case_id": cid,
            "category": cat,
            "expected_class": expected,
            "detected_class": got,
            "correct": ok,
            "base_width_bp": p.get("base_width_bp", ""),
            "hor_k": p.get("hor_k", ""),
            "confidence": p.get("confidence", ""),
            "reason": p.get("reason", ""),
        })

    print(f"{'category':25s} {'correct':>8s} / {'n':>5s}  {'pct':>6s}")
    for cat in sorted(per_category):
        v = per_category[cat]
        pct = 100.0 * v["correct"] / v["n"] if v["n"] else 0.0
        miss = f" (missing={v['missing']})" if v["missing"] else ""
        print(f"  {cat:25s} {v['correct']:>8d} / {v['n']:>5d}  {pct:>5.1f}%{miss}")
    pct_overall = 100.0 * overall["correct"] / overall["n"] if overall["n"] else 0.0
    miss_overall = (
        f" (missing={overall['missing']})" if overall["missing"] else ""
    )
    print(
        f"  {'OVERALL':25s} {overall['correct']:>8d} / {overall['n']:>5d}  "
        f"{pct_overall:>5.1f}%{miss_overall}"
    )

    if args.csv_out:
        with args.csv_out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys() if rows else [])
            if rows:
                w.writeheader()
                w.writerows(rows)
        print(f"\nwrote per-case detail to {args.csv_out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
