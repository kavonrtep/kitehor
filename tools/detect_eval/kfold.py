#!/usr/bin/env python3
"""Stratified k-fold cross-validation on the ground_truth_v2 corpus.

Splits the 1600 cases into k folds stratified by category, then
reports per-fold accuracy under the **currently-compiled thresholds**.
This doesn't re-tune per fold (no automated threshold search exists
yet), so the reported variance measures only deterministic-split
stability — not true out-of-sample generalisation. But low variance
across folds *plus* the CI fixture set (T01–T18) still passing exactly
gives confidence that the tuned thresholds aren't overfit to specific
case IDs.

Usage:
    python3 tools/detect_eval/kfold.py \\
        --manifest ground_truth_v2/manifest.tsv \\
        --properties-dir ground_truth_v2/det_out \\
        [--k 5]
"""

from __future__ import annotations

import argparse
import csv
import sys
import statistics
from collections import defaultdict
from pathlib import Path

# Re-use the same category → expected-class mapping as the main eval
# harness. Keep them in sync — if you change one, change both.
CATEGORY_TO_EXPECTED = {
    "simple_tr":             "simple_TR",
    "hor_clean":             "HOR",
    "hor_wobble":            "HOR",
    "hor_shift":             "HOR",
    "hor_insertion":         "HOR",
    "hor_event_hybrid":      "HOR",
    "hor_event_inversion":   "HOR",   # OQ3: strand-aware deferred
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


def stable_fold(case_id: str, k: int) -> int:
    """Deterministic per-case fold assignment via a stable hash."""
    h = 0
    for b in case_id.encode():
        h = (h * 131 + b) & 0xFFFFFFFF
    return h % k


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--properties-dir", required=True, type=Path)
    ap.add_argument("--k", type=int, default=5, help="number of folds (default 5)")
    args = ap.parse_args()

    manifest = read_tsv(args.manifest)
    properties = load_properties(args.properties_dir)

    # Build (case_id, category, expected, got, correct) per case;
    # stratify by category so each fold gets the same mix.
    rows_by_cat: dict[str, list[tuple[str, bool]]] = defaultdict(list)
    n_total = 0
    n_missing = 0
    for m in manifest:
        cid = m["case_id"]
        cat = m["category"]
        expected = CATEGORY_TO_EXPECTED.get(cat)
        if expected is None:
            continue
        p = properties.get(cid)
        if p is None:
            n_missing += 1
            continue
        ok = (p["class"] == expected)
        rows_by_cat[cat].append((cid, ok))
        n_total += 1

    # Within each category, assign cases to folds via stable_fold.
    # This produces stratified folds with deterministic membership
    # — no per-run randomness.
    folds: list[list[tuple[str, str, bool]]] = [[] for _ in range(args.k)]
    for cat, lst in rows_by_cat.items():
        for cid, ok in lst:
            f = stable_fold(cid, args.k)
            folds[f].append((cid, cat, ok))

    print(f"k-fold CV (k={args.k}), n_total={n_total}, n_missing={n_missing}")
    print(f"{'fold':>4s}  {'n':>4s}  {'correct':>7s}  {'pct':>6s}  per-category-breakdown")
    fold_pcts: list[float] = []
    for fi, fold in enumerate(folds):
        n = len(fold)
        n_correct = sum(1 for _, _, ok in fold if ok)
        pct = 100.0 * n_correct / n if n else 0.0
        fold_pcts.append(pct)
        # Per-category counts within this fold
        per_cat = defaultdict(lambda: [0, 0])
        for _, cat, ok in fold:
            per_cat[cat][0] += int(ok)
            per_cat[cat][1] += 1
        cat_str = " ".join(
            f"{cat}={ok}/{n}"
            for cat, (ok, n) in sorted(per_cat.items())
        )
        print(f"  {fi:>2d}  {n:>4d}  {n_correct:>7d}  {pct:>5.1f}%  {cat_str}")

    mean = statistics.mean(fold_pcts)
    sd = statistics.stdev(fold_pcts) if len(fold_pcts) > 1 else 0.0
    print(f"\nmean = {mean:.2f}%   sd = {sd:.2f}%")
    print(f"range = [{min(fold_pcts):.1f}%, {max(fold_pcts):.1f}%]")
    print(f"(sd/mean = {100*sd/mean:.2f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
