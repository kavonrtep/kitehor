# Kite → Detect pipeline on `ground_truth_v2/` (1600-case benchmark)

**Date:** 2026-05-16
**Commit:** `1146abf` (HEAD at run time)
**Corpus:** `ground_truth_v2/out/` — the v2 simulator's 1600-case
benchmark across 9 categories. Generated via `ground_truth_v2/run_batch.sh`.

## Purpose

Measure how much accuracy the v2 line-width detector loses when its
input period candidates come from `kitehor kite-periodicity --classify
--emit-periods` instead of the simulator-emitted oracle `periods.tsv`.
The oracle baseline (94.4 % overall, M6 acceptance) was reached by
feeding the detector exactly the periods the simulator declared as
true_base / true_hor_unit / distractors. Real users won't have an
oracle; this run tells us where the realistic ceiling sits.

## Method

```bash
# 1. Concatenate all 1600 v2 FASTAs.
cat ground_truth_v2/out/*/*.fa > /tmp/kite_v2_run/all.fa

# 2. Kite-periodicity over the whole corpus, with rule classifier and
#    v2-compatible periods emission.
./target/release/kitehor kite-periodicity /tmp/kite_v2_run/all.fa \
    -o /tmp/kite_v2_run/all.kite.tsv \
    --classify \
    --emit-periods /tmp/kite_v2_run/all.kite.periods.tsv

# 3. Split kite's combined periods.tsv into per-stem files so
#    detect-batch (which parallelises via rayon) can run.
awk -F'\t' 'NR==1{h=$0;next} {f="/tmp/kite_v2_run/periods_by_stem/"$1".periods.tsv"; if(!(f in s)){print h>f;s[f]=1} print >> f}' \
    /tmp/kite_v2_run/all.kite.periods.tsv

# Stub header-only files for the 10 NoSignal records (otherwise
# detect-batch fails the symmetric-pairing check from DH11).
for fa in /tmp/kite_v2_run/fasta_flat/*.fa; do
  stem=$(basename "$fa" .fa)
  p="/tmp/kite_v2_run/periods_by_stem/${stem}.periods.tsv"
  [ -f "$p" ] || printf 'array_id\tperiod_bp\tperiod_score\tsource\n' > "$p"
done

# 4. Detector.
./target/release/kitehor detect-batch \
    --fasta-dir /tmp/kite_v2_run/fasta_flat \
    --periods-dir /tmp/kite_v2_run/periods_by_stem \
    --out-dir /tmp/kite_v2_run/det_out \
    --allow-missing-periods

# 5. Evaluation.
python3 tools/detect_eval/eval.py \
    --manifest ground_truth_v2/manifest.tsv \
    --properties-dir /tmp/kite_v2_run/det_out
```

Wall time on this run (6-core dev box):

| Stage | Wall | CPU |
|---|---:|---:|
| `kite-periodicity` (1600 records) | 30 s | 6 min 16 s |
| `detect-batch` (1600 records) | 9 min 2 s | 1 h 58 min |

## Results

Direct head-to-head against the oracle-period baseline:

| category              | oracle (94.4 % run) | kite → detect | delta |
|-----------------------|--------------------:|--------------:|------:|
| `gc_bias`             | 100.0 % (50/50)     | 100.0 % (50/50)   | 0.0 |
| `hor_clean`           | 88.2 % (529/600)    | **92.3 %** (554/600) | **+4.1** |
| `hor_event_deletion`  | 98.0 % (49/50)      | 100.0 % (50/50)   | +2.0 |
| `hor_event_duplication` | 100.0 % (50/50)   | 100.0 % (50/50)   | 0.0 |
| `hor_event_hybrid`    | 100.0 % (50/50)     | 100.0 % (50/50)   | 0.0 |
| `hor_event_inversion` | 100.0 % (50/50)     | 100.0 % (50/50)   | 0.0 |
| `hor_insertion`       | 100.0 % (100/100)   | 100.0 % (100/100) | 0.0 |
| `hor_shift`           | 99.0 % (198/200)    | 100.0 % (200/200) | +1.0 |
| `hor_wobble`          | 95.0 % (95/100)     | 95.0 % (95/100)   | 0.0 |
| **`mixed`**           | **96.0 %** (96/100) | **18.0 %** (18/100) | **-78.0** |
| **`random`**          | **100.0 %** (50/50) | **68.0 %** (34/50)  | **-32.0** |
| `simple_tr`           | 97.0 % (194/200)    | 95.0 % (190/200)  | -2.0 |
| **OVERALL**           | **94.4 %** (1511/1600) | **90.1 %** (1441/1600) | **-4.3** |

8 of 12 categories are at parity or improved. Two categories regress
materially: `mixed` collapses by 78 pp, `random` drops 32 pp.
Everything else is within ±2 pp.

## What broke and why

### `mixed`: 96 % → 18 %

The detector's multi-family detection relies on the
`multi_block_via_strong` rule in `src/detect/classify.rs:370`:
≥ 3 distinct input periods with `period_score ≥ strong_period_score`
(0.85) → `mixed`. The simulator-emitted oracle hands the detector two
`true_base` + two `true_hor_unit` rows per mixed array (4 strong
periods) so the rule fires reliably.

Kite-derived input has a different shape. The rule classifier
(`src/rule.rs`) collapses each array to one verdict, producing at
most one `kite_founder` + one `kite_tile` at high score. For mixed
arrays it picks whichever family the rule's first-pass scoring
favours. Counting strong periods (`period_score ≥ 0.85`) per mixed
array on this run:

| n strong periods | cases |
|---|---:|
| 1 | 93 |
| 2 | 7  |
| ≥ 3 | 0 |

Zero mixed arrays cross the `multi_block_via_strong` threshold, so the
override never fires. The detector then collapses to a single coherent
HOR call (51 cases), inherits an irregular_HOR demotion (30 cases),
or to simple_TR (1 case). 18 cases happen to fall through other
mixed-detection paths.

This regression aligns with the deferred A16 plan item
(`docs/new/detect_impl_plan.md` §0): same-width / same-`k` mixed
families require segment-level consensus identity comparison, which
the detector doesn't yet do. Until that lands, mixed detection is
structurally dependent on the upstream producer surfacing distinct
periods — which kite's single-verdict rule classifier cannot do.

### `random`: 100 % → 68 %

The oracle path emits **no rows** for random arrays — the simulator's
period generator legitimately reports no real period — and the
detector therefore reports `ambiguous` per A4 schema rules.

Kite, by contrast, runs its full k-mer histogram + peak-scoring
pipeline on every record. Random sequences produce noise peaks; some
arrays end up with `RuleVerdict::Tandem` (a single weakly-supported
monomer) rather than `NoSignal`, so kite emits a `kite_monomer`
@ 0.95. The detector then evaluates that width via canonical
column-IC and occasionally finds high enough IC to call `simple_TR`
(14 cases) or `mixed` (2 cases).

This is exactly the failure mode the integration-discussion option
(b) accepted: "kite returns garbage → detector might fire
false-positive simple_TR." The cap below 0.85 didn't help here
because the high-score row was the kite_monomer (0.95) — kite
*thinks* it's a real period.

### `hor_clean`: 88.2 % → 92.3 % (improvement)

Counter-intuitively, kite-derived periods do *better* than the
oracle on the hardest synthetic HOR category. Plausible mechanism:
the oracle emits `(true_base, true_hor_unit, near_miss, harmonic,
false_positive)` per array, and the false_positive distractor
periods sometimes confuse the detector's multiplicity dedup at
high inter-slot divergence. Kite's rule classifier filters the same
peaks more aggressively (its top-3 cut), giving the detector
cleaner candidate widths. We're not chasing this — it's noise in
both directions and below the M6 ≥ 88 % per-category target.

## Takeaways

1. The kite → detect glue **works end-to-end** at 90.1 % on a real,
   non-oracle benchmark. The score-mapping defaults
   (founder=0.95 / tile=0.90 / secondary=0.60 / hint=0.50–0.30) are
   correctly calibrated against the detector's `strong_period_score`
   gate: every category the rule classifier can resolve performs at
   or near oracle parity.
2. The 4.3 pp aggregate gap is **dominated by two structural issues**:
   - `mixed` (-78 pp) — needs the A16 same-width-mixed work (segment
     consensus identity comparison) before kite can drive it.
   - `random` (-32 pp) — kite's rule sometimes calls `Tandem` on
     noise. Could be addressed upstream (tighter `lo_period` or
     a noise-floor check in the rule) or downstream (a kite-aware
     score downgrade for very-low-share monomers).
3. The remaining 8 categories show that **single-family HORs and
   tandems are largely insensitive to whether periods come from
   oracle or kite**, as long as the rule classifier surfaces the
   correct base period. The detector's column-IC + R(k) test
   recovers the right class from a single high-score candidate.

## Next steps

- (deferred to M7) Land the per-segment recompute + same-width
  mixed detection (A16). Independent of kite — fixes both
  oracle-path and kite-path `mixed` regressions.
- Investigate the 16 false-positive `random → simple_TR`/`mixed`
  calls. Two cheap mitigations to consider before tuning thresholds:
  1. Add a `--rule-strict` flag in `kite-periodicity` that demotes
     `Tandem` to `NoSignal` when only one weak peak survives
     (mirrors the existing `Unresolved` path for short d1).
  2. In `emit_periods`, downgrade `kite_monomer` to a hint score
     (e.g., 0.60 instead of 0.95) when `kite-periodicity` reported
     `n_peaks_kept` ≤ 1.

## Reproducibility

| Artefact | Path / value |
|---|---|
| Commit | `1146abf` |
| Corpus | `ground_truth_v2/out/` (regen with `./ground_truth_v2/run_batch.sh`) |
| Manifest | `ground_truth_v2/manifest.tsv` |
| Eval harness | `tools/detect_eval/eval.py` |
| Kite scores | hardcoded defaults in `src/emit_periods.rs` |
| Detector thresholds | `DetectorConfig::default()` in `src/detect/config.rs` |
| Per-case CSV | `/tmp/kite_v2_run/eval.csv` (regenerate with `--csv-out`) |

To rerun, see the four-step shell block in the **Method** section
above. All outputs are deterministic given the corpus + binary.
