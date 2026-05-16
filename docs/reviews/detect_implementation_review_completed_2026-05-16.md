# Completed Detector Implementation Review

Date: 2026-05-16  
Scope: current `src/detect/*`, CLI wiring, detector tests, and the
`ground_truth_v2` benchmark output.

## Executive Summary

The detector implementation has moved from a partial M3/M4 prototype to a
substantially complete M6-style implementation. Several earlier review findings
were fixed: expanded divisor widths can now drive the final call, missing
periods are hard errors by default, diagnostics JSON is emitted, batch mode has
granular viz flags, unmatched period files are detected in batch mode, and the
final-width phase-shift recomputation bug was fixed.

The current benchmark result is strong: a fresh run over `ground_truth_v2/out`
with `target/release/kitehor detect-batch` produced `1511 / 1600 = 94.4%`,
above the planned 92% target. The remaining risks are not broad failures; they
cluster around clean HOR edge cases, same-width mixed arrays, segment semantics,
and confidence/metadata accuracy.

I would call the implementation usable for synthetic benchmark iteration, but
not yet ready to present `confidence`, `inter_monomer_identity`, or mixed/segment
calls as stable biological outputs without more work.

## Verification Performed

- `cargo test --release detect --tests` could not run locally:
  - local `cargo`: `1.75.0`
  - local `rustc`: `1.75.0`
  - repo requires `rust-version = "1.85"` and `rust-toolchain.toml` pins `1.95`
  - `Cargo.lock` is v4, so Cargo 1.75 fails before compilation.
- Existing benchmark output:
  - `python3 tools/detect_eval/eval.py --manifest ground_truth_v2/manifest.tsv --properties-dir ground_truth_v2/det_out`
  - result: `1511 / 1600 = 94.4%`.
- Fresh benchmark regeneration using the available release binary:
  - generated `/tmp/kitehor_detect_review_out`
  - re-evaluated with `--csv-out /tmp/kitehor_detect_review_eval.csv`
  - result reproduced exactly: `1511 / 1600 = 94.4%`.

Per-category fresh result:

| Category | Correct | Total | Accuracy |
|---|---:|---:|---:|
| `gc_bias` | 50 | 50 | 100.0% |
| `hor_clean` | 529 | 600 | 88.2% |
| `hor_event_deletion` | 49 | 50 | 98.0% |
| `hor_event_duplication` | 50 | 50 | 100.0% |
| `hor_event_hybrid` | 50 | 50 | 100.0% |
| `hor_event_inversion` | 50 | 50 | 100.0% |
| `hor_insertion` | 100 | 100 | 100.0% |
| `hor_shift` | 198 | 200 | 99.0% |
| `hor_wobble` | 95 | 100 | 95.0% |
| `mixed` | 96 | 100 | 96.0% |
| `random` | 50 | 50 | 100.0% |
| `simple_tr` | 194 | 200 | 97.0% |

## Implemented Improvements Since Prior Review

- Final classification now considers all supported `width_features`, including
  divisor-expanded widths, not only original input periods
  (`src/detect/classify.rs:79`).
- The final-width shift recomputation now stores `pre_decision_width` before
  overwriting the selected width (`src/detect/mod.rs:164`).
- Missing period rows are hard errors unless `--allow-missing-periods` is set
  (`src/detect/io.rs:127`).
- `width_features.tsv` now contains phase separation, irregularity, and class
  hints (`src/detect/mod.rs:301`).
- `diagnostics.json` is written by default (`src/detect/mod.rs:99`,
  `src/detect/io.rs:226`).
- `detect-batch` detects periods files without matching FASTA files unless
  `--allow-extra-periods` is used (`src/detect/mod.rs:485`).
- `detect-batch` now exposes the same granular viz flags as `detect`
  (`src/cli.rs:376`).

## Findings

### High: Final `mixed` and `ambiguous` calls can inherit a pre-decision width

`run_array_m3_5()` sets `props.base_width_bp` from the highest-scored input
period used for shift analysis (`src/detect/mod.rs:292`). Later, M4 merges the
classification result with:

```rust
props_m35.base_width_bp = decision.base_width_bp.or(props_m35.base_width_bp);
```

at `src/detect/mod.rs:172`.

That means a final `mixed` or `ambiguous` decision, whose `ArrayDecision`
intentionally has `base_width_bp = None`, can still publish a heuristic
pre-classification width. This is visible in benchmark failures such as
`mixed` calls with `base_width_bp=100` and `hor_k=NA`, or `ambiguous` calls with
reason `"no width achieves ic_threshold_rescue"` but a non-NA base width.

This also causes consensus and visualization output to be emitted for some
non-resolved calls because `run_one()` checks only `props.base_width_bp`
(`src/detect/mod.rs:79`).

Recommendation: after classification, treat `decision` as authoritative for
class-defining fields. Preserve M3.5 shift fields separately, but do not retain
`base_width_bp`, `hor_k`, `hor_length_bp`, `column_conservation`, or consensus
for `mixed`/`ambiguous` unless the decision explicitly supplies them.

### High: Segment rows are boundary splits, not per-segment classifications

`segment::split()` emits rows at phase-shift boundaries, but every segment
inherits the whole-array class, base width, `k`, conservation, phase separation,
wobble, and irregularity (`src/detect/segment.rs:33`). No per-segment wrap,
classification, consensus comparison, or stratification decision is performed.

This falls short of the plan's MVP segmentation intent: threshold-based
segmentation plus per-segment regime A/B/C classification. The config fields
for stratification thresholds exist (`src/detect/config.rs:71`), but they are
not used anywhere in detector logic.

Recommendation: implement a real segment recompute path:

1. split by phase-shift/domain-transition candidates,
2. rerun width features on each segment,
3. compare segment consensuses,
4. apply the planned same/different thresholds,
5. classify same-architecture shifts as properties and different-architecture
   blocks as `mixed` or `ambiguous`.

### High: Same-width mixed arrays are still collapsed into single HOR calls

The benchmark still has four `mixed` failures, all detected as `HOR`. The
examples are same-width/same-`k` mixed-family cases such as
`mx_a200-08_b200-08_n050-050`. The current mixed detection relies mostly on
distinct high-score input periods (`src/detect/classify.rs:365`) or incompatible
candidate widths. If two adjacent families share the same base width and HOR
multiplicity, that signal disappears and the detector sees one coherent HOR.

This is biologically important because adjacent satellite families can share
period and multiplicity while differing in consensus.

Recommendation: use segment-level consensus identity, not only period structure,
for mixed detection. The already-defined `stratification_same_threshold` and
`stratification_diff_threshold` should drive this. Same width/`k` plus consensus
identity below the threshold should become `mixed`, not clean `HOR`.

### High: Irregularity currently conflates wobble with irregular HOR

`irregularity::compute()` measures block-level column-IC variance
(`src/detect/irregularity.rs:20`). `run_array_m4()` then demotes any `HOR` to
`irregular_HOR` when that score crosses a threshold (`src/detect/mod.rs:181`).

In the fresh benchmark, four `hor_wobble` cases are demoted to
`irregular_HOR`. These are not random failures: strong smooth wobble degrades
block IC and looks "irregular" to this metric, even though the architecture is
still a coherent HOR with a wobble property.

Recommendation: compute irregularity after accounting for known smooth wobble
and phase shifts. Demotion should represent local architectural inconsistency,
not smooth spacing drift. At minimum, do not demote to `irregular_HOR` when the
dominant abnormality is high `wobble_amplitude_bp` and segment architecture is
unchanged.

### Medium: `inter_monomer_identity` is still not sequence identity

The final HOR decision sets:

```rust
inter_monomer_identity: Some(c.r_lag1)
```

at `src/detect/classify.rs:396`.

`r_lag1` is k-mer row similarity at the selected base width. It is useful, but
it is not mean sequence identity between HOR slots. The spec describes
`inter_monomer_identity` as a biological regime indicator. Reporting a k-mer
composition dot product under that name will mislead downstream interpretation,
especially when rows have compositional similarity without positional identity.

Recommendation: either rename the field in a future schema revision, or compute
the actual mean pairwise identity among inferred slot consensuses. If schema
cannot change, document the current value as an approximation and include a
separate calibrated identity later.

### Medium: Confidence is not calibrated as a probability

`confidence::compute()` uses a hand-weighted sigmoid. For `mixed` and
`ambiguous`, it returns a constant sigmoid of `-1.0`, about 0.2689
(`src/detect/confidence.rs:60`). Some wrong calls have very high confidence,
including true benchmark HORs called `simple_TR` with confidence `1.0000`.

This is acceptable as an internal score, but not as a user-facing probability.
The benchmark evaluates class only; it does not validate confidence calibration.

Recommendation: either label this column as a heuristic confidence score in user
docs, or calibrate it against held-out benchmark data with reliability curves
per class. For `mixed`/`ambiguous`, confidence should reflect evidence quality,
not a class constant.

### Medium: Clean HOR edge cases remain the main accuracy gap

Most remaining errors are in `hor_clean`: 71 of the 89 total failures. Of those,
47 are called `simple_TR` and 24 are called `mixed`.

Observed patterns:

- low-copy or low-divergence clean HORs can produce both HOR and simple-TR
  candidates, then fall into `mixed`;
- some high-divergence `d40` cases are classified as regime C `simple_TR` at
  the HOR-unit width, even though the benchmark expects `HOR`.

The implementation is close enough that broad threshold lowering is risky.
The failures need targeted calibration by `k`, copy count, length, and
inter-slot divergence.

Recommendation: add a failure-reporting table to the benchmark harness that
summarizes error rates by `L`, `k`, `n`, divergence, mutation, and indel rate.
Use that table before changing thresholds like `regime_c_r1_threshold`,
`phase_separation_threshold`, or `simple_tr_r1_rescue`.

### Medium: Period score parsing silently converts invalid values to zero

`load_periods()` parses `period_score` with:

```rust
let period_score: f64 = rec.get(i_score).unwrap_or("0").parse().unwrap_or(0.0);
```

at `src/detect/io.rs:88`.

Malformed scores should be rejected, not silently changed to zero. A typo in an
upstream periods file can therefore change candidate ranking and classification
without an error. The loader also does not validate score range or finiteness.

Recommendation: parse `period_score` with context, reject NaN/infinite values,
and require a documented range, probably `[0, 1]`. Add tests for malformed,
negative, and out-of-range scores.

### Medium: Viz flags are partly accepted but not honored

The CLI exposes `--export-shift`, `--export-edges`, and `--export-ic`
(`src/cli.rs:343`, `src/cli.rs:377`). `VizFlags` stores them
(`src/detect/viz.rs:28`), but `viz::export()` writes the cheap IC, edge-rate,
`R(k)`, and shift TSVs whenever `viz_dir` is present, regardless of those flags
(`src/detect/viz.rs:72`). Only `--export-raster` changes output.

Also, `--export-edges` does not export edge matrices; only
`column_edge_rate_w*.tsv` is available.

Recommendation: choose one behavior and make it exact. Either:

- `--viz-dir` writes all cheap TSVs and remove/ignore granular cheap flags from
  the documented interface, or
- make each flag control exactly its corresponding artifact.

### Low: Several config fields are currently inert

`embedding_dim_hash`, `shift_breakpoint_window_frac`,
`stratification_same_threshold`, and `stratification_diff_threshold` are exposed
in `DetectorConfig`, but are not used in detection logic. Some are planned
future hooks, but exposing inert knobs makes calibration harder because users
can change values with no effect.

Recommendation: either wire them up now or mark them as reserved/internal until
the related behavior exists.

### Low: Single-run `detect` ignores extra unmatched period rows

Batch mode checks for `.periods.tsv` files without matching FASTA stems
(`src/detect/mod.rs:485`), but single-run `detect` does not check for leftover
`array_id`s after `join_arrays_with_periods()` consumes matching rows. A typo in
`array_id` can be caught as "missing periods" for the FASTA record, but an extra
stale row can remain unnoticed if all real records matched.

Recommendation: after joining, fail or warn on unused period groups unless an
explicit `--allow-extra-periods`-style option is set for single-run mode too.

## Benchmark Failure Pattern

Fresh benchmark failures grouped by expected category:

| Category | Failure count | Detected as |
|---|---:|---|
| `hor_clean` | 71 | 47 `simple_TR`, 24 `mixed` |
| `hor_event_deletion` | 1 | 1 `mixed` |
| `hor_shift` | 2 | 2 `mixed` |
| `hor_wobble` | 5 | 4 `irregular_HOR`, 1 `simple_TR` |
| `mixed` | 4 | 4 `HOR` |
| `simple_tr` | 6 | 3 `mixed`, 3 `ambiguous` |

This is a useful profile. Random, GC bias, insertion, inversion, hybrid, and
duplication cases are effectively solved under the current synthetic oracle.
The remaining work is mostly classification-boundary quality, not pipeline
coverage.

## Suggested Priority

1. Fix final property merge semantics so `mixed` and `ambiguous` do not inherit
   stale `base_width_bp` or emit misleading consensus records.
2. Implement real per-segment recomputation and consensus identity comparison.
3. Use segment consensus identity to catch same-width mixed arrays.
4. Refine irregularity so smooth wobble remains a property of HOR rather than a
   demotion to `irregular_HOR`.
5. Add benchmark failure stratification by length, copy count, `k`, divergence,
   mutation, and indel rate before changing more thresholds.
6. Replace or relabel `inter_monomer_identity`.
7. Calibrate confidence or document it as a heuristic score.
8. Tighten `periods.tsv` validation and make viz/config flags exact.

## Bottom Line

This is a strong implementation relative to the plan: the main benchmark target
is met, the detector writes the expected output bundle, and many previously
missing pieces have landed. The remaining issues are mostly about truthfulness
of the reported biological properties and robustness outside the synthetic
acceptance harness. I would focus next on segment-aware classification and
output semantics rather than additional broad threshold tuning.
