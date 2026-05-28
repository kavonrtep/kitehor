# Onboarding — `rescore` and `report` pipelines

This document is a single-file primer to two `kitehor` subcommands that
operate on tandem-repeat array FASTA input:

- **`kitehor rescore`** — augments kite's peaks TSV with nucleotide-level
  identity statistics, shift diagnostics, two derived flags (`phantom`,
  `subrepeat`), a spatial-coherence statistic (`spatial_contrast`), and
  the per-record founder period (`founder_period`). Per-peak refinement
  of kite's k-mer signal.
- **`kitehor report`** — runs kite, peak clustering, SSR scan, and
  irregularity scan on a FASTA and emits a single TSV with raw
  measurements per record. Whole-array observation-only mode; no
  categorical verdicts.

Detailed CLI references live in [`docs/rescore.md`](rescore.md) and
[`docs/report.md`](report.md). This document is the orientation map.

## TL;DR — which one do I run?

| You want… | Run |
|---|---|
| To refine kite's per-peak signal with pairwise nucleotide identity | `rescore` |
| To know which kite peaks are real HOR-unit periods vs. monomer harmonics vs. sub-period artifacts vs. internal sub-motifs | `rescore` |
| A whole-array snapshot — kite peaks, SSR coverage, irregularity metrics — without any classification cascade | `report` |
| Filterable per-record measurements with no opinionated calls | `report` |
| The opinionated cascade (`combined_class`, rule-classify verdicts) | `analyze` (not covered here; see [`rule_proto.md`](rule_proto.md)) |

`rescore` and `report` are **siblings**, not stages in a pipeline. They
consume FASTA independently. `rescore` additionally consumes kite's
peaks TSV (so the typical workflow is `kite-periodicity` → `rescore`),
while `report` runs kite internally.

## Status at a glance

| Stage | Maturity | Recently changed | Known caveats |
|---|---|---|---|
| `rescore` | **stable** — eight features shipped (Tier 1 + phantom + subrepeat Step A & B + auto-band + founder gate + spatial-contrast gate + period/founder ratio gate). 72 unit tests + 2 integration tests. | spatial-contrast gate (`--subrepeat-spatial-contrast-min`, default 0.40) + period/founder ratio gate (`--subrepeat-period-founder-max-ratio`, default 0.25); `founder_period` column exposed | `--max-period 5000` default skips very long candidates; subrepeat detection floor at K=200 is ≈ 10 % array footprint |
| `report` | **stable** — observation-only; no recent changes to columns or semantics. | none (last touch: `irreg_dropout_rate_per_pair` added in v0.12) | tandem-validate deliberately excluded; depends on rank-1 kite peak for downstream stage input |

Both are deterministic given inputs and CLI flags. Both run records in
parallel via rayon.

## Pipeline 1 — `rescore`

### Purpose

Kite scores periodicity by k-mer set overlap in the neighbour-distance
histogram. That signal can fail to separate:

- HOR-unit period vs. monomer period in well-formed HORs (the k-mer pool
  is shared between the two scales);
- a true period vs. a sub-period harmonic (e.g. a 56-bp candidate that's
  actually the boundary-shifted real 62-bp period);
- a real period vs. a localized subrepeat (a short tandem motif inside
  the founder monomer).

`rescore` runs banded semi-global edit distance on sampled adjacent
tile pairs to get a nucleotide-level identity per peak, then derives two
boolean flags (`phantom`, `subrepeat`) from the per-pair distribution and
alignment-shift statistics.

The output is the input peaks TSV with **13 columns appended**. Nothing
upstream is mutated. The metric is additive — downstream stages
(rule-classify, analyze) still drive decisions from kite's `score2_norm`.

### Inputs

- **FASTA file(s)** (positional) — tandem-repeat arrays.
- **`--peaks <peaks.tsv>`** — long-format kite peaks TSV (output of
  `kitehor kite-periodicity`). Header must contain `case_id`, `rank`,
  `period`.

### Output schema (full reference)

`<prefix>.peaks.tsv` is the input file with 13 columns appended:

| # | column | type | meaning |
|---|---|---|---|
| 1 | `identity_med` | float | median pairwise identity across K sampled adjacent-tile pairs |
| 2 | `identity_iqr` | float | interquartile range of pair identities (Q75 − Q25) |
| 3 | `identity_p25` | float | 25th-percentile identity — worst-quartile sentinel |
| 4 | `identity_n` | int | effective sample count (≤ K after N-rejection) |
| 5 | `shift_med` | int (bp) | median alignment shift between best DP path and natural mapping, over high-identity pairs only |
| 6 | `shift_consistency` | float | fraction of high-identity pairs with shift within ±1 bp of `shift_med` |
| 7 | `phantom` | bool | candidate period is a sub-tile of a longer real period |
| 8 | `subrepeat` | bool | candidate period is a localized short motif inside the founder monomer |
| 9 | `coverage_frac` | float | fraction of pairs with identity ≥ `--coverage-threshold`; independent diagnostic |
| 10 | `spatial_contrast` | float | max − min per-bin hit fraction across 10 anchor-offset bins; high (≈ 1) = localised subrepeat, low (≈ 0) = near-founder harmonic. `NA` when fewer than 2 bins meet the per-bin minimum |
| 11 | `founder_period` | int (bp) | per-record founder period used by the founder gate (same value across all rows of one record); `NA` when no row met `identity_med ≥ founder_id_min` and `phantom != true` |
| 12 | `kmer_autocorr_founder` | float | Pearson autocorrelation of the period-P k-mer pair density profile at lag = `founder_period`. Range `[−1, +1]`. High when density(x) oscillates with the founder period (real nested subrepeat). **Observational** — does not gate the `subrepeat` flag. `NA` when founder unknown |
| 13 | `kmer_phase_contrast` | float | Max-contiguous-half-fraction excess of midpoints folded by `(mid mod founder_period)` into 12 phase bins. Range `[0, 0.5]`. High when midpoints prefer one half of the founder cycle (TRC_104-style). **Observational** — does not gate the `subrepeat` flag. `NA` when founder unknown |

`NA` cells in any of the additive columns indicate either kite-filtered
rows (rank/period out of range or record missing from FASTA) or
kernel-NA (array too short, all-N sample windows). The first 9 columns
of the kite input pass through verbatim.

**Flag invariant**: `phantom = true` and `subrepeat = true` are mutually
exclusive — the founder gate enforces it.

### What the flags mean

- **`phantom = true`** — the kernel finds high identity, but the optimal
  alignment is systematically shifted by `shift_med` bp from the natural
  position. Evidence that the claimed period P is actually a sub-tile of
  a longer real period P + shift. Worked example: TRC_755 P=56 in the
  IPIP corpus, where `shift_med = +6` correctly identifies the real
  62 bp period.
- **`subrepeat = true`** — bimodal identity distribution (some sampled
  pairs hit hard, others miss) **AND** the high-identity hits are
  spatially clustered (`spatial_contrast ≥ 0.40`) **AND** the period
  is short enough to tile multiple times inside the founder
  (`period ≤ founder_period · 0.25` by default — tiles ≥ 4 times).
  Together these identify a short motif tiling only part of the
  array, not a near-founder harmonic. The selected founder period
  is exposed in the `founder_period` column. Worked example:
  TRC_104 P=36 inside founder P=180 (ratio 0.20).

### Algorithm (sketch — see `rescore.md` for full detail)

For each (record, period) row of the input peaks TSV:

1. Sample K=200 anchor offsets across the array. Pair the anchor tile
   with the next tile (with ±slop boundary slack).
2. Skip-pair: drop pairs whose combined N fraction exceeds
   `--max-n-frac`; redraw up to `--max-retries` times.
3. Compute banded semi-global edit distance per pair. Band auto-scales
   to `max(20, 2·slop, ⌈0.02·P⌉)` — short periods stay at the slop
   floor, long monomers get a band proportional to the period so the
   DP doesn't saturate on realistic indel rates.
4. Aggregate per-pair identities (median, IQR, p25) and shifts
   (median over high-identity pairs).
5. Apply heuristic gates to derive `phantom` and `subrepeat`.
6. Post-pass: founder gate forces `subrepeat → false` on any row with
   `period ≥ founder period`, where the founder is the lowest-rank row
   with `identity_med ≥ subrepeat_founder_id_min` and `phantom != true`.

### Calibration

On the synthetic 1600-case `ground_truth_v2/` corpus (defaults: K=200,
auto-band, founder-gated subrepeat):

| metric | value |
|---|---|
| Phantom FPs on true HOR-unit periods | 0 / 1313 |
| Phantom FPs on true monomer periods | 8 / 1576 (0.5 %) |
| Subrepeat FPs on true HOR-unit periods | 0 / 1313 |
| Subrepeat FPs on true monomer periods | 29 / 1576 (1.8 %) |
| HOR-unit identity > monomer identity (clean HORs) | 100 % |
| HOR-unit identity > monomer identity (mixed HORs) | 67 % |

On the real-data IPIP 2026-04-14 corpus (3024 centromeric arrays,
305 MB, ~17 s kite + ~165 s rescore):

| | count | of rescored |
|---|---|---|
| Rescored rows | 21,715 | 100 % |
| `phantom = true` | 306 | 1.4 % |
| `subrepeat = true` | 909 | 4.2 % |
| Either flag | 1,215 | 5.6 % |
| Both true (must be 0) | 0 | — |

### CLI knobs that matter most

- `--samples K` (default 200) — sampled pairs per row. Linear in cost,
  diminishing returns past K=200 for typical periods.
- `--top-n N` (default 10) — limits rescoring to top-N peaks per record.
  Combined with `--max-period 5000`, this caps wall time on long-tail
  peak lists.
- `--band 0` (auto) — period-relative band. Set to a non-zero integer
  to override the formula on the whole run.
- All threshold knobs for flag tuning (`--subrepeat-*`, `--shift-*`) have
  sensible defaults; calibrate by re-running the IPIP eval and comparing
  flag rates after a change.

### When to use it

- **Always after `kite-periodicity` on a new dataset** if you want a
  cheap correctness check on the kite peak ranking.
- **Diagnostic on suspected mis-calls** — read the four numeric identity
  columns first, then the two flags, then `shift_med` if `phantom` fires
  or `coverage_frac` if you're uncertain about a borderline row.

## Pipeline 2 — `report`

### Purpose

A defensive, observation-only mode that lays out four families of
measurements for one FASTA record:

1. kite peaks (raw and clustered);
2. SSR motif coverage (TideCluster-style);
3. irregularity / indel-event metrics;
4. array-level metadata.

No verdicts, no `combined_class`, no rule-classify labels. Designed for
when you don't trust the cascade's classifications and want raw numbers
to filter on yourself.

### Inputs

- **FASTA file** (positional, single file).

That's it — kite and all downstream stages run internally with their
own defaults (configurable via the CLI).

### Output schema (full reference — 20 columns)

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `kite_n_peaks` | int | number of peaks kite kept |
| 4 | `kite_peaks` | list-str | `period:score2_norm` entries, `;`-separated, sorted by score desc |
| 5 | `kite_n_clusters` | int | number of period clusters after `cluster_peaks(tol=0.015)` |
| 6 | `kite_clusters` | list-str | `rep_period:total_score:n_peaks` entries |
| 7 | `ssr_total_coverage_pct` | float | sum of all SSR motifs' coverage (cap 100 %) |
| 8 | `ssr_dominant_motif` | str | top SSR motif by coverage (canonical, upper) |
| 9 | `ssr_dominant_motif_length` | int | length in bp of the dominant motif |
| 10 | `ssr_dominant_motif_repeats` | int | total dominant-motif repeat units in the array |
| 11 | `ssr_dominant_coverage_pct` | float | dominant motif's coverage (% of array) |
| 12 | `ssr_top_motifs` | list-str | top-3 motifs: `motif:pct;motif:pct;motif:pct` |
| 13 | `irreg_flag` | str | scan outcome — `ok`, `too_short`, `no_kite`, … |
| 14 | `irreg_n_kmer_groups` | int | independent phase-bin groups used by the scan |
| 15 | `irreg_indel_event_count` | int | discrete indel events detected |
| 16 | `irreg_indel_burden_pct` | float | `Σ\|step_size\| / array_length × 100` |
| 17 | `irreg_indel_max_shift_bp` | float | size of the largest detected event |
| 18 | `irreg_indel_drift_bp_per_kb` | float | net cumulative drift per kb |
| 19 | `irreg_dropout_rate_per_pair` | float | fraction of expected k-mer occurrences missing |
| 20 | `irreg_baseline_jitter_bp` | float | per-array MAD-based residual scale |

Lists inside a cell use `;` between entries and `:` between fields
within an entry. Empty cells (not `NA`) indicate missing/non-applicable
values; the empty-string convention keeps column types numeric where
possible.

### Internal stages

`report` calls the following library functions in order (no parallelism
across stages; rayon parallelises across records):

1. `kite::analyze(record)` — k-mer immediate-neighbour histogram, peak
   detection, kite's filtering.
2. `rule_classify::cluster::cluster_peaks(peaks, tol)` — single-linkage
   clusters of peaks within `±cluster_tol` (default 0.015 = ±1.5 %).
3. `ssr::scan::scan_record(record, top_period)` — TideCluster-style
   motif scan; top_period comes from the rank-1 kite peak.
4. `irregularity::analyse_record(record, top_period)` — register-locked
   k-mer phase-shift scan for indel events.

The rank-1 peak (highest `score2_norm`) drives both SSR and irregularity
stages, matching `analyze.rs`'s choice. **Tandem-validate is
deliberately excluded** — its output is a categorical verdict, not a
measurement, and `report` is observation-only.

### Calibration

`report` is observation-only — there's no truth label to score against.
Field calibration is done by spot-checking against dotplots: the
`kite_peaks` / `kite_clusters` columns should match dotplot diagonal
spacings; SSR coverage should match the visible short-motif content;
irregularity metrics should match visible boundary shifts in tile
spacing.

Per-record wall time on the IPIP corpus: ~20 ms median, dominated by
kite (k-mer histogram) on long arrays.

### When to use it

- **As the default reporting mode** when downstream consumers want to
  filter on raw numbers rather than accept the cascade's calls.
- **As a quick health check** on a new array — read `kite_n_clusters`,
  `ssr_total_coverage_pct`, and `irreg_flag` to triage.

## Reading rows — interpretation guide

### From `report.tsv` — first triage

1. Look at `kite_n_clusters` — small means clear periodicity, large
   means kite finds many candidates.
2. Look at `ssr_total_coverage_pct` — > 80 % means the array is
   dominated by short SSR motifs; the kite signal is probably from the
   SSR period.
3. Look at `irreg_flag`:
   - `ok` — irregularity scan ran; check the indel metrics.
   - `too_short` — array_length / rank-1 period < 10; not enough copies.
   - `no_kite` — kite emitted no rank-1 peak; nothing to scan against.
4. If `kite_clusters` has multiple high-`total_score` clusters at
   non-harmonic spacings, drop into `rescore`.

### From `rescore.peaks.tsv` — peak-level diagnostics

For each rescored row:

1. **`identity_med ≥ 0.85` + narrow IQR** (≤ 0.05) → clean real period.
2. **`identity_med ≥ 0.85` + `phantom = true`** → kite saw a sub-period
   harmonic; the real period is `period + shift_med` bp.
3. **Moderate `identity_med` (0.5–0.7) + wide IQR + `subrepeat = true`**
   → localized subrepeat motif inside the founder monomer.
4. **`identity_med ≈ 0.3–0.4` + narrow IQR** → noise period; not real.
5. **High `coverage_frac` (≥ 0.85)** independently confirms a real
   period regardless of the flags.

### Cross-pipeline reading

A typical sequence on a new FASTA:

```bash
# whole-array snapshot
kitehor report sample.fa -o sample
less sample.report.tsv

# if the report is interesting, get peak-level detail
kitehor kite-periodicity sample.fa -o /tmp/k.tsv
kitehor rescore sample.fa --peaks /tmp/k.tsv.peaks.tsv -o sample
less sample.peaks.tsv
```

The `kite_peaks` / `kite_clusters` columns in the report tell you which
periods to focus on; rescore tells you which of those are real,
phantom, or subrepeat artifacts.

## Source layout

| Code | Purpose |
|---|---|
| `src/rescore/mod.rs` | `Config`, `PhantomConfig`, `SubrepeatConfig`, `run_subcommand`, `rescore_one`, `enforce_subrepeat_founder_gate`, unit tests |
| `src/rescore/aligner.rs` | banded semi-global edit distance kernel; `Scratch` + `ScoringConfig` |
| `src/rescore/sample.rs` | per-record deterministic anchor sampler with N-rejection |
| `src/rescore/io.rs` | peaks-TSV reader (column-preserving) |
| `src/report/mod.rs` | `Config`, `ReportRow`, `build_row`, `run_subcommand` |
| `src/report/io.rs` | report TSV writer |
| `tests/rescore_smoke.rs` | integration test for rescore |

Internal stages reused by report:

| Library | Used by report for |
|---|---|
| `src/kite.rs` | rank-1 period + all kept peaks |
| `src/rule_classify/cluster.rs` | period clustering at `±cluster_tol` |
| `src/ssr/scan.rs` | SSR motif coverage |
| `src/irregularity/mod.rs` | indel-event metrics |

## Known limitations

### rescore

- Detection floor at default K=200 is ≈ 10 % array footprint for the
  subrepeat flag. Smaller footprints need `--samples 500` or higher.
- `--max-period 5000` default skips very long candidates. For
  centromeric satellite arrays with HOR units > 5 kb, set `--max-period
  0` (unlimited) and expect ~3× wall time.
- Subrepeat false positives on true periods: ~2 % on synthetic ground
  truth, mostly when the founder gate's `identity_med ≥ 0.70` threshold
  is borderline (e.g. divergent satellites).
- The flag `coverage_frac` is computed but not currently used for any
  decision; it's a standalone diagnostic for downstream consumers.

### report

- Depends on rank-1 kite peak for SSR and irregularity stage input. When
  rank-1 is a phantom or noise peak (~1 % of cases on IPIP), the
  downstream metrics are computed against the wrong period.
- `irreg_flag = no_kite` rows have all `irreg_*` columns empty;
  downstream filters should treat this as "no signal", not "passed".
- No tandem-validate output by design — pair with `analyze` if you need
  the localization signal.

## Where to read next

- [`docs/rescore.md`](rescore.md) — full CLI reference, algorithm
  detail, calibration on synthetic ground truth, worked examples.
- [`docs/report.md`](report.md) — column-by-column semantics for the
  20-column TSV, irregularity flag taxonomy, design rationale for
  excluding tandem-validate.
- [`docs/rule_proto.md`](rule_proto.md) — the upstream rule-classify
  pipeline that `analyze` runs and `report` deliberately avoids.
- [`docs/irregularity_and_subrepeat_v0_12.md`](irregularity_and_subrepeat_v0_12.md)
  — the irregularity algorithm used by `report` columns 13–20.
