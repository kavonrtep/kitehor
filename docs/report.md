# `kitehor report` — observation-only output

A defensive-analysis mode that runs kite, peak clustering, SSR scan,
and irregularity scan against an input FASTA and emits a single
tab-separated table with **measurements only** — no categorical
verdicts, no `combined_class`, no rule-classify HOR/simple_tr/
unresolved labels. Designed as a clean-slate alternative to
`analyze` + `summary-merge` when the cascade's classifications are
not trustworthy or you want raw numbers to filter on yourself.

Designed as a sibling to `analyze`. May replace it later.

> **Onboarding**: read [`docs/onboarding_pipelines.md`](onboarding_pipelines.md)
> for a side-by-side view of `report` and `rescore`. This document is
> the per-column reference.
>
> **Related**: for **per-peak** identity diagnostics (rather than the
> whole-array snapshot this stage produces), see
> [`docs/rescore.md`](rescore.md). A typical workflow uses `report`
> for the array-level triage and `rescore` to drill into specific
> peaks of interest.

## Status

Stable. The 20-column schema and irregularity flag taxonomy have not
changed since the v0.12 release; the only addition since then was
`irreg_dropout_rate_per_pair` (column 19). Field-calibrated against
dotplots; no automated calibration corpus.

Per-record wall time on the IPIP 2026-04-14 corpus (3024 records,
305 MB): ~20 ms median, kite-dominated.

## Usage

```bash
kitehor report <fasta> -o <prefix>
```

Writes `<prefix>.report.tsv`. Tunables:

| flag | default | what |
|---|---|---|
| `--cluster-tol` | 0.015 | peak-clustering relative-period tolerance |
| `--irreg-step-min-frac-of-p` | 0.05 | irregularity event step floor (fraction of P) |
| `--irreg-min-copies-for-scan` | 10 | irregularity `too_short` cutoff (array_length / P) |
| `--threads` | 0 | rayon worker count (0 = auto) |

Reads FASTA only. Internal stages: `kite::analyze` →
`rule_classify::cluster_peaks` → `ssr::scan::scan_record` →
`irregularity::analyse_record`. Records are processed in parallel via
rayon.

## Output schema

20 tab-separated columns, one row per FASTA record. Lists inside a
cell use `;` between entries and `:` between fields within an entry.
Empty cells (not `NA`) indicate missing / non-applicable values.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier (everything after `>` up to first whitespace) |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `kite_n_peaks` | int | number of peaks kite kept (after its own filtering) |
| 4 | `kite_peaks` | list-str | `period:score2_norm` entries, `;`-separated, sorted by score desc — **all** kept peaks |
| 5 | `kite_n_clusters` | int | number of period clusters after `cluster_peaks(tol=0.015)` |
| 6 | `kite_clusters` | list-str | `rep_period:total_score:n_peaks` entries, `;`-separated, sorted by total_score desc — **all** clusters |
| 7 | `ssr_total_coverage_pct` | float | sum of all SSR motifs' coverage, as percent of `array_length` |
| 8 | `ssr_dominant_motif` | str | top SSR motif by coverage (upper canonical, e.g. `AT`) |
| 9 | `ssr_dominant_motif_length` | int | length in bp of the dominant motif |
| 10 | `ssr_dominant_motif_repeats` | int | total dominant-motif repeat units summed across the array |
| 11 | `ssr_dominant_coverage_pct` | float | dominant motif's coverage, as percent of `array_length` |
| 12 | `ssr_top_motifs` | list-str | top-3 motifs by coverage: `motif:pct;motif:pct;motif:pct` (percent of array_length) |
| 13 | `irreg_flag` | str | scan outcome (see [Irregularity flags](#irregularity-flags) below) |
| 14 | `irreg_n_kmer_groups` | int | independent phase-bin groups used by the scan (0 when `flag != ok`) |
| 15 | `irreg_indel_event_count` | int | discrete indel events detected; empty when `flag != ok` |
| 16 | `irreg_indel_burden_pct` | float | `(Σ |step_size|) / array_length × 100` |
| 17 | `irreg_indel_max_shift_bp` | float | size of the largest detected event |
| 18 | `irreg_indel_drift_bp_per_kb` | float | net cumulative drift in bp per kb |
| 19 | `irreg_dropout_rate_per_pair` | float | fraction of expected k-mer occurrences missing (substitution / divergence signal — orthogonal to indels) |
| 20 | `irreg_baseline_jitter_bp` | float | per-array MAD-based residual scale |

### Column 4 — `kite_peaks`

The raw output of kite's peak detector after its own filtering. Each
entry:

```
period:score2_norm
```

- `period`: integer bp.
- `score2_norm`: kite's normalised `score · log2(period)` value. Sums
  to 1 across all kept peaks within one record. Used by rule-classify
  to rank candidate tile/founder peaks.

Sorted descending by `score2_norm`. ALL peaks are emitted — no top-N
truncation.

Example (TRC_1:chr1_109549_122053):

```
50:0.359798;100:0.087102;150:0.035368;200:0.026694;1037:0.018104;…
```

The 50-bp founder dominates, with harmonic ladder at 100, 150, 200,
and a separate signal at ~1037 bp.

### Column 6 — `kite_clusters`

Output of `rule_classify::cluster::cluster_peaks` at
`tol = --cluster-tol` (default 0.015 = ±1.5%). Single-linkage
clusters peaks whose consecutive periods are within `tol` relative.
Each entry:

```
rep_period:total_score:n_peaks
```

- `rep_period`: cluster's weighted-mean period (score-weighted across
  member peaks), 2 decimals.
- `total_score`: sum of `score2_norm` across member peaks. This is
  the aggregated cluster score that rule-classify uses for tile /
  founder candidate ranking.
- `n_peaks`: how many raw peaks the cluster contains. `1` for a
  singleton; `≥2` indicates the clusterer merged adjacent peaks.

Sorted descending by `total_score`. ALL clusters emitted.

Example:

```
50:0.359798:1;100:0.087102:1;1186.67:0.015:2;…
```

The third cluster shows that two raw peaks near period 1186 were
merged into a single cluster centred at 1186.67 bp.

### Columns 7–12 — SSR

All percentages are relative to `array_length` in bp.

| field | scope |
|---|---|
| `ssr_total_coverage_pct` | sum of per-motif coverage (after merging overlapping motif regions; capped at 100%) |
| `ssr_dominant_motif` / `_length` / `_repeats` / `_coverage_pct` | the single motif with the most total coverage |
| `ssr_top_motifs` | up to 3 motifs, sorted by coverage desc, format `motif:pct;…` |

**Note on the top-3 percentages**: per-motif coverages are NOT capped
at 100% individually because the same bp can be counted under
multiple overlapping canonical motifs (e.g. an `(AT)n` region is
also `(ATAT)n`). Only the aggregated `ssr_total_coverage_pct` is
capped at 100%.

`ssr_dominant_motif` and `_motif_length` are empty when no SSR
survived the per-motif-length minimum-repeats filter. The empty cell
is preferred over `NA` so column types stay numeric where possible.

### Columns 13–20 — Irregularity

The irregularity scan quantifies indel-like phase-shift events using
register-locked k-mers as phase markers (see
`docs/irregularity_and_subrepeat_v0_12.md` §3 for the full
algorithm). Two metric families:

- **Indel signal** (`indel_event_count`, `indel_burden_pct`,
  `indel_max_shift_bp`, `indel_drift_bp_per_kb`): coordinate shifts
  detected across multiple independent k-mer groups. A "real" indel
  pushes many register markers together.
- **Dropout signal** (`dropout_rate_per_pair`): fraction of expected
  consecutive-pair k-mer observations that were missing (i.e.
  distance = 2·P/3·P/… instead of P, indicating one or more
  occurrences were knocked out by substitution). Orthogonal to indel
  signal.

### Irregularity flags

| flag | meaning | metric cols filled? |
|---|---|---|
| `ok` | scan completed; all metric columns valid | yes |
| `no_period` | kite produced no peaks for this record | no |
| `too_short` | `array_length < min_copies_for_scan × P` (default 10×) | no |
| `no_register_lock` | < 3 register-locked k-mer groups at the kite period | no |
| `too_long` | post-aggregation signal would exceed 50,000 copy indices | no |

`no_register_lock` is the most common non-`ok` flag in real data —
arrays whose top kite period doesn't have phase-locked k-mer markers
compatible with `m·P` or `P/m` for small integer `m`. Long-period
arrays where kite picks a sub-monomer scale typically land here.

## Cell encoding notes

- **`;` separates list entries**, **`:` separates fields within an entry**.
  Both delimiters were chosen to be distinct from tab.
- **Floats**: trailing zeros trimmed (e.g. `0.5` not `0.500000`).
  Trailing decimal point also trimmed (e.g. `100` not `100.`).
- **Empty cells** are emitted as the empty string, never the literal
  `NA`. Distinguishes "no signal" (empty) from "value is the string
  NA" (which never occurs).
- **Negative zero** (`-0`) may appear in irregularity drift columns
  from float arithmetic. Treat as 0.

## Worked example

Smoke fixture's `hor_k3` record (24 kb HOR array, monomer = 100 bp,
tile = 300 bp, 80 tile copies):

```
record_id                      hor_k3
array_length                   24000
kite_n_peaks                   55
kite_peaks                     300:0.357;100:0.129;600:0.099;200:0.098;900:0.026;1200:0.010;…
kite_n_clusters                52
kite_clusters                  300:0.357:1;100:0.129:1;600:0.099:1;200:0.098:1;900:0.026:1;…
ssr_total_coverage_pct         0
ssr_dominant_motif             (empty)
ssr_top_motifs                 (empty)
irreg_flag                     ok
irreg_n_kmer_groups            31
irreg_indel_event_count        0
irreg_indel_burden_pct         0
irreg_indel_max_shift_bp       0
irreg_indel_drift_bp_per_kb    0
irreg_dropout_rate_per_pair    0.611
irreg_baseline_jitter_bp       1
```

Reading the row:
- Kite found 55 peaks, top-3 are the tile (300 bp), founder (100 bp),
  and 2× tile (600 bp).
- After clustering at ±1.5%, 52 distinct period clusters remain.
- No SSR content (clean ACGT-balanced array).
- Irregularity: 31 phase-bin groups, 0 detected indel events
  (expected — clean HOR), baseline jitter 1 bp.
- `dropout_rate_per_pair = 0.611` means 61% of expected k-mer
  consecutive pairs were "missing one" — typical for HOR arrays
  where many k-mers occur only in some founder copies due to
  inter-founder divergence.

## What `report` deliberately does NOT include

- No `combined_class` / `hor_verdict` / `tv_decision` columns
  (those are labels, not measurements).
- No `tandem-validate` columns — tandem-validate's `tv_density` /
  `tv_decision` are inputs to the v0.12 cascade but were excluded
  from the report by design (Petr's spec: keep only measurements
  that don't depend on cascade-style judgments).
- No `confidence` / `share` numbers from rule-classify (compound of
  the cascade's biases).

If you need the cascade-style verdict alongside the values, run
`kitehor analyze` separately — both can be run on the same FASTA;
they don't conflict.

## Performance

Per-record cost is dominated by kite. The full 3024-record IPIP
corpus (305 MB FASTA) processes in ~1m 45s on 8 threads (~28 records/
sec per core after the parallel slack). Memory usage stays under 1 GB
RSS — kite's k-mer histogram is the largest transient allocation.

## Source

- `src/report/mod.rs` — driver + cell encoders
- `src/report/io.rs` — TSV writer + kite-summary reader
- CLI: `src/cli.rs::ReportArgs`, dispatch in `src/main.rs`

Currently NOT wired into `kitehor analyze`. The two commands are
independent — `analyze` is the cascade-style pipeline, `report` is
the measurements-only sibling.
