# `kitehor irregularity` — indel-event scan

Per-record scan that quantifies indel-like phase-shift events inside a
tandem array using **register-locked k-mer phase markers**. Reports
discrete indel events, indel burden, drift, and a substitution-only
"dropout" rate orthogonal to the indel signal. Port of
`tools/rule_proto/irregularity_v2.py` (Approach 6 — distance-residual
+ phase-bin clustering).

The full algorithm is documented in
[`docs/irregularity_and_subrepeat_v0_12.md`](irregularity_and_subrepeat_v0_12.md)
§§ 2–3. This page is the per-column / per-flag reference for the CLI
subcommand.

## Status

Stable. v0.12 schema (14 columns), no changes since the irregularity-v2
Rust port landed.

## Usage

```bash
kitehor irregularity <fasta> --kite <kite_summary.tsv> -o <prefix>
```

Writes `<prefix>.irregularity.tsv`. `--kite` supplies the per-record
top period (`monomer_size`) that anchors the scan — typically the
output of `kitehor kite-periodicity` (no `--classify` needed).

### Flags

| flag | default | what |
|---|---|---|
| `--k <K>` | `6` | k-mer length (must match the kite run) |
| `--top-kmers <N>` | `100` | top-N most frequent k-mers per record for the phase-bin scan |
| `--min-copies-for-scan <N>` | `10` | minimum `array_length / period` ratio; below → `too_short` |
| `--step-min-frac-of-p <F>` | `0.05` | event step floor (fraction of `P`) |
| `--min-kmer-groups <N>` | `3` | minimum independent phase-bin groups required to run; below → `no_register_lock` |

## Output schema — `<prefix>.irregularity.tsv` (14 columns)

One row per FASTA record. Cells use `NA` (not empty) for missing
values; the `notes` column carries free-text diagnostics that the flag
alone doesn't convey.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `period_P` | int / `NA` | kite-supplied period the scan ran against |
| 4 | `n_kmer_groups` | int | independent phase-bin groups the scan used |
| 5 | `n_pairs_total` | int | total consecutive-pair k-mer observations across groups |
| 6 | `baseline_jitter_bp` | float / `NA` | per-array MAD-based residual scale (bp); below-jitter events are filtered out |
| 7 | `indel_event_count` | int / `NA` | discrete indel events detected (clustered phase shifts) |
| 8 | `indel_burden_pct` | float / `NA` | `(Σ \|step_size\|) / array_length × 100` |
| 9 | `indel_max_shift_bp` | float / `NA` | size of the largest detected event (bp) |
| 10 | `indel_drift_bp_per_kb` | float / `NA` | net cumulative drift in bp per kb of array |
| 11 | `dropout_event_count` | int / `NA` | sum of `n_skipped` across all consecutive pairs (positions where the expected k-mer occurrence was missing — orthogonal to indels) |
| 12 | `dropout_rate_per_pair` | float / `NA` | `dropout_event_count / n_pairs_total` — fraction of expected pairs that dropped an intermediate occurrence |
| 13 | `flag` | str | scan outcome (see below) |
| 14 | `notes` | str | free-text diagnostic — e.g. "no kite top period", "k-mer groups < min" |

### Flag values

| flag | meaning | metric cols filled |
|---|---|---|
| `ok` | scan completed | columns 6–12 valid |
| `no_period` | kite produced no top period for the record | `NA` for metric cols |
| `too_short` | `array_length < min_copies_for_scan · period_P` | `NA` for metric cols |
| `no_register_lock` | fewer than `--min-kmer-groups` phase-locked k-mer groups at the kite period — typical when kite picked a sub-monomer scale | `NA` for metric cols |
| `too_long` | post-aggregation signal would exceed 50,000 copy indices (safety guard) | `NA` for metric cols |

### Indel vs dropout — what the split means

- **Indel signal** (`indel_event_count`, `indel_burden_pct`,
  `indel_max_shift_bp`, `indel_drift_bp_per_kb`): coordinate shifts
  detected across multiple independent k-mer groups simultaneously.
  A "real" indel pushes many register markers together by the same bp
  offset.
- **Dropout signal** (`dropout_rate_per_pair`): fraction of expected
  consecutive-pair k-mer observations that were missing (distances of
  `2·P`, `3·P`, ... instead of `P`, indicating intermediate
  occurrences were knocked out by substitution / divergence). This is
  the substitution-rate proxy that does **not** disturb the phase
  register.

The split is the load-bearing v2 fix: the two metrics are orthogonal,
so a divergent-but-non-indel array shows up as high `dropout_rate`
with `indel_event_count = 0`, and vice versa.

## Used downstream by

- `kitehor report` — columns 13–20 of `report.tsv` (with the
  `report.md` column naming convention; see
  [`docs/report.md`](report.md)).
- `kitehor analyze` — does **not** invoke this stage. Run it
  separately when you need the metrics.

## Performance

Per-record cost is dominated by the phase-bin scan over the top-N
k-mers. On a 30 kb array with `--k 6 --top-kmers 100` the scan runs in
≈ 5–10 ms. Records are processed in parallel via rayon.

## Source

- `src/irregularity/mod.rs` — driver
- `src/irregularity/scan.rs` — phase-bin clustering + event detection
- `src/irregularity/io.rs` — TSV writer + kite-summary reader
- CLI: `src/cli.rs::IrregularityArgs`, dispatch in `src/main.rs`
