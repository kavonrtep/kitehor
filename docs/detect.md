# `kitehor detect` / `detect-batch` — v2 line-width detector

A sequence-agnostic threshold-rule classifier that operates on **row
embeddings** of the input FASTA (rather than the k-mer periodograms
that drive the `rule-classify` cascade). Reads period candidates per
array, expands each to a 2-D matrix at the candidate width, and
extracts a fixed feature bundle: column conservation, row
autocorrelation, edge field, shift signal, irregularity score, and a
consensus monomer / HOR unit.

Per-class accuracy on the 1600-case `ground_truth_v2/` benchmark is
**94.4 %** (M6 acceptance gate, alt-aware). See
[`docs/new/detect_spec.md`](new/detect_spec.md) for the full algorithm
design and [`docs/new/detect_impl_plan.md`](new/detect_impl_plan.md)
for the implementation contract.

## Status

Stable. Five output files, schemas frozen at M0–M7.3:

- `<prefix>.properties.tsv` — primary per-array result (20 cols)
- `<prefix>.segments.tsv` — per-segment breakdown (13 cols; M7.2)
- `<prefix>.width_features.tsv` — per-tested-width diagnostic (17 cols)
- `<prefix>.consensus.fa` — per-array monomer + optional `hor_unit`
- `<prefix>.diagnostics.json` — structured diagnostics (`schema_version=2`, M7.3)

Auto-mode (`detect` without `--periods`) additionally writes
`<prefix>.periods.tsv`.

## Usage

```bash
# Auto-mode — kite-periodicity + rule classifier run internally,
# their period output is persisted and consumed by detect.
kitehor detect <fasta> -o <prefix>

# With explicit period candidates (e.g. from a tuned kite run)
kitehor detect <fasta> --periods <periods.tsv> -o <prefix>
                       --allow-missing-periods  # if some FASTA records have no period rows

# Batch over a directory of <stem>.fa + <stem>.periods.tsv pairs
kitehor detect-batch --fasta-dir <dir> --periods-dir <dir> --out-dir <dir>
```

### Auto-mode behaviour

Without `--periods`, `detect`:

1. Runs `kite-periodicity --classify` with kite defaults.
2. Maps each classifier verdict to up to 3 period rows
   (founder = 0.95, tile = 0.90, other top-3 = 0.60; `Unresolved`
   hints at 0.50 / 0.40 / 0.30).
3. Persists those rows to `<prefix>.periods.tsv`.
4. Implicitly applies `--allow-missing-periods` so QC-rejected
   records end up classified as `ambiguous` instead of erroring.

To use tuned kite parameters, run kite explicitly and feed the
resulting periods via `--periods`.

### Key flags

| flag | what |
|---|---|
| `--periods <path>` | period-candidate TSV (`array_id`, `period_bp`, `period_score`, optional `source`); when omitted, auto-mode runs |
| `--config <path>` | TOML override of `DetectorConfig` defaults |
| `--viz-dir <dir>` | per-array matrix/signal export root |
| `--export-raster` / `--export-shift` / `--export-edges` / `--export-ic` | per-channel visualisation toggles (require `--viz-dir`) |
| `--allow-missing-periods` | downgrade "no rows for record" errors to warnings |
| `--allow-extra-periods` | downgrade "periods TSV has rows for unknown record IDs" errors to warnings |
| `--threads <N>` | rayon worker count (0 = auto) |

## Output schemas

### `<prefix>.properties.tsv` (20 columns)

One row per FASTA record. The frozen schema is in
`src/detect/types.rs::PROPERTIES_HEADER`.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | FASTA record identifier |
| 2 | `length_bp` | int | sequence length in bp |
| 3 | `class` | str | summary class — `HOR`, `irregular_HOR`, `simple_TR`, `mixed`, `ambiguous` |
| 4 | `base_width_bp` | int / `NA` | inferred base monomer width |
| 5 | `hor_k` | int / `NA` | inferred HOR multiplicity (when `class` ∈ {HOR, irregular_HOR}) |
| 6 | `hor_length_bp` | int / `NA` | inferred HOR unit length (`hor_k · base_width_bp`) |
| 7 | `n_complete_copies` | int / `NA` | number of complete HOR-unit (or monomer) copies in the array |
| 8 | `column_conservation` | float / `NA` | background-corrected mean column IC at the chosen width (regime A indicator) |
| 9 | `phase_separation` | float / `NA` | best-lag autocorrelation gap — separates HOR slot phases (regime B) |
| 10 | `mean_shift_bp` | float / `NA` | mean per-row alignment shift from the natural width (drift detector) |
| 11 | `wobble_amplitude_bp` | float / `NA` | sinusoidal envelope amplitude of the shift signal |
| 12 | `wobble_periodicity_bp` | float / `NA` | dominant period of the shift signal (bp), when sinusoidal |
| 13 | `n_phase_shifts` | int | discrete phase-shift breakpoints detected |
| 14 | `phase_shift_positions` | list-str / `NA` | breakpoint positions (bp), comma-separated |
| 15 | `phase_shift_offsets` | list-str / `NA` | per-breakpoint shift offsets (bp, signed), comma-separated |
| 16 | `irregularity_score` | float / `NA` | scalar wobble + indel severity (0 = clean) |
| 17 | `inter_monomer_identity` | float / `NA` | **approximation** — `R(1)` at base width (k-mer-composition row similarity, not pairwise sequence identity). Useful as a regime A vs B/C indicator; not a calibrated identity |
| 18 | `confidence` | float / `NA` | heuristic confidence score (not a calibrated probability) — see [`docs/new/detect_impl_plan.md`](new/detect_impl_plan.md) §0 A6 |
| 19 | `n_segments` | int | number of segments emitted (1 if no segmentation triggered) |
| 20 | `reason` | str | free-text decision-path tag |

### `<prefix>.segments.tsv` (13 columns)

One row per segment, only when `n_segments > 1` or the same-width mixed
override fires. `class` is the whole-array final class applied to the
sub-range (light per-segment override only; heavy per-segment
classification is deferred to M8+).

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | FASTA record identifier |
| 2 | `segment_id` | int | 0-based ordinal within the record |
| 3 | `start_bp` | int | inclusive start in bp |
| 4 | `end_bp` | int | exclusive end in bp |
| 5 | `class` | str | class label applied to this sub-range |
| 6 | `base_width_bp` | int / `NA` | base width for this segment |
| 7 | `hor_k` | int / `NA` | HOR multiplicity for this segment |
| 8 | `column_conservation` | float / `NA` | column conservation in this segment |
| 9 | `phase_separation` | float / `NA` | phase separation in this segment |
| 10 | `wobble_amplitude_bp` | float / `NA` | wobble amplitude in this segment |
| 11 | `irregularity_score` | float / `NA` | irregularity score in this segment |
| 12 | `consensus_identity_to_reference` | float / `NA` | Hamming identity (N-skip) of this segment's consensus to the medoid block's consensus (M7.2; filled for same-width mixed-override segments only — phase-shift-derived segments leave it `NA`) |
| 13 | `consensus_identity_coverage` | float / `NA` | fraction of consensus length used in the identity comparison (positions where neither side was N) |

### `<prefix>.width_features.tsv` (17 columns)

One row per `(array, tested-width)` pair. Diagnostic dump driving the
M2–M3 features.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | FASTA record identifier |
| 2 | `width_bp` | int | tested width in bp |
| 3 | `rows` | int | number of complete rows at this width |
| 4 | `column_IC` | float / `NA` | background-corrected mean column information content |
| 5 | `fraction_conserved_columns` | float / `NA` | fraction of columns with IC ≥ threshold |
| 6 | `row_lag1_similarity` | float / `NA` | row-to-row autocorrelation at lag 1 (`R(1)`) |
| 7 | `best_lag` | int / `NA` | best autocorrelation lag (in rows) |
| 8 | `best_lag_score` | float / `NA` | autocorrelation value at `best_lag` |
| 9 | `phase_separation` | float / `NA` | best_lag_score − R(1) (HOR slot phase contrast) |
| 10 | `vertical_edge_rate` | float / `NA` | density of column-edge transitions |
| 11 | `column_edge_autocorr_k` | int / `NA` | best lag (rows) of the column-edge-rate autocorrelation |
| 12 | `column_edge_autocorr_score` | float / `NA` | column-edge autocorrelation value at `_k` |
| 13 | `mean_shift_bp` | float / `NA` | mean per-row shift from natural width |
| 14 | `wobble_amplitude_bp` | float / `NA` | sinusoidal envelope amplitude of the shift signal |
| 15 | `n_phase_shifts` | int | discrete phase-shift breakpoints detected at this width |
| 16 | `irregularity_score` | float / `NA` | wobble + indel severity at this width |
| 17 | `class_hint` | str | per-width class hint (`simple_TR_base`, `HOR_base`, `HOR_unit`, `unsupported`) |

### `<prefix>.consensus.fa`

One or more FASTA records per array:

```
>{array_id}.monomer  length={base_width_bp}
ACGT...
>{array_id}.hor_unit  length={hor_length_bp}  k={hor_k}
ACGT...
```

- `monomer` is always present for arrays with a resolved
  `base_width_bp`. For `class = mixed`, multiple `monomer.{segment_id}`
  records may be emitted (one per analysis block).
- `hor_unit` is present only when `class ∈ {HOR, irregular_HOR}` and a
  valid `hor_k` was inferred. Built directly from the HOR-unit width's
  matrix, not by repeating the monomer consensus.

### `<prefix>.periods.tsv` (auto-mode only)

The period rows the detector consumed, persisted for reproducibility.
4 columns, schema shared with `kitehor synth`:

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | FASTA record identifier |
| 2 | `period_bp` | int | candidate period |
| 3 | `period_score` | float | candidate score (rank-derived: founder 0.95, tile 0.90, …) |
| 4 | `source` | str | origin tag (`kite_classifier_founder`, `kite_classifier_tile`, `kite_top3`, `unresolved_hint`) |

### `<prefix>.diagnostics.json`

Structured per-array dump intended for programmatic consumption and
plotting. `schema_version = 2` since M7.3. Top-level keys:

- `array_id`, `length_bp`, `class`, `reason`, `confidence`
- `properties` — array-level Properties record (same fields as
  `properties.tsv`)
- `segments` — list of segment records
- `width_features` — list of `(width_bp, WidthFeatures)` entries
- `consensus` — list of `ConsensusRecord` enums (monomer + optional
  hor_unit + optional per-segment monomers)
- `detailed` — internal diagnostic block (R(k) curves, shift signal
  samples, edge-rate samples, breakpoint signals)

### Visualisation (when `--viz-dir` is set)

Per-array TSV + (optional) PNG matrix dumps. Channels:

- `--export-raster` → row-major one-hot matrices (PNG only)
- `--export-shift` → per-row shift signal (TSV + PNG)
- `--export-edges` → per-column edge rates (TSV + PNG)
- `--export-ic` → per-column information content (TSV + PNG)

File naming: `<viz_dir>/<array_id>__<channel>__width=<w>.{tsv,png}`.
PNG output requires the `viz` Cargo feature (on by default).

## Tradeoffs and limits

- Detection is **width-driven**: at least one period candidate must
  land in the detector's expand window (`[period_bp · 0.5,
  period_bp · 1.5]` by default). Outside this window, the candidate
  is silently dropped.
- `inter_monomer_identity` is an approximation (see column 17 above).
  Treat as a regime indicator, not a calibrated identity.
- `confidence` is a heuristic, not a calibrated probability. Use the
  class label + the underlying features (cols 8–16) for downstream
  decisions.

## Source

- `src/detect/types.rs` — `Properties`, `Segment`, `WidthFeatures` definitions + frozen headers
- `src/detect/mod.rs` — driver (`run_one`)
- `src/detect/io.rs` — TSV writers + period reader
- `src/detect/classify.rs` — array-level decision tree
- `src/detect/segment.rs` — segmentation
- `src/detect/consensus.rs` — column-vote monomer + HOR unit consensus
- `src/detect/viz.rs` — visualisation dumps
