# `scan` port plan — shifted-self-alignment scan into rescore

Implementation contract for porting the Python prototype
[`tools/subrepeat_scan/scan.py`](../../tools/subrepeat_scan/scan.py)
into Rust as part of the `kitehor rescore` pipeline. After this port
the scan becomes the **primary subrepeat-detection workhorse** —
the alignment-based bimodality flags stay observational; downstream
decisions can either be gated on the scan's `scan_occupancy_frac`
directly or on a combination (decision deferred until after we see
real-data calibration).

## 0. Scope and non-goals

**In scope:**
- Per-peak scan computed for **every kite-reported `(record,
  period)` row that gets rescored** — no separate period filter.
  The input is already kite peaks (kite vetted them as real
  periodic signals), so the existing rescore filters
  (`--top-n`, `--min-period`, `--max-period`) define the scan's
  scope. The plan retains a single tunable, the identity
  threshold.
- Two new TSV columns emitted as observational data — the existing
  `subrepeat` flag is not modified in this PR.
- CLI flag + main wire-up + docs.

### Note on long-period peaks

When the scan runs at a long candidate period (e.g. `period =
founder_period` for a clean tandem array): every position in the
array satisfies `seq[i] == seq[i+period]`, so the scan reports
`occupancy_frac ≈ 1.0`. **This is correct**: the array *is* a
tandem at that period. The existing `subrepeat` flag's
period/founder ratio gate already prevents these rows from firing
`subrepeat=true`. So `scan_occupancy_frac` carries useful signal
at both scales:

- short period rows: `scan_occupancy_frac` = "is there a nested
  short tandem at this period inside this array?"
- long period rows: `scan_occupancy_frac` ≈ 1.0 confirms the
  founder-scale tandem (downstream consumers can subtract the
  expected baseline if interested in nested-only signal).

**Out of scope (future PRs):**
- Promoting `scan_occupancy_frac ≥ threshold` to a gate on the
  `subrepeat` flag. Will land after we calibrate against a labelled
  IPIP subset.
- Period-aliasing resolution across rows of the same record (each
  row currently reports its own period's scan independently).
- Banded local self-alignment refinement of interval boundaries
  (the prototype skips this and shows the simple per-base scan is
  already informative).
- Reporting per-interval coordinates in the TSV (would change the
  schema from "one row per kite peak" to "one row per interval";
  out of scope for this integration).

## 1. Algorithm — what the scan computes

For one `(seq, period)` pair:

1. Build `match[i] = 1 if seq[i] == seq[i+period] else 0` for
   `i ∈ [0, L − period)` — a `u8`/`bool` array of length `L − P`.
   N bases are encoded so they never match themselves (reuse the
   existing convention in `kmer_scan::encode_base`).
2. Compute the windowed match rate
   `rate[i] = mean(match[i .. i+period])` for
   `i ∈ [0, L − 2·period]`. Use the cumulative-sum trick —
   `O(L)` not `O(L · P)`.
3. Find every maximal contiguous run of indices where
   `rate[i] ≥ id_threshold`, keeping only runs whose length is
   ≥ `(min_copies − 1) · period`.
4. For each surviving run `[s, e)` (over window-start indices),
   the corresponding sequence-coordinate interval is
   `[s, min(L, e + period − 1)]`. Accumulate
   `occupied_bp = Σ (interval_end − interval_start)` taking the
   union to avoid double-counting overlaps within the same
   period scan.
5. Return:
   - `n_intervals`: int — surviving runs
   - `occupied_bp`: int — union of interval lengths
   - `occupancy_frac`: f64 ∈ [0, 1] — `occupied_bp / L`

Per-row cost: `O(L)`. The K-mer position cache built for the
existing autoF / phaseC diagnostics is **not** used here; the scan
operates directly on the raw byte sequence.

## 2. Module layout

New file: `src/rescore/scan.rs` mirroring `kmer_scan.rs`'s
structure.

```rust
// src/rescore/scan.rs

//! Shifted self-alignment scan for nested short tandem repeats.
//!
//! For one (sequence, period) pair, computes the windowed match
//! rate at lag `period` and identifies contiguous runs above an
//! identity threshold. Per-row diagnostic emitted from
//! `kitehor rescore` as `scan_occupancy_frac` and `scan_n_intervals`.

pub struct ScanResult {
    pub n_intervals: usize,
    pub occupied_bp: usize,
    pub occupancy_frac: f64,
}

pub fn scan_one_period(
    seq: &[u8],
    period: usize,
    id_threshold: f64,
    min_copies: usize,
) -> Option<ScanResult>;

// Internal helpers (private, but exposed in the test module):
fn build_match(seq: &[u8], period: usize) -> Vec<u8>;
fn windowed_rate(match_buf: &[u8], period: usize) -> Vec<f64>;
fn find_runs_above(rate: &[f64], threshold: f64, min_length: usize)
    -> Vec<(usize, usize)>;
fn union_length(intervals: &[(usize, usize)]) -> usize;
```

Wire in `src/rescore/mod.rs`:

```rust
pub mod aligner;
pub mod io;
pub mod kmer_scan;
pub mod sample;
pub mod scan;        // NEW
```

## 3. Data model

### `ScanConfig` in `src/rescore/mod.rs`

Mirrors the existing `KmerSpatialConfig` pattern.

```rust
#[derive(Debug, Clone, Copy)]
pub struct ScanConfig {
    /// Minimum per-window match rate for a window-start index to
    /// qualify as "tandem evidence" at this period. Calibrated on
    /// the synth corpus + IPIP TRC_104/TRC_666 to 0.55 (covers
    /// real divergence; rejects pure noise where match ≈ 0.25).
    pub id_threshold: f64,

    /// Minimum number of tandem copies. Translates to a minimum
    /// qualifying-window-run of `(min_copies − 1) · period`
    /// consecutive indices. Default 3 — matches the biological
    /// definition of a "nested short tandem".
    pub min_copies: usize,

    /// Skip the scan entirely when false; both columns emit `NA`.
    /// Default true. CLI: `--no-scan` flips this off.
    pub enabled: bool,
}

impl Default for ScanConfig {
    fn default() -> Self {
        Self {
            id_threshold: 0.55,
            min_copies: 3,
            enabled: true,
        }
    }
}
```

Add to `Config`:

```rust
pub struct Config {
    // ... existing fields ...
    pub scan: ScanConfig,           // NEW
}
```

### `RowStats` additions

Append two fields:

```rust
pub struct RowStats {
    // ... existing fields through kmer_phase_contrast ...

    /// Number of contiguous tandem-positive runs found at this row's
    /// period by `scan::scan_one_period`. `None` when the period is
    /// out of `ScanConfig::period_*` bounds, when scan is disabled,
    /// or when the array is too short for the kernel.
    pub scan_n_intervals: Option<usize>,

    /// Fraction of the array occupied by tandem-positive runs at
    /// this row's period — `occupied_bp / array_length` ∈ [0, 1].
    /// `None` under the same conditions as `scan_n_intervals`.
    pub scan_occupancy_frac: Option<f64>,
}
```

## 4. Pipeline integration

Two integration points in `src/rescore/mod.rs::run_subcommand`:

1. **Inside the existing parallel `rescore` loop** (right after
   `rescore_one` returns and the inline `kmer_spatial_contrast`
   call). The scan needs only the raw sequence and the row's
   period — no founder dependency, so it runs alongside the
   per-row work:

```rust
} else if let Some(record) = records.get(&row.case_id) {
    let band = cfg.resolved_band(row.period);
    let mut r = rescore_one(/* … */);

    if cfg.scan.enabled {
        if let Some(sr) = scan::scan_one_period(
            &record.seq,
            row.period,
            cfg.scan.id_threshold,
            cfg.scan.min_copies,
        ) {
            r.scan_n_intervals = Some(sr.n_intervals);
            r.scan_occupancy_frac = Some(sr.occupancy_frac);
        }
    }
    (r, true)
}
```

2. **The kmer-position post-pass** (where `kmer_autocorr_founder`
   and `kmer_phase_contrast` are filled after the founder gate) is
   not affected — the scan doesn't need a founder.

### No separate period filter

The scan reuses rescore's existing row-filtering machinery: rows
with `rank > top_n`, `period < min_period`, or `period >
max_period` already get `RowStats::na()`, which leaves both
scan columns as `None`. No additional period bound — the scan
runs on every kite period that rescore is already processing.

`scan_one_period` itself returns `None` (→ NA in the TSV) when
the array is too short to host any qualifying run:
`seq.len() < 2 · period + (min_copies − 1) · period`.

### Parallelism

The scan is per-row and stateless — slots straight into the
existing rayon `par_iter()` loop. No per-thread scratch
needed (a `Vec<u8>` of length `L − P` is allocated per call;
allocation cost is negligible vs the work).

## 5. TSV schema

Append two columns at the end of the existing 13:

```
… | spatial_contrast | founder_period | kmer_autocorr_founder | kmer_phase_contrast | scan_n_intervals | scan_occupancy_frac
```

New count: 15 appended columns. Existing column indices
unchanged.

- `scan_n_intervals` — int / `NA`. `0` means "scan ran but found
  no qualifying runs." `NA` means "scan did not run for this row"
  (out of period bounds, disabled, or too-short array).
- `scan_occupancy_frac` — `%.4f` ∈ `[0, 1]` / `NA`. Companion to
  `n_intervals`. `0.0000` and `NA` carry the same distinction as
  above.

## 6. CLI surface

In `src/cli.rs::RescoreArgs`:

```rust
/// Minimum windowed match rate for a window to qualify as
/// "tandem evidence" at the row's period. Default 0.55 (tuned
/// to admit real subrepeats at ~12 pct per-copy divergence;
/// rejects pure noise where match ≈ 0.25).
#[arg(long, default_value_t = 0.55)]
pub scan_id_threshold: f64,

/// Minimum number of tandem copies required for a positive call;
/// the minimum qualifying-window run is (min_copies − 1) · period.
/// Default 3 — matches the biological definition.
#[arg(long, default_value_t = 3)]
pub scan_min_copies: usize,

/// Disable the shifted-self-alignment scan entirely. Both scan
/// columns will emit NA for every row. Default: enabled.
#[arg(long, default_value_t = false)]
pub no_scan: bool,
```

Wire-up in `src/main.rs`:

```rust
scan: kitehor::rescore::ScanConfig {
    id_threshold: args.scan_id_threshold,
    min_copies: args.scan_min_copies,
    enabled: !args.no_scan,
},
```

## 7. Tests

### Unit tests in `src/rescore/scan.rs`

```rust
#[test] fn build_match_perfect_tandem_is_all_ones() {…}
#[test] fn build_match_random_is_random() {…}
#[test] fn windowed_rate_perfect_tandem_is_one() {…}
#[test] fn windowed_rate_random_is_near_quarter() {…}
#[test] fn find_runs_skips_too_short() {…}
#[test] fn scan_perfect_tandem_returns_full_occupancy() {…}
#[test] fn scan_no_tandem_returns_empty() {…}
#[test] fn scan_nested_subrepeat_returns_partial_occupancy() {…}
#[test] fn scan_with_noise_still_fires_at_threshold() {…}
#[test] fn scan_handles_short_array_gracefully() {…}
#[test] fn scan_at_period_minus_one_doesnt_panic() {…}
#[test] fn scan_n_handling_treats_each_n_as_mismatch() {…}
```

12 new tests targeted.

### Integration smoke

Update `tests/rescore_smoke.rs` header assertion to 15 appended
columns. Add an assertion on `tandem_pure` (60 bp tandem): expect
`scan_n_intervals = 1` and `scan_occupancy_frac > 0.95` at
period 60. (Wait — scan-period-max default is 45, so 60 is out of
range — the row will emit NA. Either drop the assertion or add a
second smoke fixture with a 36-bp tandem.)

Cleaner: add a `--scan-period-max 100` flag override in the test
so we can assert against the existing fixture without adding new
test data.

### Validation against Python prototype

Add `tests/scan_python_parity.rs` (`#[ignore]`-flagged) that:
1. Runs the Rust scan on every record of
   `test_data/subrepeat_sim/subrepeat_sim.fasta` at the same
   parameters as the Python prototype.
2. Runs the Python prototype on the same input via subprocess.
3. Asserts per-record `scan_occupancy_frac` matches within ±0.005.

Mirrors the existing `tests/tandem_validate_python_parity.rs`.

## 8. Validation strategy

### Synth corpus (`test_data/subrepeat_sim/`)

After implementation, re-run `kitehor rescore --top-n 0` on the
synth FASTA + truth TSV. Expected behaviour:

| case | truth | Python prototype | target Rust | tolerance |
|---|---|---|---|---|
| C01/C02/C10/C11 | 0 % | 0 % | 0 % | exact |
| C03 (borderline) | 0 % | 0 % | 0 % | exact |
| C04 (40 %) | undetectable with min_copies=3 | 0 % | 0 % | exact |
| C05 (60 %) | 58 % expected | 58 % | within ±2 % | smoothed-noise-tolerant |
| C06 (80 %, TRC_104-like) | 80 % | 69 % | within ±3 % | |
| C07 (degenerate, 100 %) | 100 % | 89 % | within ±5 % | |
| C08 (2 copies of 60) | undetectable | 0 % | 0 % | exact |
| C09 (HOR-nested) | partial | partial | match Python | within ±5 % |

### IPIP regression spot-checks

After integration, the rescore TSV on the full IPIP corpus must
reproduce these expected values:

- **TRC_104:chr3_411509737** P=36 (real nested subrepeat):
  `scan_occupancy_frac ≈ 0.15`, `n_intervals ≈ 16`.
- **TRC_104:chr3_411509737** P=180 (the founder row): scan runs
  at period 180; expect `scan_occupancy_frac ≈ 1.0` (the array
  IS a 180-bp tandem) — confirms the long-period behaviour
  described in §0.
- **TRC_115:chr7_353599568** P=1955 (near-founder FP): scan runs;
  expect `scan_occupancy_frac ≈ 0.0` because at P=1955 against
  a 2018-bp founder, `seq[i] == seq[i+1955]` is rare (positions
  are off-register).
- **TRC_115** P=2018 (the true founder): scan runs; expect
  `≈ 1.0`.
- **TRC_666:chr4_523278767** P=36 (real, high divergence): at
  threshold 0.55, `occupancy_frac ≈ 0.46`, `n_intervals ≈ 150`.
- **TRC_14:chr1_315645785** P=163 (rescore=true but no founder,
  SSR-rich content): scan runs; expect
  `scan_occupancy_frac ≈ 0.0` (no qualifying tandem run at
  P=163). This is the "rescore over-flagged" case the scan
  correctly rejects.

These are not unit-test assertions (they require the IPIP FASTA);
they're documented as the manual smoke check in `docs/rescore.md`.

## 9. Documentation updates

1. `README.md` — bump `rescore` row's "13 columns appended" → 15.
2. `docs/rescore.md`:
   - Update the headline column list (13 → 15).
   - Add per-column descriptions for `scan_n_intervals` +
     `scan_occupancy_frac` in the "Output schema" section.
   - Add a new sub-section **"Interpreting `scan_occupancy_frac`"**
     under the existing "Interpreting `kmer_autocorr_founder` +
     `kmer_phase_contrast`" guide. Patterns to document:
     - `0.05 – 0.20` low occupancy: nested subrepeat exists but
       only in some founder copies (TRC_104-class).
     - `0.20 – 0.50` moderate: real subrepeat in most founders.
     - `≥ 0.50` high: nearly-uniform nested tandem; verify it's
       not the degenerate "founder = period" case.
     - `0` and `n_intervals = 0`: no contiguous tandem run met
       the threshold — strong evidence against nested subrepeat
       at this period in this array.
   - Update "Status" feature list — 9 feature drops now.
3. `docs/onboarding_pipelines.md`:
   - Update column count `11 → 15` and the per-column table.
4. `docs/release.md` — add a `feat(rescore)` bullet to the v0.12
   section (or a new release section if we cut v0.12.1).

## 10. Performance budget

Per-row cost: O(L) — one pass to build `match`, one cumulative
sum, one threshold + run-finding pass. For L = 30 kb arrays and
~10 rescored rows per record on IPIP, that's ~300 k ops per
record, ~1 G ops total → estimated < 5 s extra on top of the
existing ~165 s rescore time on IPIP. Negligible.

Memory: peak ~2L per row (the `match` u8 array + the `rate` f64
array, deallocated when the row finishes). With default rayon
parallelism on 16 threads and 30 kb arrays, peak ~16 × 240 kB ≈
4 MB transient — irrelevant.

## 11. Order of operations / commit plan

Single commit per natural unit:

1. `feat(rescore): scan module + ScanConfig + scan column wiring`
   - Adds `src/rescore/scan.rs`, types, `Config.scan`, and the
     inline call in the parallel rescore loop.
   - 12 unit tests + smoke integration test update.
2. `feat(rescore-cli): --scan-* flags + main wire-up`
   - CLI args in `src/cli.rs`, plumb into `KmerSpatialConfig`
     wait no, into `ScanConfig` via `src/main.rs`.
3. `docs(rescore): describe scan_n_intervals + scan_occupancy_frac`
   - `docs/rescore.md`, `docs/onboarding_pipelines.md`,
     `README.md`, `docs/release.md`.

Optional 4th commit if we add the Python-parity test:

4. `test(rescore): python-parity for scan_one_period`

## 12. Hand-off criteria — "scan integration done"

- `cargo test --release --locked` passes including the new scan
  unit tests.
- `cargo fmt --all --check` and `cargo clippy ... -D warnings`
  clean.
- Running `kitehor rescore` on the synth corpus reproduces the
  Python prototype's per-record `scan_occupancy_frac` within
  ±2 % on every case (C05–C09).
- Running on the IPIP corpus reproduces the spot-check values
  for TRC_104 / TRC_115 / TRC_666 listed in §8.
- The TSV schema is 15 appended columns; existing column indices
  unchanged.
- `kitehor rescore --no-scan` emits NA for both new columns,
  same speed as v0.12.

## 13. Decisions (locked)

1. **Scope** — scan runs on every kite-reported `(record, period)`
   row that rescore already processes. No separate period
   bounds. Single tunable is the identity threshold.
2. **Subrepeat gating** — scan stays observational in this PR.
   Two new columns (`scan_n_intervals`, `scan_occupancy_frac`)
   are emitted alongside the existing flag, which is unchanged.
   Promotion to a gate is a future PR after we look at more
   IPIP cases via dotplot.
3. **Default `--scan-id-threshold`** — 0.55. Calibrated on
   TRC_666 (real subrepeat at ~ 12 % per-copy divergence).
   Tunable via CLI.
4. **Per-interval coordinates** — out of scope. The TSV is
   one row per peak; emitting per-interval rows would change
   the schema's grain. A sidecar `<prefix>.scan_intervals.tsv`
   is a clean future addition if downstream consumers want it.
