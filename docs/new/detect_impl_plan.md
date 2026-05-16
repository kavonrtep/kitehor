# Detector implementation plan — kitehor

A kitehor-crate-specific implementation plan for the v2 line-width
detector described in the companion docs in this folder:

- `detect_spec.md` — detector design (regime A/B/C theory, property
  vector, classification logic, milestones M1–M6).
- `taxonomy.md` — structural taxonomy + truth-file schema.
- `simulator_impl_plan.md` — sibling document for the simulator that
  produces the detector's input.
- `simulator_schema.json` — YAML schema for the simulator (the
  detector consumes its output, not this schema directly).

This document is the **kitehor implementation contract** for the
detector: where the code lives, which Rust modules and types
implement which spec section, what the CLI surface is, and the
milestone acceptance gates. Read the companion docs first for *what*
each construct means biologically; this doc is about *how* it lands
in this crate.

## 0. Decisions on the open questions

Locked in before implementation:

| Q | Decision |
|---|---|
| Q1 | Subcommand prefix = `detect`. Existing `kite-periodicity` (rule + ML) stays as the production path until the detector is validated on real data. |
| Q2 | Coexists with `kite-periodicity` — same pattern as `synth*` and `simulate*`. The detector is opt-in via a new subcommand. |
| Q3 | Visualisation **off by default**; enable per-run with `--viz-dir <dir>`. Granular flags (`--export-raster`, `--export-shift`, ...) layer on top for cherry-picking specific artefacts. |
| Q4 | Output format: 2D matrices → PNG (visual) + TSV (numeric); 1D signals → TSV only. PNG behind a Cargo feature flag (`viz`, default-on) so the `image` dep is optional for builds that don't need it. **The CLI surface is the same in both builds**: a no-viz build accepts `--viz-dir`/`--export-*` and runtime-errors with `PNG visualisation support was not compiled in` if PNG output is actually requested. TSV diagnostics still work without `image`. (Amendment from review §13 finding M3 — avoids two clap surfaces.) |
| Q5 | HMM-based segment classification **deferred to v2**. MVP uses threshold-based segmentation at phase-shift breakpoints + per-segment regime-A/B/C classification. Post-MVP HMM trigger: revisit only if M4 produces repeated false splits/merges on phase-shift, stratification, or mixed-family cases, **or** if M6 benchmark failures cluster around segmentation. (Amendment from review §13 finding Q5.) |

### Decisions resolving the parking-lot open questions

From the review's Recommended parking-lot answers (§12 originally
listed these as OPEN; they're now decided):

| OQ | Decision |
|---|---|
| OQ1 — k-mer dim | **k = 4 oriented 4-mers (256 dims)** for MVP, exposed in `DetectorConfig.embedding_k`. Sweep `{3, 4, 5}` only after M4 works; don't block M1–M3 on this. (See §6.3 — "oriented", not "canonical": input is assumed canonically oriented and strand-aware comparison is v2.) |
| OQ2 — min_hor_units | Default **`min_hor_units = 5`**. Below threshold, call HOR only if both base- and unit-period evidence is strong; otherwise prefer `ambiguous` over forcing `simple_TR`. Pinned in `DetectorConfig`. |
| OQ3 — inversion handling | Defer strand-aware recognition. For MVP, **T12's expected detector class is `irregular_HOR` or `ambiguous`** (whichever segmentation can support) — *not* clean `HOR`. Encoded in `tests/detect_expectations.tsv` (see §9.2). |
| OQ4 — stratification vs phase-shift threshold | First-pass segment-consensus rule: same `base_width`/`k` and consensus identity ≥ 0.90 → phase shift property; < 0.80 → `mixed`; 0.80–0.90 → `ambiguous` pending calibration. Thresholds in `DetectorConfig.stratification_*`. |
| OQ5 — viz scope | Render best width only by default. Add `--export-all-widths` as a separate flag later rather than overloading `--export-raster`. |
| OQ6 — HMM revisit | Tied to Q5 trigger above. |

### Amendments arising from review

This subsection consolidates the in-place changes triggered by
`docs/reviews/detect_impl_plan_review_2026-05-16.md`. Each item maps
to where in this plan the change applied.

- **A1 — detector oracle ≠ simulator truth** (review high-finding 1).
  Affects §9 (Test corpus) and §10 (M0/M4 acceptance). New
  `tests/detect_expectations.tsv` is the CI oracle; synth `truth.tsv`
  remains generative metadata only. Regime-C example: T07 generates
  as HOR but the detector should call `simple_TR` at HOR-unit width.
- **A2 — two-stage shift recovery** (high-finding 2). Affects §6.6.
  Local drift/wobble uses `S = 5`; breakpoint offset recovery uses a
  separate pass over a wider window (`s ∈ [-w/2, +w/2]` or circular
  cross-correlation) on rows flanking each candidate breakpoint.
- **A3 — width prioritisation before cap** (high-finding 3). Affects
  §6.1. Explicit ranking: input periods → divisors of top-N → near-
  misses → harmonics. Every original period and every valid divisor
  of the top-N period candidates is kept before any near-miss.
- **A4 — `inter_monomer_identity` in `properties.tsv`** (high-finding
  4). Affects §4 (data model already had the field) and §10 (output
  schema must declare it).
- **A5 — wrap/N/trailing-row semantics pinned** (medium-finding 2).
  Affects §6.2. Uppercase input; non-ACGT → N; trailing partial row
  dropped for feature extraction but counted in `length_bp`; minimum
  complete rows per width (`min_rows_per_width = 8` default); N
  excluded from IC denominators.
- **A6 — oriented 4-mers, not canonical** (medium-finding 1). Affects
  §6.3.
- **A7 — viz runtime error, not clap parse error** (medium-finding 3,
  also §0 Q4 above). Affects §8.4 and §13.
- **A8 — M4 acceptance pinned to specific core fixtures** (medium-
  finding 4). Affects §10. Replaces "≥ 90 % of CI fixtures" with
  exact-pass on T01, T05, T06, T07, T10, T13, T17, T18; percentage
  thresholds apply only to the v2 benchmark.
- **A9 — M5 consensus tested against synth diagnostics or inline-
  sequence fixtures** (medium-finding 5). Affects §10. HOR-unit
  consensus must be built from the HOR-unit width, not by repeating
  the base-width consensus.
- **A10 — new M0 milestone** (review §11). Affects §10. IO + config +
  schema scaffolding only; M1 then runs on top of a known-good IO
  layer.
- **A11 — `DetectorConfig` explicit** (review §11). New §6.0 lists
  every default; CLI flags map 1:1 onto config fields.
- **A12 — periods.tsv input validation** (review §11). Affects §7
  (CLI) and §6.0 (config).
- **A13 — `detect-batch` and multi-record FASTA** (review §11).
  Stem-pairing is fine only when each `.fa` has one record; otherwise
  pair by `array_id` joined against an array-id column in
  `periods.tsv`.
- **A14 — drop `ndarray`** (review §11). Affects §2. Row-major
  `Vec<u8>` + `(rows, cols)` is enough for M1–M3.
- **A15 — review-2026-05-16 surgical fixes**
  (`docs/reviews/detect_implementation_review_completed_2026-05-16.md`).
  Affects §5 (pipeline), §10 (output schema notes), §8 (viz flags),
  §7 (CLI), §6.11 (confidence):
  - mixed/ambiguous decisions no longer inherit pre-decision
    `base_width_bp`/`hor_k`/`hor_length_bp`/`column_conservation`/
    `phase_separation`/`inter_monomer_identity`, and no consensus
    or viz is emitted for them (review #1, high);
  - irregularity demotion now exempts HOR cases whose dominant
    abnormality is smooth wobble (`wobble_amplitude_bp / w ≥ 0.05`
    AND no phase shifts) — those stay `HOR` with wobble flagged
    in `reason` (review #4, high);
  - `properties.tsv::inter_monomer_identity` documented as an
    **approximation** carrying `R(1)` at base width (k-mer
    composition similarity), not pairwise sequence identity
    (review #5, medium);
  - `properties.tsv::confidence` documented as a **heuristic
    score, not a calibrated probability**; mixed/ambiguous now
    derive logit from `n_complete_copies` + `n_phase_shifts`
    rather than a class constant (review #6, medium);
  - `periods.tsv::period_score` parsing is strict: empty,
    malformed, NaN/infinite, and out-of-`[0,1]` values are hard
    errors instead of being coerced to `0.0` (review #8, medium);
  - viz flags made exact: `--viz-dir` alone emits all cheap TSVs
    (back-compat); any granular flag switches to per-flag gating,
    and `--export-edges` now also writes a per-row `edge_matrix`
    TSV (review #9, medium);
  - single-run `detect` now mirrors batch-mode DH11: periods rows
    whose `array_id` matches no FASTA record are a hard error
    unless `--allow-extra-periods` is set (review #11, low).
- **A18 — kite emit-periods review fixes**
  (`docs/reviews/kite_emit_periods_integration_review_2026-05-16.md`).
  Tightens A17:
  - **Top-3 cap** is enforced on the *source set*, not the *output
    count*. Secondaries are drawn only from `kr.peaks[0..3]` minus
    founder/tile — never from rank 4+ (review #1, high). Closes
    the ambiguity in the original "other top-3 peaks" phrasing.
  - **ML classifier conflict.** `--emit-periods` is now mutually
    exclusive with `--use-ml-classifier` at the clap level —
    silently falling back to raw kite peaks while the ML path
    held founder/tile in memory was confusing (review #2, medium).
  - **QC-skipped records** documented as a second case where
    `--allow-missing-periods` is required on the detector side
    (alongside `NoSignal`). README + CLAUDE.md updated; CLI help
    text spells out both cases (review #3, medium).
  - **Tandem source label** renamed `kite_founder` → `kite_monomer`
    (review #5, low). Detector ignores `source` so this is a
    user-facing relabel only.
  - **End-to-end test** added at `tests/detect_kite_emit.rs`:
    runs `kite-periodicity --classify --emit-periods` then
    `detect --periods` on the committed smoke fixture, asserts
    v2-schema header, periods-row sanity, and a 20-column
    properties.tsv with at least one resolved class. Also
    asserts the `--use-ml-classifier` conflict (review #4, medium).
- **A17 — kite → detector integration**
  (settled in the 2026-05-16 integration discussion). Affects §7 (CLI).
  `kitehor kite-periodicity --emit-periods <path>` writes a v2
  `periods.tsv` directly so the same FASTA can run end-to-end as
  `kite-periodicity --emit-periods` → `detect --periods`. Score
  mapping in `src/emit_periods.rs`, chosen relative to
  `DetectorConfig::strong_period_score` (0.85):
  - founder → 0.95 (`kite_founder`),
  - tile → 0.90 (`kite_tile`, only when ≠ founder),
  - other top-3 kite peaks → 0.60 (`kite_secondary`),
  - `Unresolved` verdict → top-3 hints @ 0.50 / 0.40 / 0.30
    (`kite_peak`),
  - `NoSignal` → no rows (use `detect --allow-missing-periods`
    to keep the record in the output as `ambiguous`),
  - no `--classify` → raw top-3 @ 0.60 hints (`kite_peak`).
  Tradeoff: separate emit-periods → detect glue keeps each stage
  individually debuggable (you can inspect the periods file
  before running detect). A combined `kitehor analyze` subcommand
  is a possible follow-up once the integration is stable.
- **A16 — deferred from review-2026-05-16** (high findings #2 + #3).
  Not yet addressed; documented here so the trade-off is explicit:
  - **Per-segment recompute (#2).** `segment::split()` still emits
    boundary splits inheriting whole-array class/base_width/k.
    A real implementation needs per-segment wrap, R(k), shift,
    consensus, and stratification-threshold comparison
    (`stratification_same_threshold` / `stratification_diff_threshold`
    are already in `DetectorConfig` but inert). Planned for M7.
  - **Same-width mixed detection (#3).** Today's mixed detection
    relies on distinct high-score input periods or incompatible
    candidate widths; same-width / same-`k` mixed families
    (e.g., `mx_a200-08_b200-08_n050-050`) collapse to a single
    HOR. Fix requires segment-level consensus identity comparison
    using the stratification thresholds — strictly depends on the
    per-segment recompute above. Planned for M7.

## 1. Scope, naming, language

- **Language**: Rust (consistent with the rest of kitehor; the
  upstream `detect_spec.md` §14 calls for Rust).
- **Subcommand prefix**: `detect`. New commands:
  - `kitehor detect <fasta> --periods <tsv> -o PREFIX` — one or many
    arrays.
  - `kitehor detect-batch --fasta-dir DIR --periods-dir DIR --out-dir
    DIR` — batch over many FASTA + period-TSV pairs.
- **Coexists with existing pipeline**: `kite-periodicity` (rule.rs +
  classifier.rs + classify.rs) stays untouched. The new detector is
  additive.
- **Inputs**: FASTA of pre-extracted tandem-repeat arrays + per-array
  period-candidate TSV (the schema `kitehor synth` already emits as
  `{prefix}.periods.tsv`). No assumptions about which upstream tool
  produced the candidates.

## 2. Cargo dependencies to add

The crate already has: clap, needletail, rayon, rustfft, serde,
serde_json, serde_yaml, toml, csv, ahash, anyhow, thiserror, log,
env_logger, rand, rand_chacha, rand_distr, jsonschema. Add:

```toml
image = { version = "0.25", optional = true }   # PNG export (Q4)

[features]
default = ["viz"]
viz = ["dep:image"]
```

`image` is feature-gated; `cargo build --no-default-features`
produces a slim binary without the PNG codec. The CLI surface stays
identical between the two builds: `--viz-dir` and `--export-*` are
always accepted, but a no-viz build returns a runtime error if PNG
output is actually requested (see A7 in §0). TSV diagnostics work in
both builds.

**No `ndarray` for MVP** (A14). Row-major `Vec<u8>` plus an
explicit `(rows, cols)` pair is sufficient for M1–M3; revisit only
when a concrete bottleneck demands the slicing/broadcasting
machinery. `rustfft` is already in the crate (used by the wobble
detector); the detector reuses it for `best_shift(r)` FFT
(spec §7.6).

## 3. Module layout

New module tree under `src/detect/`:

```
src/detect/
  mod.rs               public API: detect(fasta, periods) -> DetectorOutput
  io.rs                FASTA loader, periods.tsv reader, output writers
  widths.rs            candidate width expansion (period ± N + divisors)
  wrap.rs              wrap sequence to 2D + background-corrected column IC
  embed.rs             k-mer row embeddings (default k=4, 256-dim, L2-normalised)
  autocorr.rs          R(k) over row embeddings
  edges.rs             diff_x, diff_y, column-edge profile + lag autocorr
  shift.rs             per-row best_shift(r); drift / wobble / breakpoint analysis
  phase.rs             phase_separation, primitive-multiplicity correction
  irregularity.rs      block-level features → irregularity_score
  segment.rs           breakpoint segmentation + per-segment recompute
  classify.rs          array-level decision (regime A/B/C)
  confidence.rs        sigmoid logit formula (spec §9)
  consensus.rs         column-vote consensus monomer + HOR unit
  viz.rs               §8 — matrix/signal export (PNG + TSV)
  types.rs             shared types: ArrayId, Properties, WidthFeatures, Segment
```

Wiring: `pub mod detect;` in `src/lib.rs`. CLI dispatch in
`src/main.rs` mirrors the `synth*` pattern.

## 4. Data model

Central types in `src/detect/types.rs`:

```rust
pub struct Properties {
    pub array_id:             String,
    pub length_bp:            usize,
    pub class:                Class,
    pub base_width_bp:        usize,
    pub hor_k:                Option<usize>,
    pub hor_length_bp:        Option<usize>,
    pub n_complete_copies:    usize,
    pub column_conservation:  f64,
    pub phase_separation:     f64,
    pub mean_shift_bp:        f64,
    pub wobble_amplitude_bp:  f64,
    pub wobble_periodicity_bp:Option<f64>,
    pub n_phase_shifts:       usize,
    pub phase_shift_positions:Vec<usize>,
    pub phase_shift_offsets:  Vec<i64>,
    pub irregularity_score:   f64,
    pub inter_monomer_identity: Option<f64>,
    pub confidence:           f64,
    pub n_segments:           usize,
    pub reason:               String,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Class { SimpleTR, HOR, IrregularHOR, Mixed, Ambiguous }

pub struct WidthFeatures {
    pub width_bp:            usize,
    pub n_rows:              usize,
    pub column_ic:           f64,
    pub fraction_conserved:  f64,
    pub r_lag1:              f64,
    pub best_lag:            usize,
    pub best_lag_score:      f64,
    pub phase_separation:    f64,
    pub vertical_edge_rate:  f64,
    pub col_edge_autocorr_k: usize,
    pub col_edge_autocorr_s: f64,
    pub mean_shift_bp:       f64,
    pub wobble_amplitude_bp: f64,
    pub n_phase_shifts:      usize,
    pub irregularity_score:  f64,
    pub class_hint:          ClassHint,
}

pub enum ClassHint {
    SimpleTRBaseWidth,
    HORBaseWidth { k: usize },
    HORUnitWidth,
    UnsupportedWidth,
}

pub struct Segment {
    pub array_id:    String,
    pub segment_id:  usize,
    pub start_bp:    usize,
    pub end_bp:      usize,
    pub class:       Class,
    pub base_width_bp: usize,
    pub hor_k:       Option<usize>,
    // ...properties subset, per spec §12 segments.tsv schema
}

pub struct DetectorOutput {
    pub properties:       Vec<Properties>,
    pub segments:         Vec<Segment>,
    pub width_features:   Vec<(String, WidthFeatures)>,
    pub consensus:        Vec<ConsensusRecord>,
    pub diagnostics:      DetailedDiagnostics,
}
```

Internal per-width intermediate state (built and consumed inside
`run_one`; not part of the public output unless `--viz-dir` is set):

```rust
struct WidthState {
    width_bp:        usize,
    rows:            Vec<Vec<u8>>,        // wrapped 2D bytes (n_rows × w)
    column_ic:       Vec<f64>,            // per-column information content
    embeddings:      Vec<Vec<f32>>,       // k-mer L2-normed (n_rows × 256)
    r_k:             Vec<f64>,            // R(k) for k = 1..K
    diff_x:          Vec<u8>,             // edge-field flags (n_rows × w)
    diff_y:          Vec<u8>,
    shift:           Vec<i32>,            // best_shift(r), length n_rows
    breakpoints:     Vec<usize>,
}
```

## 5. Pipeline (`detect/mod.rs::run_one`)

```rust
pub fn run_one(
    array: &Array,
    periods: &[PeriodCandidate],
    cfg: &DetectorConfig,
    viz: Option<&VizSink>,
) -> Result<DetectorRecord> {
    // 1. Expand widths
    let widths = widths::expand(periods, cfg, array.len());

    // 2. Per-width feature extraction (parallel)
    let states: Vec<WidthState> = widths
        .par_iter()
        .map(|w| extract_width(*w, array, cfg))
        .collect();

    // 3. Width-level classification → array-level evidence
    let (best, secondary) = phase::pick_best_width(&states, cfg);

    // 4. Refine width via mean_shift (re-extract at refined w if needed)
    let best = shift::refine_width(best, array, cfg)?;

    // 5. Segment at phase-shift breakpoints; recompute per segment
    let segments = segment::split(&best, array, cfg)?;

    // 6. Combine into array-level properties
    let props = classify::decide(&states, &segments, cfg)?;
    let confidence = confidence::compute(&props, &best);

    // 7. Consensus monomer + HOR unit (column-vote)
    let consensus = consensus::build(&best, props.hor_k);

    // 8. Visualisation, if requested
    if let Some(viz) = viz {
        viz::export(viz, &states, &best, &segments, &consensus, array.id())?;
    }

    Ok(DetectorRecord { props, segments, states, consensus })
}
```

Per-stage notes:

- **Step 2** uses rayon over widths. Per-width work is independent.
- **Step 4** width refinement re-extracts the affected width only.
- **Step 5** segments use the same per-segment feature extraction;
  per-segment work is parallel too.
- **Step 8** is a no-op when `viz` is `None`.

## 6. Critical algorithms — restate the non-obvious ones

Detailed algorithms live in `detect_spec.md` §7. This section pins
down kitehor-specific implementation choices.

### 6.0 `DetectorConfig` — every threshold in one place (A11)

```rust
pub struct DetectorConfig {
    // Widths
    pub min_width: usize,                  // 20 bp
    pub max_width: usize,                  // 5_000 bp
    pub max_widths_per_array: usize,       // 40
    pub neighborhood_n: usize,             // ±3 around each candidate
    pub max_hor_k: usize,                  // 30

    // Wrap / IC
    pub min_rows_per_width: usize,         // 8 (below this, width is unsupported)
    pub ic_threshold_min: f64,             // 0.5 bits — minimum mean column IC
    pub ic_threshold_hor_base: f64,        // 0.4 — laxer because HOR base columns
                                           //  are multimodal across slots
    pub ic_threshold_hor_unit: f64,        // 0.7
    pub ic_threshold_simple_tr: f64,       // 0.7

    // Embeddings
    pub embedding_k: usize,                // 4 (oriented 4-mers → 256 dims)
    pub embedding_dim_hash: Option<usize>, // None unless feature-hashing kicks in

    // Phase / multiplicity
    pub phase_separation_threshold: f64,   // 0.15
    pub primitive_correction_delta: f64,   // 0.05
    pub min_hor_units: usize,              // 5 (OQ2)

    // Shift
    pub shift_local_range_bp: i32,         // ±5  (Pass A: drift / wobble)
    pub shift_breakpoint_threshold: i32,   // 3   (|Δ best_shift| triggers Pass B)
    pub shift_breakpoint_window_frac: f64, // 0.5 (Pass B window = ±frac × w bp)

    // Irregularity / segmentation
    pub block_size_rows_min: usize,        // 100  or  L/50  (whichever is larger)
    pub stratification_same_threshold: f64, // 0.90 (≥ → phase shift property)
    pub stratification_diff_threshold: f64, // 0.80 (< → `mixed`)
                                            // in between → `ambiguous` (OQ4)

    // Confidence (sigmoid logit weights, calibrated on ground_truth_v2/)
    pub confidence_weights: ConfidenceWeights,
}
```

CLI flags map 1:1 onto these fields where exposed; a `--config
detect.toml` flag loads overrides. All defaults match the upstream
spec §7–§9.

### 6.1 Candidate width expansion (`widths.rs`)

For each input period `p`: include `p ± N` (default `N = 3`), plus
divisors of `p` ≥ `min_width` (default 20 bp) and ≤ `max_width`
(default 5 kb). Divisors matter because the strongest period from
the generator may be the HOR unit length while the true base unit is
a divisor.

**Prioritisation (A3) — ranking before the `max_widths_per_array`
truncation**:

1. **Every original input period**, sorted by `period_score` desc.
2. **Every valid divisor of the top-`N` input periods** (default
   `N = 5`), where "valid" means `min_width ≤ d ≤ max_width` and
   `d != p` (the divisor must be a proper sub-period).
3. **Near-misses** (±`N` neighborhood) around the items already in
   tiers 1 and 2.
4. **Harmonics** (`2·p`, `3·p`) and any remaining low-score extras.

Cap at `max_widths_per_array` only after tier 3 fills. Tiers 1 + 2 are
**never** dropped by the cap — if they collectively exceed the cap,
raise the cap for that array and log a warning. This guarantees the
true base-width divisor reaches feature extraction even when many
distractor periods are in the input.

### 6.2 Wrap + background-corrected column IC (`wrap.rs`)

**Wrap semantics (A5)** — pin every edge case explicitly so that IC,
embeddings, edges, and consensus all see the same alphabet:

1. Uppercase the input sequence on load.
2. Encode `A` → 0, `C` → 1, `G` → 2, `T` → 3, everything else
   (including `N`, IUPAC codes, lower-case leftovers) → 4 (`N`).
3. At width `w`, full rows are `n_rows = length_bp / w` (integer
   division). The **trailing partial row is dropped** from feature
   extraction; `length_bp` in the property table still reflects the
   original FASTA length.
4. If `n_rows < min_rows_per_width` (default 8), mark the width
   `UnsupportedWidth` and skip the per-width compute path.
5. **N is excluded from IC denominators**: per-column,
   `p_b = count_b / (n_rows − n_N)` for `b ∈ {A,C,G,T}`. When
   `n_N == n_rows`, IC is 0 by definition.
6. Array-wide background `q_b` is computed over A/C/G/T only (Ns
   excluded). If any `q_b == 0` (no instance of that base across the
   array), fall back to the uniform `q_b = 0.25` to avoid `log(0)`.

Then per the spec §7.2: array-wide background `q_b` first, then
per-column `IC = Σ p_b · log2(p_b / q_b)`. Tolerates AT-rich
satellites without falsely inflating IC against a uniform baseline.

Width-level outputs: `mean_column_IC`, `fraction_conserved_columns`
(thresholded at `ic_threshold_min`, default 0.5 bits).

### 6.3 K-mer row embeddings (`embed.rs`)

**Oriented 4-mers (A6)** — not reverse-complement-canonicalised.
Default `k = 4` gives 4⁴ = 256 oriented k-mers, L2-normalised. Row
similarity is the dot product. Input is assumed canonically oriented
upstream and strand-aware comparison is v2 (`detect_spec.md` §15);
oriented k-mers keep the embedding sensitive to inversions when that
mode arrives. Tolerant of intra-row indels because k-mer composition
is position-invariant within the row.

If memory matters for very long arrays (`L > 5 Mb`), feature-hash to
64–128 dims via `ahash` (already a dep), controlled by
`DetectorConfig.embedding_dim_hash`.

### 6.4 R(k) (`autocorr.rs`)

```rust
R[k] = mean_i dot(emb[i], emb[i+k])    for k = 1..K
```

K capped at `max_hor_k` (default 30). O(n_rows × K × d).

### 6.5 Edge fields (`edges.rs`)

`diff_x[r, c]` and `diff_y[r, c]` per `detect_spec.md` §7.5.
Aggregate `horizontal_edge_rate`, `vertical_edge_rate`, plus the
per-column `column_edge_rate[c]` whose autocorrelation along `c`
gives an *independent* vote on HOR multiplicity (two-source
confidence boost).

### 6.6 Shift signal (`shift.rs`) — two-stage (A2)

A single shift search over `s ∈ [-5, +5]` is enough to recover drift
and wobble but **cannot estimate large phase-shift offsets** (T10's
85 bp shift on a 171 bp monomer is well outside the small window).
Split the logic into two passes:

**Pass A — local drift / wobble** (per spec §7.6 with `S =
shift_local_range_bp`, default 5). For each adjacent row pair:

```text
match(r, s)  = fraction of c where S[r·w + c] == S[(r+1)·w + c + s]
best_shift(r) = argmax_s match(r, s)        for s ∈ [-S, +S]
```

Use `best_shift(r)` only for:

- `mean(best_shift)` → drift → width refinement.
- `std(best_shift)` after removing breakpoint segments → wobble amplitude.
- 1D FFT of detrended `best_shift` (via `rustfft`) → wobble periodicity.
- **Breakpoint *detection*** (not offset recovery):
  `|Δ best_shift(r)| ≥ shift_breakpoint_threshold` (default 3) or a
  single-row drop in `match(r, best_shift(r))` flags a candidate
  breakpoint at row `r`.

**Pass B — phase-shift offset recovery** (the new piece). For each
candidate breakpoint at row `r`, take a row window before
(`rows r−B..r`) and after (`rows r+1..r+1+B`) where `B` is small
(e.g. 8 rows or `block_size_rows_min / 8`). Compute the inter-window
shift via **circular cross-correlation** of the column-mean profile
of each window (length `w`):

```text
cc[s] = Σ_c  pre_col_mean[c] * post_col_mean[(c + s) mod w]
phase_shift_offset = argmax_s cc[s]    for s ∈ [-w/2, +w/2]
```

This delivers the offset in bp on the same scale as `offset_bp` in
the simulator's truth. Wide-window search via FFT is O(w log w) per
breakpoint — cheap even for kb-scale monomers.

Reported in `properties.tsv`:

- `n_phase_shifts` = number of Pass-A breakpoints.
- `phase_shift_positions` = breakpoint row positions converted to bp.
- `phase_shift_offsets` = Pass-B offsets (one per breakpoint).

This is the most complex single module; landing it as two PRs (Pass
A then Pass B) is the natural split.

### 6.7 Phase separation + primitive correction (`phase.rs`)

`phase_separation = R(best_k) − median R(k')` for k' near best_k but
not a multiple. Threshold default 0.15.

Primitive correction (`detect_spec.md` §7.8): if `k = 12` and `k = 6`
both show strong R, prefer the smaller. `δ` margin default 0.05.

### 6.8 Irregularity (`irregularity.rs`)

Per-block: column IC, best lag, phase separation, mean shift. Block
size `B = max(100, L / 50)` rows. Aggregate as weighted block-to-block
variance.

### 6.9 Segmentation (`segment.rs`)

If `n_phase_shifts > 0`, split at shift positions and recompute
§7.2–7.9 per segment. If all segments share `base_width` and `k`, the
array inherits that class; phase shifts become properties. If
segments differ in `base_width` or `k`, class becomes `mixed`.

(MVP: threshold-based only. HMM deferred per §0 Q5.)

### 6.10 Classification logic (`classify.rs`)

Verbatim from `detect_spec.md` §8. The regime-A/B/C decision is the
load-bearing piece — HOR is called only when *both* the base period
and the HOR-unit period are statistically valid.

### 6.11 Confidence (`confidence.rs`)

Sigmoid logit with weights α..ζ from `detect_spec.md` §9. Calibrate
on `ground_truth_v2/` such that:

- Clean HOR cases (categories 02, 04, 05, 11) → confidence ≥ 0.9.
- Ambiguous cases (random, regime-C boundary) → ≈ 0.5.
- Negative controls (08_random) → ≤ 0.2.

Weights stored as `pub const ALPHA: f64 = ...;` in `confidence.rs`
with documented derivation comments. Make every threshold a CLI flag
with the documented default.

### 6.12 Consensus (`consensus.rs`)

Column-vote on `n_rows × w` matrix at the chosen width. Plurality
base per column; tie-breaking by alphabetic order. Output:

```
>{array_id}.monomer  length={w}
ACGT...
>{array_id}.hor_unit length={w·k}  k={k}
ACGT...
```

## 7. CLI surface

```rust
// in cli.rs
pub enum Command {
    // ... existing variants ...
    Detect(DetectArgs),
    DetectBatch(DetectBatchArgs),
}

pub struct DetectArgs {
    /// Input FASTA (one or many records).
    pub fasta: PathBuf,
    /// Period candidates TSV (matches `kitehor synth` output schema).
    #[arg(long)]
    pub periods: PathBuf,
    /// Output prefix; writes PREFIX.properties.tsv, .segments.tsv,
    /// .width_features.tsv, .consensus.fa, .diagnostics.json.
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Visualisation: emit per-array matrices to this directory.
    #[arg(long)]
    pub viz_dir: Option<PathBuf>,
    /// Override default config (TOML).
    #[arg(long)]
    pub config: Option<PathBuf>,
    /// Number of rayon worker threads (0 = auto).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
    /// Granular viz flags (any subset; implies --viz-dir).
    #[arg(long)]
    pub export_raster: bool,
    #[arg(long)]
    pub export_shift: bool,
    #[arg(long)]
    pub export_edges: bool,
    #[arg(long)]
    pub export_ic: bool,
}

pub struct DetectBatchArgs {
    #[arg(long)]
    pub fasta_dir: PathBuf,
    #[arg(long)]
    pub periods_dir: PathBuf,
    #[arg(long)]
    pub out_dir: PathBuf,
    #[arg(long)]
    pub viz_dir: Option<PathBuf>,
    #[arg(long, default_value_t = 0)]
    pub threads: usize,
}
```

`detect-batch` pairs each `<fasta_dir>/<stem>.fa` with
`<periods_dir>/<stem>.periods.tsv`; mismatched stems are a hard
error (not silently skipped). Parallelises over arrays with rayon.

**Multi-record FASTA semantics (A13)**: stem-pairing only works when
each `.fa` contains exactly one record. For multi-record inputs
(e.g., upstream `mafft` outputs), `periods.tsv` must contain an
`array_id` column that matches the FASTA record IDs and the detector
joins on that — `<stem>.periods.tsv` then contains rows for every
record in `<stem>.fa`. The loader errors out if a record in the
FASTA has no matching `array_id` rows (or vice versa).

**`periods.tsv` validation (A12)**:

| Column        | Required? | Notes |
|---|---|---|
| `array_id`    | Yes (multi-record) / No (single-record per file, defaults to file stem) | Must match a FASTA record ID for that array |
| `period_bp`   | Yes | Integer in `[min_width, max_width]` |
| `period_score`| Yes | f64; only used for ranking (§6.1 A3) |
| `source`      | Optional | Documentation-only |

Duplicate `(array_id, period_bp)` rows → warning, kept as a single
entry with `period_score` = max of duplicates. Missing `array_id` in
a multi-record context → hard error. Out-of-range `period_bp` →
warning, dropped.

## 8. Visualisation & matrix export (`viz.rs`)

**Goal**: let the user `eog the.png` an array and see what the
detector saw. Diagnostics for paper figures, manual labelling, and
calibration-set inspection.

### 8.1 What gets emitted

Per array, into `<viz_dir>/<array_id>/`:

| File | Stage | Format | Default |
|---|---|---|---|
| `raster_w{w}.png`              | wrapped sequence at width *w*               | PNG (4-colour A/C/G/T) | best width only |
| `raster_w{w}.tsv`              | same matrix, numeric (A=0, C=1, G=2, T=3, N=4) | TSV                  | with `--export-raster` |
| `column_ic_w{w}.png`           | per-column IC, heatmap (1-pixel-tall strip) | PNG (greyscale)        | best width only |
| `column_ic_w{w}.tsv`           | per-column IC, numeric                      | TSV                    | always (cheap) |
| `edges_w{w}.png`               | `diff_y` field at best width                | PNG (binary)           | with `--export-edges` |
| `column_edge_rate_w{w}.tsv`    | per-column vertical edge rate               | TSV                    | always |
| `rk_w{w}.tsv`                  | R(k) curve, k = 1..K                        | TSV                    | always |
| `shift_w{w}.tsv`               | best_shift(r), one row per array row        | TSV                    | always |
| `segments.tsv`                 | per-segment props (if `n_segments > 1`)     | TSV                    | always |

PNG only when `--viz-dir` is set or any `--export-*` flag is set. TSV
artefacts that are cheap to compute (column_ic, column_edge_rate,
rk, shift) are always emitted when `--viz-dir` is set, regardless of
PNG flags — they're useful even without visual inspection.

For batch mode, `<viz_dir>/<array_id>/` is keyed off the per-record
`array_id`, mirroring the simulator's per-array prefix scheme.

### 8.2 PNG colour scheme

Wrapped raster at width *w*:

```
A → #2ca02c (green)
C → #1f77b4 (blue)
G → #d62728 (red)
T → #9467bd (purple)
N / other → #999999 (grey)
```

These are the matplotlib defaults; one of them per base. Rows are
laid top-to-bottom in array order. Width = `w` pixels; height =
`n_rows` pixels. For large arrays (`n_rows > 4096`), downsample
vertically by averaging blocks of `n_rows / 4096` rows; record the
downsample ratio in the PNG metadata (PNG tEXt chunk).

Column-IC strip: greyscale (`IC → 255·min(IC/2, 1)`), 8 px tall.

Edge field: binary (white = edge, black = no edge), same dimensions
as the raster.

### 8.3 TSV format

All TSV matrix outputs share a header row:

```
# array_id={…}  width_bp={w}  n_rows={n}  schema_version=1
row\tcol_1\tcol_2\t...\tcol_w
0\t0\t1\t...
1\t0\t1\t...
```

Numeric encoding: `A=0, C=1, G=2, T=3, N=4`. Same scheme used by
existing kitehor diagnostic outputs for consistency.

1D signals (`column_ic`, `column_edge_rate`, `shift`, `rk`) are
two-column TSVs: `index\tvalue`.

### 8.4 Cargo feature flag (A7)

PNG output is gated behind the `viz` Cargo feature (default-on):

```toml
[features]
default = ["viz"]
viz = ["dep:image"]
```

CI builds with `--no-default-features` to verify the lean path; the
release pipeline ships with `viz` on. Modules that touch PNG live
inside `#[cfg(feature = "viz")]` blocks; **TSV output paths are
always compiled**, including `--viz-dir` discovery / directory
creation / TSV writers.

**CLI surface is identical between the two builds**: `--viz-dir`,
`--export-raster`, `--export-edges`, `--export-ic` are always
accepted by clap. In a no-viz build, the runtime emits a single
clear error if any PNG-producing path fires (raster, edges, IC
heatmap):

```text
error: PNG visualisation support was not compiled in
       (`viz` Cargo feature disabled). Re-build with `--features viz`
       or drop --viz-dir / --export-raster / --export-edges / --export-ic.
       TSV diagnostics still work without re-building.
```

This avoids maintaining two clap surfaces and gives a clear path
forward for users on the slim build.

### 8.5 Tests for viz

- **Round-trip**: `raster_w{w}.tsv` → reload → bytes match the
  original wrapped sequence at width *w*.
- **Snapshot**: a small fixture's `raster_w100.png` is binary-equal
  across runs (PNG encoder is deterministic given the same input).
- **Slim build**: `cargo build --no-default-features` succeeds and
  the binary lacks `--viz-dir` (clap rejects the flag at parse time).

## 9. Test corpus + evaluation

### 9.1 Two ground-truth sources

1. **`tests/synth_configs/`** — 22-fixture CI corpus.
2. **`ground_truth_v2/`** — 1,600-case benchmark corpus.

The simulator's `truth.tsv` is **generative metadata** — it records
what the simulator built, not what the detector is expected to call.
Several cases differ:

- **Regime C** (e.g. T07: HOR with `div = 0.7`). Generated as HOR;
  the detector should call `simple_TR` at the HOR-unit width.
- **Inversion** (T12). Generated as HOR-with-INVERSION; pending
  strand-aware recognition (v2, OQ3), the detector's expected output
  is `irregular_HOR` or `ambiguous`, not clean `HOR`.
- **Random control** (T17). Generated as `random`; detector should
  call `ambiguous`.

### 9.2 Detector oracle: `tests/detect_expectations.tsv` (A1)

The CI oracle is a separate manifest joined to the simulator
fixtures by `array_id`. Schema:

| column                    | type            | notes |
|---|---|---|
| `array_id`                | str             | matches the synth fixture stem |
| `expected_class`          | str             | one of {`simple_TR`, `HOR`, `irregular_HOR`, `mixed`, `ambiguous`} |
| `expected_base_width_bp` | int / `NA`      | `NA` if `expected_class ∈ {mixed, ambiguous}` |
| `expected_hor_k`          | int / `NA`      | `NA` if class != HOR / irregular_HOR |
| `base_width_tol_bp`       | int             | absolute tolerance (e.g. 3 bp) |
| `expected_reason_contains`| str (or empty)  | substring assertion against the `reason` field — e.g. `"regime C"` for T07 |
| `notes`                   | str (or empty)  | free-form, ignored by CI |

The eval harness (§9.3) joins `properties.tsv` against this oracle
and reports per-fixture pass/fail. CI tests under `tests/detect_*.rs`
import the same file via a small loader.

The simulator's `truth.tsv` is **still useful** — joined alongside
the oracle, it gives the eval harness access to generative
parameters (mutation rate, divergence, etc.) for stratified
reporting.

### 9.3 Eval harness

Small Python tool under `tools/detect_eval/`. Joins:

- `manifest.tsv` (v2 corpus: generative parameters)
- `tests/detect_expectations.tsv` (CI: expected detector output)
- detector's `properties.tsv`

Reports:

- Class accuracy stratified by **expected** class (HOR / simple_TR /
  irregular_HOR / mixed / ambiguous).
- `base_width_bp` error: median absolute error, exact-match rate,
  within-tolerance rate.
- `hor_k` recovery rate.
- Per-category breakdown (the 9 v2 corpus categories).
- Per-regime breakdown on the divergence sweep (T08 in CI + the
  category 02 sub-sweep in v2).

## 10. Milestones & acceptance gates

Mirrors `detect_spec.md` §13, with kitehor-side acceptance. Adds a
new **M0** for IO/config scaffolding (A10) and tightens M4 and M5.

| # | Milestone | CI fixtures | v2 benchmark | Cargo deliverable |
|---|---|---|---|---|
| **M0** | IO + config + schema scaffolding (A10) | every fixture loads | every benchmark config loads | `kitehor detect` reads FASTA + `periods.tsv`, validates `DetectorConfig`, writes empty-but-schema-correct `properties.tsv` / `segments.tsv` / `width_features.tsv`. No detection logic yet — but every column in §10.1 is present with `NA` values. |
| M1 | Widths + column IC | T01, T05 | category 01, 02 base widths | `width_features.tsv` populated with `column_ic`, `fraction_conserved`; no class yet |
| M2 | Row embeddings + R(k) | T05 (k=12), T06 (regime A) | category 02 sweep | `r_lag1`, `best_lag`, `best_lag_score` in width features; oriented 4-mers per A6 |
| M3 | Edge field + shift signal (Pass A) | T03 (wobble), T10 (phase shift breakpoint detection only) | category 03, 04 | `mean_shift_bp`, `wobble_amplitude_bp` accurate within tolerance; `n_phase_shifts` recovers count on T10 (offsets land in M3.5) |
| M3.5 | Shift Pass B — phase-shift offset recovery (A2) | T10 | category 04 | `phase_shift_offsets` matches simulator truth within ±5 bp on the entire shift sweep |
| **M4** | Array classification + segmentation | **exact pass on T01, T05, T06, T07, T10, T13, T17, T18** (A8) | categories 01–07, 09 | every listed core fixture produces the `expected_class`/`expected_base_width_bp`/`expected_hor_k` from `detect_expectations.tsv`; reason field substring-matches when required (regime C → `"regime C"`); divergence sweep T08 a→f hits the expected regime transitions |
| M5 | Consensus + diagnostics (A9) | T05 (HOR clean), plus a deterministic inline-sequence fixture (new — added under `tests/synth_configs/T19_inline_consensus.yaml`) | category 02 | `consensus.fa` round-trips: **(a)** monomer consensus at `base_width` matches the simulator's diagnostics slot consensuses within 5 % per-base disagreement; **(b)** HOR-unit consensus is built from the HOR-unit-width column votes, **not** by repeating the base-width consensus |
| M6 | Property-level eval harness | — | full 1600-case | `tools/detect_eval/` reports a single summary tsv; weights in `confidence.rs` calibrated; **CI core fixtures must pass exactly**; v2 benchmark class accuracy ≥ 92 % overall and ≥ 88 % within each of the 9 categories |

Each milestone closes a PR. CI gate per milestone: `cargo test
--release detect::` plus an integration test that runs M*N* against
its assigned CI fixtures and asserts the documented property.

### 10.1 `properties.tsv` output schema (frozen at M0, per A4)

```
array_id  length_bp  class  base_width_bp  hor_k  hor_length_bp
n_complete_copies  column_conservation  phase_separation  mean_shift_bp
wobble_amplitude_bp  wobble_periodicity_bp  n_phase_shifts
phase_shift_positions  phase_shift_offsets  irregularity_score
inter_monomer_identity  confidence  n_segments  reason
```

`inter_monomer_identity` is `NA` unless `class = HOR` (or
`irregular_HOR` with a meaningful estimate); for `simple_TR`,
`mixed`, `ambiguous` it stays `NA`. Derived from `R(1)` at
`base_width` after accounting for the embedding's similarity floor
(per `detect_spec.md` §4).

> **Note (A15).** Today's implementation publishes `R(1)` directly
> as `inter_monomer_identity` — i.e., k-mer-composition row
> similarity, NOT mean pairwise sequence identity between inferred
> slot consensuses. The two values are correlated but not the same;
> downstream consumers should treat this column as a regime
> indicator rather than a calibrated biological identity. Schema
> is frozen so the column stays; a future major version will rename
> or recompute it.

> **Note (A15).** `confidence` is a heuristic per-class signal-
> quality score, NOT a calibrated probability. Wrong calls with
> high `confidence` are not bugs — they reflect the score
> reporting "the evidence for this class was strong" when the
> class itself was misidentified upstream. For `mixed`/`ambiguous`,
> `confidence` now varies with `n_complete_copies` and
> `n_phase_shifts` instead of being a class constant.

Schema is **frozen at M0**: every later milestone fills more columns
with real values, but no column is added/removed after M0 lands
without bumping a `schema_version` field.

## 11. Testing strategy

Test pyramid mirrors the simulator's.

**Unit tests** (`#[cfg(test)] mod tests` per module):

- `widths::expand`: divisor + neighbourhood logic on synthetic
  period lists; cap behaviour at `max_widths_per_array`.
- `wrap::column_ic`: known matrix → expected IC values.
- `embed::row_similarity`: identical rows → similarity 1.0;
  orthogonal random rows → similarity ≈ 0.
- `autocorr::r_of_k`: clean HOR matrix → R(k) peak at multiplicity.
- `shift::best_shift`: synthetic shifted-row pair → exact `s` recovered.
- `shift::breakpoints`: step function with one jump → one breakpoint
  at the right row.
- `phase::primitive_correction`: deliberate k=6 and k=12 peaks →
  prefer k=6.
- `classify::regime_a_b_c`: clean A, B, C inputs each get the right
  class.

**Integration tests** (`tests/detect_*.rs`):

- `detect_t01_simple_tr`: load `tests/synth_configs/T01_simple_tr.yaml`,
  run synth, run detect, assert `class=simple_TR`, `base_width=170`.
- `detect_t05_hor_clean`: same flow → `class=HOR`, `k=12`,
  `base_width=171`.
- `detect_t10_phase_shift`: → `n_phase_shifts=1`, position within
  ±100 bp.
- Plus T03, T04, T07, T13, T16 spot-checks.

**Determinism tests**:

- Same FASTA + same periods + same config → byte-identical
  `properties.tsv`.

**Property-level eval** (run manually, not in CI; too slow):

- Drive `tools/detect_eval/` over `ground_truth_v2/out/` and assert
  thresholds from §10 M6.

**Visualisation tests**:

- Per §8.5: TSV ↔ raster round-trip; PNG snapshot equality; slim
  build excludes viz code.

## 12. Open questions (parking lot)

OQ1–OQ6 are now resolved (see §0 decision table). What remains:

| # | Question |
|---|---|
| OQ7 | When to graduate the eval harness from `tools/detect_eval/` (Python) to a Rust subcommand (`kitehor detect-eval`)? Cheap follow-up after M6 if it gets used routinely. |
| OQ8 | Confidence-weight learning vs hand-pinned. M6 calibrates them against `ground_truth_v2/`; if the calibrated weights converge stably across reruns, freeze them. If they drift with corpus changes, consider a logistic-regression fit baked into the binary. |
| OQ9 | Should `--viz-dir` use one-array-per-subdir (`<viz>/<array_id>/...`, current spec) or flat-with-prefixed-names (`<viz>/<array_id>.raster.png`, ...)? Subdir is cleaner for batches >100 arrays; flat is friendlier for one-off inspection. Settle when first batch >100 arrays runs. |
| OQ10 | If M4 produces repeated false splits/merges on phase-shift / stratification cases — that's the HMM trigger from §0 Q5. Capture failing CI fixtures here as they appear. |

## 13. Hand-off criteria — "detector done"

1. All CI integration tests (`tests/detect_*.rs`) pass under
   `cargo test --release`. The core fixtures T01, T05, T06, T07,
   T10, T13, T17, T18 pass **exactly** against
   `tests/detect_expectations.tsv` (A8).
2. `cargo test --release detect::` lib tests are green.
3. Determinism: re-running `kitehor detect` on a fixed FASTA + periods
   yields byte-identical `properties.tsv`.
4. `tools/detect_eval/` summary on `ground_truth_v2/out/` meets the
   M6 thresholds (class acc ≥ 92 % overall, ≥ 88 % per category;
   base_width median abs err ≤ half the smallest monomer length
   tested).
5. `--viz-dir` produces the §8 file set; PNGs open in a standard
   image viewer; TSVs parse in pandas/awk.
6. `cargo build --release --no-default-features` succeeds. The
   resulting binary **accepts** `--viz-dir` / `--export-*` flags at
   the CLI surface (one clap definition shared with the default
   build, A7); requesting any PNG-producing artefact emits the
   runtime error documented in §8.4. TSV-only viz output is fully
   functional in the no-default-features build.
7. README and CLAUDE.md updated with the new `kitehor detect*`
   commands; `docs/new/detect_impl_plan.md` (this file) updated with
   any amendments arising from review.

## 14. Out of scope for MVP

Explicit non-goals (mirror `detect_spec.md` §15 plus kitehor-side):

- Hilbert curve / 2D FFT / Hough/Radon representations.
- CNN or any learned classifier.
- Dense L² dotplot.
- Whole-genome tandem-array discovery (input is pre-annotated).
- Reverse-complement orientation search (input is canonical).
- Variant-HOR reconstruction (which monomer occupies which slot).
- Species-specific centromere annotation.
- Live / interactive visualisation (the PNG export is enough for MVP).
- Real-time streaming detection.
- HMM-based segment classification (`§0 Q5`, deferred to v2).

These can be added in later iterations if the line-width features
prove predictive on real data.
