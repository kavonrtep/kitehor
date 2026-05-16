# Simulator implementation plan — kitehor

A kitehor-crate-specific implementation plan for the v2 YAML-driven
tandem-repeat simulator described in the companion docs in this folder:

- `taxonomy.md` — structural categories, grammar (§2), truth schema (§5.2).
- `simulator_plan.md` — language-agnostic implementation plan.
- `simulator_schema.json` — formal config schema.
- `detect_spec.md` — the detector that will consume this simulator's output.

This document is the **kitehor implementation contract**: where the code
lives, which Rust modules and types implement which spec section, what
the CLI surface is, and the milestone acceptance gates. Read the
companion docs first for *what* each construct means biologically; this
doc is about *how* it lands in this crate.

## 0. Decisions on the open questions

The plan as a whole is approved directionally. Open questions Q1–Q8
from the first draft are resolved here. Future readers: this section
is canonical; the original §11 list now points back to here.

| Q | Decision |
|---|---|
| Q1 | Subcommand prefix = `synth`. Existing `simulate`/`simulate-grid` unchanged. |
| Q2 | One YAML per array; multi-array runs via `synth-batch` over a directory. |
| Q3 | Nested HOR (taxonomy A3) **deferred to v2**. T09 stays in `tests/synth/configs/` as `T09_nested_hor.deferred.yaml` (placeholder doc, not generated). Do **not** fake nested HOR as a flat HOR with `k = k_inner · k_outer`. |
| Q4 | Implement `INVERSION` at the sequence level (reverse-complement of a range) now. Detector-side recognition is a separate v2 effort. |
| Q5 | Reject `\|offset_bp\| > preceding_block.monomer_length / 2` at validation time. "Preceding block" = the immediately previous `HOR`/`SIMPLE_TR` in `structure`. A `SHIFT` that follows a non-repeat block (another `SHIFT` or `INSERTION`) is an error. |
| Q6 | Keep `simulate` and `simulate-grid` alive until both (a) detector M6 lands and (b) the rule.rs eval flow migrates off `params.tsv`. |
| Q7 | `docs/new/simulator_schema.json` is the canonical schema. `src/synth/simulator.schema.json` is an embedded build asset. `tests/schema_drift.rs` asserts byte-equality. |
| Q8 | `source: file` stays in the JSON Schema, but the validator **rejects it in MVP** with `"source: file is not implemented in MVP — use 'random' or 'sequence'"`. No silent acceptance of configs that would fail at generation time. |

### Amendments arising from review

These resolve contract gaps identified in the first round of review.
Each amendment lists where it applies in the rest of the plan.

- **A1 — event block targeting** (affects §4, §6.4, schema). Every
  post-generation event (`HYBRID`, `INVERSION`, `DUPLICATION`,
  `DELETION`) carries a **required** `block: <usize>` field
  (zero-indexed into `structure`). The targeted block must be `HOR`
  or `SIMPLE_TR`. `at_copy` / `start_copy` are 1-indexed copies
  **within that block**. Without this field, multi-block configs
  (T13 coexisting periods, T14 mixed families) have ambiguous event
  coordinates. **Schema change required, lands in M1 alongside the
  validator.**
- **A2 — diagnostics is CLI-only** (affects §5, §7). `global.diagnostics`
  is *not* in the schema and will not be added. `run_one()` takes a
  `diagnostics: bool` parameter set from the `--diagnostics` CLI flag.
- **A3 — output precedence** (affects §7, schema). `-o/--out` is
  **required** for `synth`. `global.output` in YAML is silently
  ignored in MVP (the schema permits it for backwards compatibility
  with the upstream `simulator_plan.md`; the validator emits a
  non-fatal warning if present). All output paths come from the CLI.
- **A4 — invariant: one FASTA record per YAML** (forever). Multi-record
  output is not on the roadmap. Matrix runs use a directory of YAMLs
  + `synth-batch`. `truth.tsv` is one row per file. `array_id` derives
  from the file stem unless `global.array_id` is set.
- **A5 — inversion YAML shape** (affects §6.6, §8). Taxonomy's
  `H([M_1..M_12],100) + INV(H([M_1..M_12],10)) + H([M_1..M_12],100)`
  is realised as **one** `HOR` block with `n_copies = 210` plus an
  `INVERSION` post-generation event targeting `block: 0, start_copy:
  101, length_copies: 10`. The literal three-block decomposition is
  **not** used.
- **A6 — test-count language** (affects §8, §9). The corpus is **18
  conceptual tests** but the YAML-fixture count is larger because T08
  is a six-point sweep. Expected fixtures in `tests/synth/configs/`:
  22 active + 1 deferred placeholder.

## 1. Scope, naming, language

- **Language**: Rust. The kitehor crate is Rust; the existing
  `simulate` and `simulate-grid` subcommands are Rust; keeping one
  language across the project beats the marginal iteration-speed gain
  of a Python sidecar. Upstream plan §10 allows either.
- **Coexists with existing simulator**: do **not** modify or remove
  `src/simulate.rs` or `src/simulate_grid.rs`. The current rule.rs
  evaluation flow regenerates `ground_truth/sequences.fasta` from
  `ground_truth/params.tsv` via `simulate-grid` and we have not yet
  migrated that flow to YAML configs. Both simulators live side by
  side until the detector work is far enough along to drop the old
  path.
- **Subcommand prefix**: `synth` — short, distinct from `simulate`,
  no name collision. New commands:
  - `kitehor synth <config.yaml> -o PREFIX` — one array.
  - `kitehor synth-batch --config-dir DIR -o DIR` — every `*.yaml`
    in `DIR`, parallel over configs.
  - `kitehor synth-validate <config.yaml>` — schema-validate only.
  - `kitehor synth-schema --print` — emit the JSON Schema to stdout.

## 2. Cargo dependencies to add

The crate already has clap, needletail, rayon, rustfft, serde,
serde_json, toml, csv, ahash, anyhow, thiserror, log, env_logger,
proptest, tempfile. Add:

```toml
serde_yaml = "0.9"   # YAML config parsing
jsonschema = "0.18"  # validate config against simulator.schema.json
rand       = "0.8"   # PRNGs
rand_distr = "0.4"   # Normal noise for random_walk wobble
```

No new container-side tools required (cargo + rustc are pre-installed
in `bioinfo-agent.sif`). Use `htool <name>` only if a future need
arises (e.g. plotting).

## 3. Module layout

New module tree under `src/synth/`:

```
src/
  cli.rs              # add Synth, SynthBatch, SynthValidate, SynthSchema variants
  synth/
    mod.rs            # public API: run_one(cfg, out_prefix, seed) -> Result<()>
    config.rs         # serde YAML loader + jsonschema validation
    rng.rs            # FNV-1a-derived sub-stream seeds; thin wrapper over rand
    templates.rs      # HOR_slots & monomer instantiation; cache by template ID
    blocks.rs         # structure expansion: HOR, SIMPLE_TR, SHIFT, INSERTION
    wobble.rs         # residual-accumulator integer edits; sinusoidal + random_walk
    events.rs         # HYBRID, INVERSION, DUPLICATION, DELETION
    coords.rs         # CoordMap: logical (block, copy, slot) -> realised bp
    noise.rs          # final mutation + indel pass
    grammar.rs        # serialise §2 structural_expression string
    truth.rs          # truth.tsv writer; events_json builder
    periods.rs        # period candidate generator (true + distractors)
    fasta.rs          # FASTA record writer
    diagnostics.rs    # optional JSON diagnostics
    simulator.schema.json   # copy of docs/new/simulator_schema.json, embedded via include_str!
```

Embed the schema at build time:

```rust
// synth/config.rs
const SCHEMA: &str = include_str!("simulator.schema.json");
```

CI must check that `src/synth/simulator.schema.json` is byte-identical
to `docs/new/simulator_schema.json` (a small `tests/schema_drift.rs`
test does this).

## 4. Data model

Central serde types in `synth/config.rs`:

```rust
#[derive(Debug, Deserialize)]
pub struct Config {
    pub schema_version: u32,
    #[serde(default)] pub seed: u64,
    #[serde(default)] pub global: Global,
    #[serde(default)] pub templates: HashMap<String, Template>,
    pub structure: Vec<Block>,
    #[serde(default)] pub modifiers: Vec<Modifier>,
    #[serde(default)] pub post_generation: Vec<Event>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
pub enum Template {
    #[serde(rename = "HOR_slots")]
    HorSlots {
        monomer_length_bp: usize,
        k: usize,
        #[serde(default = "default_random")] source: Source,
        #[serde(default)] sequence: Option<String>,
        #[serde(default)] file: Option<PathBuf>,
        #[serde(default = "default_gc")] gc_content: f64,
        #[serde(default)] inter_slot_divergence: f64,
    },
    #[serde(rename = "monomer")]
    Monomer { /* analogous fields, k = 1 implicit */ },
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
pub enum Block {
    HOR        { template: String, n_copies: usize },
    SIMPLE_TR  { template: String, n_copies: usize },
    SHIFT      { offset_bp: i64 },
    INSERTION  { length_bp: usize, kind: InsertionKind },
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type")]
pub enum Event {
    HYBRID {
        block: usize,                    // 0-indexed into `structure`; target must be HOR or SIMPLE_TR
        at_copy: usize,                  // 1-indexed within the targeted block
        slot: usize,
        source_slots: [usize; 2],
        #[serde(default = "default_split")] split_fraction: f64,
    },
    INVERSION   { block: usize, start_copy: usize, length_copies: usize },
    DUPLICATION { block: usize, start_copy: usize, length_copies: usize },
    DELETION    { block: usize, start_copy: usize, length_copies: usize },
}
```

Internal runtime types:

```rust
pub struct CoordEntry {
    pub block_idx: usize,
    pub copy_idx: usize,        // 1-indexed, matches YAML
    pub slot_idx: usize,        // 1-indexed
    pub realised_start_bp: usize,
    pub realised_len_bp: usize,
}

pub struct CoordMap {
    entries: Vec<CoordEntry>,   // sorted by realised_start_bp
}

pub struct InstantiatedTemplate {
    pub slots: Vec<Vec<u8>>,    // realised slot consensus bytes (ACGT only)
    pub realised_inter_slot_divergence: f64,
}

pub struct SimState {
    pub sequence: Vec<u8>,
    pub coord_map: CoordMap,
    pub templates: HashMap<String, InstantiatedTemplate>,
}

pub struct TruthRow {
    pub array_id: String,
    pub length_bp: usize,
    pub truth_class: TruthClass,
    pub base_width_bp: usize,
    pub hor_k: Option<usize>,
    pub hor_length_bp: Option<usize>,
    pub n_complete_copies: usize,
    pub wobble_amplitude_bp: f64,
    pub wobble_periodicity_bp: Option<f64>,
    pub n_phase_shifts: usize,
    pub phase_shift_positions: Vec<usize>,
    pub phase_shift_offsets: Vec<i64>,
    pub n_segments: usize,
    pub reason: String,
    pub structural_expression: String,
    pub schema_version: u32,
    pub events_json: String,
}
```

## 5. Pipeline (synth/mod.rs)

```rust
pub fn run_one(
    cfg_path: &Path,
    prefix: &Path,
    seed_override: Option<u64>,
    diagnostics: bool,                // from --diagnostics CLI flag (A2)
) -> Result<()> {
    // 1. Load + schema-validate
    let cfg = config::load_and_validate(cfg_path)?;
    let top_seed = seed_override.unwrap_or(cfg.seed);

    // 2. Sub-stream seeds (FNV-1a of "top_seed:name")
    let mut streams = rng::derive_streams(top_seed);

    // 3. Instantiate templates (cached by name)
    let templates_inst = templates::instantiate(&cfg.templates, &mut streams.templates)?;

    // 4. Expand structure: build sequence + coord_map
    let mut state = blocks::expand(&cfg.structure, &templates_inst, &mut streams.structure)?;

    // 5. Modifiers: wobble (integer edits via residual accumulator)
    wobble::apply(&mut state, &cfg.modifiers, &mut streams.wobble)?;

    // 6. Post-generation events (each event names its target block — A1)
    events::apply(&mut state, &cfg.post_generation, &mut streams.events)?;

    // 7. Final noise pass (mutation + indel)
    let noise_log = noise::apply(&mut state, &cfg.global, &mut streams.noise)?;

    // 8. Emit outputs (diagnostics controlled by CLI flag, not config — A2)
    fasta::write(prefix, &cfg.global, &state.sequence)?;
    truth::write(prefix, &cfg, &state, &noise_log)?;
    periods::write(prefix, &cfg, &state)?;
    if diagnostics {
        diagnostics::write(prefix, &cfg, &state, &noise_log)?;
    }
    Ok(())
}
```

Order is normative (upstream §4.1): mutations & indels go last so the
structural truth (which template a base came from, which HOR copy,
which slot) is established on a clean sequence and only the realised
bp positions need adjustment afterwards.

## 6. Critical algorithms

### 6.1 RNG sub-streams

Use FNV-1a `(top_seed, name)` derivation, matching the convention
already established in the parent project (per parent CLAUDE.md:
*"Per-case seed is derived deterministically from master_seed:case_id
via FNV-1a hash"*). Stream names: `templates`, `structure`, `wobble`,
`events`, `noise`. This isolates the effect of each stage — changing
the noise rate does not perturb template generation, which is
essential for reproducible A/B detector calibration.

### 6.2 HOR slot derivation

For an `HOR_slots` template with `k` slots and target divergence `d`:

1. Draw slot 1 randomly with requested GC.
2. Slots 2..k: mutate slot 1 at per-base rate `d/2` (independent draws
   per slot).
3. Realised pairwise divergence is approximately `d` (small
   under-shoot at short monomer lengths). Record the *realised* mean
   pairwise divergence in diagnostics; do not back-correct (upstream
   §11.1 recommendation).

Cache instantiated templates by name. Two `HOR` blocks referencing the
same template name share slot consensuses **byte-for-byte**. This is
the load-bearing distinction between *one phase-shifted HOR* and *two
unrelated HORs that happen to have the same k*.

### 6.3 Wobble realisation (residual accumulator)

Both `sinusoidal` and `random_walk` models produce a smooth `δ(r)`
curve at row resolution:

- `sinusoidal`: `δ(r) = amplitude_bp · sin(2π · r / period_rows)`.
- `random_walk`: cumulative sum of `N(0, σ²)` noise, σ chosen so the
  smoothed series has std ≈ `amplitude_bp`. Smoothing window =
  `period_rows / 4` if `period_rows > 0`, else default 50 rows.

Convert to integer base-level edits:

```
residual = 0
for r in 0..rows:
    residual += δ(r) - δ(r-1)
    while residual >= 1:
        insert 1 base at row-r boundary; residual -= 1
    while residual <= -1:
        delete 1 base at row-r boundary; residual += 1
```

Inserted bases are sampled uniformly from the last 50 bp of
`output_sequence` (fallback: uniform random if prefix is too short).
Update `coord_map` after each edit. Record realised wobble amplitude
(std of realised `δ` series) and recovered periodicity in the truth
file.

### 6.4 Logical → realised coordinate translation

Events specify positions in logical coordinates (`at_copy: 27` means
the 27th HOR copy of the surrounding block). Translation:

1. Block expansion records `(start_bp_realised, len_bp, block_idx,
   copy_idx, slot_idx)` per slot emitted — this is the initial
   `coord_map`.
2. Wobble may insert/delete bases inside a slot; coord_map updates
   incrementally during the wobble pass.
3. Each post-generation event applies a downstream bp offset (`+len`
   for INS/DUP, `-len` for DEL, 0 for INV); subsequent events see
   updated coordinates.
4. Final noise pass produces a `noise_indel_log: Vec<(pos, delta)>`
   so truth-file positions reported in the YAML's logical terms can
   be back-computed against the final FASTA.

MVP representation: `Vec<CoordEntry>` sorted by `realised_start_bp` +
a running "cumulative indel net" delta function. Naive complexity
O(L · events). For typical inputs (a few Mb, <100 events) this is well
under 1 s; switch to interval-tree or rope only if profiling demands
(upstream §11.3).

### 6.5 Structural expression emission

Walk the original `structure` list and serialise to a §2 grammar
string. Examples:

```
single HOR (n=100, k=12, div=0.15):
  H([M_1..M_12],100,div=0.15)

phase-shifted concatenation:
  H([M_1..M_12],100,div=0.15)|shift(85)|H([M_1..M_12],100,div=0.15)

with retro-like insertion:
  H([M_1..M_12],50)+INS(5000,retro_like)+H([M_1..M_12],50)
```

Modifiers (`.mut(p).indel(q).wobble(A,T=t)`) appended to the whole
expression. Post-generation events are *not* in
`structural_expression` — they live in `events_json`. Implement as
`synth/grammar.rs::to_expression(&Config) -> String`, a pure function
over the config.

### 6.6 Inversion

`INVERSION { block, start_copy, length_copies }` is realised at the
sequence level: identify the bp range via coord_map (block index +
copy range), reverse-complement that slice in place (reverse the
order **and** complement each base), then shift `coord_map` entries
inside the range to reflect that slot order is reversed. The
detector-side *recognition* of inversions is deferred (taxonomy
§8.2), but the simulator must emit them correctly so when v2
detection lands, the test set is already there.

**YAML shape (A5)**: Taxonomy's
`H([M_1..M_12],100) + INV(H([M_1..M_12],10)) + H([M_1..M_12],100)`
is realised as a **single** `HOR` block with `n_copies = 210`
(100 forward + 10 to be inverted + 100 forward), plus one
`INVERSION` post-generation event targeting
`block: 0, start_copy: 101, length_copies: 10`. The three-block
decomposition from the grammar is *not* used — there is one
underlying HOR with a locally inverted segment, and the simulator
expresses that with one block + one event.

### 6.7 Period candidate generation (`periods.tsv`)

For each array, emit at least:
- True `base_width` with score 0.94 (`source = true_base`).
- True HOR-unit length if `k > 1`, score 0.88 (`true_hor_unit`).
- A near-miss at `base_width ± rand{2..4}`, score 0.71 (`near_miss`).
- A harmonic at `2 × base_width` or `3 × base_width`, score 0.65
  (`harmonic`).
- Optionally a random false positive between 100 and 5000 bp not
  matching any real period, score 0.42 (`false_positive`).

Score values are documentation-only; detector should not depend on
ranking. The point is to feed the detector a realistic mix.

## 7. CLI surface

Add to `cli.rs` (sketch):

```rust
#[derive(Debug, Subcommand)]
pub enum Command {
    // ... existing variants (KitePeriodicity, Simulate, SimulateGrid, ...) ...
    Synth(SynthArgs),
    SynthBatch(SynthBatchArgs),
    SynthValidate(SynthValidateArgs),
    SynthSchema,
}

#[derive(Debug, Args)]
pub struct SynthArgs {
    /// YAML config file.
    pub config: PathBuf,
    /// Output prefix (PREFIX.fa, PREFIX.truth.tsv, PREFIX.periods.tsv).
    /// REQUIRED. There is no fallback to `global.output` in the YAML
    /// (that field is silently ignored in MVP — see §0 A3).
    #[arg(short, long, required = true)]
    pub out: PathBuf,
    /// Override the YAML's `seed`.
    #[arg(long)]
    pub seed: Option<u64>,
    /// Also emit PREFIX.diagnostics.json. CLI-only — there is no
    /// equivalent YAML field (see §0 A2).
    #[arg(long)]
    pub diagnostics: bool,
}

#[derive(Debug, Args)]
pub struct SynthBatchArgs {
    #[arg(long)] pub config_dir: PathBuf,
    #[arg(long)] pub out_dir:    PathBuf,
    #[arg(long, default_value_t = 0)] pub seed_offset: u64,
    #[arg(long)] pub diagnostics: bool,
}
```

Batch mode iterates every `*.yaml` under `config_dir` and writes to
`out_dir/<stem>.{fa,truth.tsv,periods.tsv}`. Per-config seed =
`top_seed XOR fnv1a(filename) + seed_offset`. Parallelise with rayon
over configs.

## 8. Test config corpus

The corpus is **18 conceptual tests** from taxonomy §5.4, realised as
**22 active YAML fixtures + 1 deferred placeholder** in
`tests/synth/configs/`:

```
T01_simple_tr.yaml
T02_simple_tr_indel.yaml
T03_wobble_aperiodic.yaml
T04_wobble_periodic.yaml
T05_hor_clean.yaml
T06_regime_A.yaml
T07_regime_C.yaml
T08a_div_0.00.yaml \
T08b_div_0.05.yaml  \
T08c_div_0.15.yaml   |  divergence sweep at k=4 — six fixtures, one
T08d_div_0.35.yaml   |  conceptual test (T08)
T08e_div_0.55.yaml  /
T08f_div_0.70.yaml /
T09_nested_hor.deferred.yaml     # NOT generated — placeholder doc
                                 # (nested HOR deferred per §0 Q3)
T10_phase_shift.yaml             # uses A5 single-HOR-plus-event shape
T11_insertion.yaml
T12_inversion.yaml               # uses A5 single-HOR-plus-event shape
T13_coexisting_periods.yaml
T14_mixed_families.yaml
T15_stratification.yaml
T16_hybrid.yaml
T17_random_negative.yaml
T18_at_rich.yaml
```

Keep each config ≤300 kb so the full batch runs in seconds. The
validator skips any file whose extension is `.deferred.yaml`; T09 is
checked in with a header comment explaining the deferral so it
remains visible as future work without breaking `synth-batch`.

## 9. Milestones & acceptance gates

Milestone numbering matches upstream §9, with kitehor-side gates.

| # | Milestone | Configs passing | Cargo deliverable |
|---|---|---|---|
| M1 | Schema + validator | — | `synth-validate` accepts/rejects expected configs; `synth-schema --print` matches `docs/new/simulator_schema.json` byte-for-byte |
| M2 | Templates + HOR + SIMPLE_TR blocks (no noise, no events) | T01 simple-TR, T05 clean-HOR, T13 coexisting periods | Generates correct length & structure (manual check: `k × n_copies × monomer_length`) |
| M3 | Final noise pass | T01, T02 | Realised mutation/indel counts within 5% of requested rates |
| M4 | SHIFT + INSERTION + periods.tsv | T10, T11 | periods.tsv contains true base, true HOR-unit, near-miss, harmonic; SHIFT round-trips through coord_map |
| M5 | Wobble (residual accumulator, both models) | T03, T04 | std of realised δ within 10% of `amplitude_bp`; FFT recovers `period_rows` to ±10% |
| M6 | Post-generation events (HYBRID, INVERSION, DUP, DEL) | T12, T16; `events_json` populated | HYBRID at_copy=27 slot=4 resolves to the expected bp range; INV preserves array length |
| M7 | Diagnostics + batch + full corpus | T01–T18 except T09 (deferred) | `kitehor synth-batch --config-dir tests/synth/configs --out-dir /tmp/out` produces **22 valid bundles** in <30 s; T09 `.deferred.yaml` placeholder is skipped |

One PR per milestone. CI gate per milestone: `cargo test --release
synth::` plus a smoke `synth-batch` of that milestone's configs.

## 10. Testing strategy

Test pyramid mirrors upstream §8.

**Unit tests** (`#[cfg(test)] mod tests` per module):

- `templates`: realised divergence within tolerance over
  `k ∈ {3, 6, 12}`, `monomer_len ∈ {50, 171, 500}`, sampling 1000+
  base pairs of pairwise comparison.
- `wobble`: realised amplitude within 10% of target on a 1 Mb dummy
  sequence; periodicity recovery via rustfft (already a crate dep).
- `coords`: logical→realised→logical round-trip identity for a
  representative config with wobble + post-gen events.
- `events::INVERSION`: `rc(b"ACGT") == b"ACGT"` (palindrome), `rc(b"AAAA")
  == b"TTTT"` (one-base test catches reverse-only or complement-only
  bugs).
- `noise`: counts within statistical tolerance for L=100 kb, p=0.02,
  using a fixed seed and a Monte Carlo over 32 seeds for the rare-event
  tail.

**Integration tests** (`tests/synth_*.rs`):

- Each T-config generates without error.
- FASTA length matches expected (block sum ± realised indels ± noise).
- `truth.tsv` columns populated; hardcode expected values for T01,
  T05, T10.

**Determinism**:

- `(config, seed)` → byte-identical FASTA and TSV across two runs
  (compare via `sha256sum`).
- Different seeds → different FASTA but truth-property fields within
  tolerance.

**Property-level cross-check** (the bookkeeping-bug catcher):

For each T-config, recompute properties from the realised FASTA by an
independent path (count monomers by direct match against template,
detect indels by walking the sequence) and verify they match
`truth.tsv`. Coordinate-bookkeeping bugs hide here and nowhere else.

## 11. Open questions

Resolved. See §0 (Decisions on the open questions) for the canonical
answers to Q1–Q8 and the six amendments (A1–A6) arising from review.
This section is intentionally left as a pointer.

## 12. Hand-off criteria — "simulator done, start detector"

Per upstream §12 plus kitehor-specific:

1. All 18 T-configs (or 17 if Q3 defers nested HOR) generate without
   error via `kitehor synth-batch`.
2. `cargo test --release` is green for `tests/synth_*` and `synth::*`.
3. Determinism: re-running `synth-batch` produces byte-identical
   outputs across two consecutive runs.
4. `truth.tsv` schema matches taxonomy §5.2 verbatim — manually
   verified for T01, T05, T10.
5. `docs/new/simulator_schema.json` is the source of truth: CI checks
   `kitehor synth-schema --print` equals the file content.
6. A `kitehor synth-batch --config-dir tests/synth/configs --out-dir
   ground_truth_v2/` run completes in <30 s on the container (6 CPUs,
   64 GB).
7. The generated `truth.tsv` files are the input to the upcoming
   detector M1–M6.

## 13. Out of scope for this plan

- Detector implementation (`detect_spec.md` M1–M6) — separate plan.
- Migration of the existing rule.rs evaluation flow off `params.tsv`
  to YAML configs — separate effort, do not block detector work.
- Visual diagnostics (PNG rasters, FFT plots) — `detect_spec.md` §5
  marks these as detector outputs, not simulator.
- Real biological FASTA in `tests/synth/configs/` — kitehor data
  policy forbids it.
