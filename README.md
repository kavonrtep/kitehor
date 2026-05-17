# kitehor

De novo Higher-Order Repeat (HOR) detector for tandem-repeat array FASTA
sequences. Sequence-agnostic — no monomer consensus, no species library,
no neural network.

`kitehor` takes a FASTA of pre-extracted tandem arrays and reports, per
record, whether the array contains a HOR, its multiplicity (`k`), founder
period, tile period, and a calibrated probability score.

## Pipeline

```
FASTA → kite k-mer-distance histogram → top peaks (d1..d10)
      → feature row (kite peaks + diversity + family + homology)
      → random-forest HOR-probability + Platt calibration
      → 4-category verdict: hor / tandem / unresolved / no_signal
      → k-recovery from k-predictor when family-fit is ambiguous
```

The kite stage is a Rust port of the periodicity detector in
`TideCluster/tarean/kite.R`. The classifier is trained on synthetic data
(see `tools/training/`) and ships as JSON-serialised random forests in
`models/`.

## Build

Requires Rust ≥ 1.75.

```bash
cargo build --release
```

The resulting binary is `target/release/kitehor`.

## Smoke test

A 3-record synthetic fixture (87 KB) is shipped under `test_data/smoke/`:
one pure tandem and two HORs (k=3, k=5). To verify the build:

```bash
./target/release/kitehor kite-periodicity \
    test_data/smoke/sequences.fasta \
    -o /tmp/smoke.tsv \
    --classify --no-hor-call

# Expected verdicts under the default (rule) classifier:
#   tandem_pure    -> tandem
#   hor_k3         -> hor   (founder=100, k=3, tile=300, share~0.45)
#   hor_k5         -> hor   (founder=150, k=5, tile=750, share~0.80)
```

The fixture itself was produced from `test_data/smoke/params.tsv` via
`kitehor simulate-grid` and is fully deterministic; you can regenerate
it at any time:

```bash
./target/release/kitehor simulate-grid \
    --params test_data/smoke/params.tsv \
    --outdir test_data/smoke --seed 42
```

## Usage

```bash
# Periodicity scan only (no classifier).
kitehor kite-periodicity input.fasta -o periods.tsv

# Periodicity + HOR call (rule-based, default).
kitehor kite-periodicity input.fasta -o predictions.tsv --classify

# Periodicity + legacy ML classifier (opt-in).
kitehor kite-periodicity input.fasta -o predictions.tsv --classify --use-ml-classifier

# v2 line-width detector (sequence-agnostic; consumes a periods.tsv).
kitehor detect input.fasta --periods periods.tsv -o out_prefix
kitehor detect-batch \
    --fasta-dir corpus/ --periods-dir corpus/ --out-dir det_out/

# Combined pipeline: kite candidates → v2 detector.
kitehor kite-periodicity input.fasta -o predictions.tsv --classify \
    --emit-periods kite.periods.tsv
kitehor detect input.fasta --periods kite.periods.tsv -o det_out \
    --allow-missing-periods   # if kite returned NoSignal for any record

# Simulate one synthetic array.
kitehor simulate --monomer-size 171 --multiplicity 12 --copies 100 -o sim.fa

# Simulate the full ground-truth grid (1,600 cases) from a params TSV.
kitehor simulate-grid --params ground_truth/params.tsv --outdir out/sim

# v2 simulator: YAML-driven, structurally rich (wobble, phase shifts,
# insertions, hybrid monomers, inversions, dup/del events). See
# docs/new/ for the design contract.
kitehor synth tests/synth_configs/T05_hor_clean.yaml -o /tmp/t05
kitehor synth-batch \
    --config-dir tests/synth_configs --out-dir /tmp/synth_corpus
kitehor synth-validate path/to/config.yaml
kitehor synth-schema > simulator_schema.json
```

### Default `--classify` columns (rule-based)

| Column | Meaning |
|---|---|
| `verdict`       | `hor` / `tandem` / `unresolved` / `no_signal` |
| `founder`       | Inferred founder period (bp), only for `hor` |
| `multiplicity`  | k (1 for tandem) |
| `tile`          | HOR tile period (bp) |
| `share`         | `min(s_founder, s_tile) / max(...)` — diagnostic only |

The rule is documented in detail in [docs/rule.md](docs/rule.md). It
trusts kite peak detection (every kite peak has already passed the
`peak > background` filter) and calls HOR when `d1` is a `k≥2`
integer multiple of a top-3 kite peak within tolerance.

### Supplementary coverage QC (`--coverage`)

Pass `--coverage` together with `--classify` to add nine columns
quantifying how well the rule's `tile` period actually tiles the
array. Levenshtein identity between the first tile and each
subsequent tile-aligned window:

```bash
kitehor kite-periodicity input.fasta -o predictions.tsv --classify --coverage
```

Adds: `cov_mean`, `cov_pass_70`, `cov_pass_80`, `cov_pass_90`,
`cov_first_half`, `cov_second_half`, `cov_min`, `cov_max`,
`cov_n_tiles`. Non-HOR rows get NA. Supplementary only — it does not
enter the HOR decision. See [docs/rule.md](docs/rule.md) for what the
score does and does not catch.

### v2 line-width detector (`kitehor detect*`)

A sequence-agnostic, threshold-rule classifier that operates on row
embeddings of the input FASTA. Reads a `periods.tsv` (v2 schema:
`array_id\tperiod_bp\tperiod_score\tsource`) and writes a fixed
property bundle per record:

| Output | Contents |
|---|---|
| `{prefix}.properties.tsv`     | Per-record class + base width, k, IC, phase_sep, wobble, n_phase_shifts, irregularity, confidence |
| `{prefix}.width_features.tsv` | One row per tested width with R(k), IC, edge rates, class hint |
| `{prefix}.segments.tsv`       | Per-segment rows when n_phase_shifts > 0 |
| `{prefix}.diagnostics.json`   | Structured per-array reason + all of the above |
| `{prefix}.consensus.fa`       | Monomer + HOR-unit consensuses (only for resolved classes) |

Classes: `HOR`, `irregular_HOR`, `simple_TR`, `mixed`, `ambiguous`.
M6 calibration baseline: **94.4%** per-class accuracy on the
1600-case `ground_truth_v2/` benchmark (≥ 92% target met). Design
contract: [`docs/new/detect_impl_plan.md`](docs/new/detect_impl_plan.md).

### Combined pipeline (kite → detector)

`kite-periodicity --emit-periods` writes a v2-compatible
`periods.tsv` so kite candidates can drive the line-width detector
in one shell pipeline. Score mapping
(`src/emit_periods.rs`):

| Rule verdict | Rows written |
|---|---|
| `Hor{founder, tile}` | founder @ 0.95 (`kite_founder`); tile @ 0.90 (`kite_tile`) if distinct; remaining top-3 peaks @ 0.60 (`kite_secondary`) |
| `Tandem{monomer}`    | monomer @ 0.95 (`kite_monomer`); remaining top-3 peaks @ 0.60 (`kite_secondary`) |
| `Unresolved`         | top-3 peaks @ 0.50 / 0.40 / 0.30 (`kite_peak`) |
| `NoSignal`           | no rows — pass `--allow-missing-periods` to detect |
| QC-skipped record    | no rows — pass `--allow-missing-periods` to detect |
| no `--classify`      | top-3 peaks @ 0.60 (`kite_peak`) |

The emitter never looks past Kite's top-3 peaks — anything at rank
4 or deeper is discarded regardless of verdict. Mutually exclusive
with `--use-ml-classifier`.

Scores are chosen relative to the detector's
`strong_period_score = 0.85` gate: values ≥ that floor can fire HOR
rescue paths; below-floor scores act as hints to the canonical
column-IC test only.

### Legacy ML classifier (`--use-ml-classifier`)

A random-forest + Platt-scaled classifier with k-recovery and homology
features. Output adds: `hor_score`, `hor_score_raw`, `k_pred`,
`recovered`, `h_d1`, `h_founder`. Useful when working with data drawn
from the same distribution as the synthetic training set; otherwise
the rule-based default is more reliable on real centromeric arrays
(see [docs/rule.md](docs/rule.md) for the empirical comparison).

The legacy thresholds (`t_low = 0.15`, `t_high = 0.71`) live in
`config/classifier.toml`. Override them or the Platt coefficients with
`--classifier-config <path.toml>`.

## v2 simulator (`synth*`)

A richer, YAML-driven simulator that coexists with
`simulate`/`simulate-grid`. Built to feed the upcoming line-width
detector ([`docs/new/detect_spec.md`](docs/new/detect_spec.md)),
`synth` expresses the full structural taxonomy: arbitrary-`k` HORs
with tunable inter-slot divergence; continuous wobble (sinusoidal or
AR(1) random walk via residual-accumulator integer edits); discrete
phase shifts; non-tandem insertions (`random`/`AT_rich`/`GC_rich`/
`retro_like`/`segdup_like`); post-generation events (HYBRID,
INVERSION, DUPLICATION, DELETION).

```bash
# One array.
kitehor synth tests/synth_configs/T05_hor_clean.yaml -o /tmp/t05

# Whole corpus in parallel.
kitehor synth-batch \
    --config-dir tests/synth_configs --out-dir /tmp/synth_corpus

# Schema-validate a config without generating.
kitehor synth-validate path/to/config.yaml

# Dump the canonical JSON Schema.
kitehor synth-schema > simulator_schema.json
```

Per-array outputs: `{prefix}.fa`, `{prefix}.truth.tsv` (property
vector + structural-expression + events_json), `{prefix}.periods.tsv`
(period candidates as an upstream generator would emit). Add
`--diagnostics` for `{prefix}.diagnostics.json` with per-stage
provenance (RNG sub-stream seeds, realised template slots,
per-block coordinates, wobble/noise counts).

The canonical YAML schema lives at
[`docs/new/simulator_schema.json`](docs/new/simulator_schema.json) and
is embedded into the binary (drift-tested in CI). Two corpora ship:

- **`tests/synth_configs/`** — 23-fixture CI corpus covering T01–T20
  (T08 a six-point divergence sweep; T09 nested-HOR is deferred).
- **`ground_truth_v2/`** — 1,600-case benchmark corpus across 9
  categories (simple_TR, hor_clean, hor_wobble, hor_shift,
  hor_insertion, hor_event, mixed, random, gc_bias). Spec only;
  generated FASTA bundles are gitignored. Re-run with
  `./ground_truth_v2/run_batch.sh` (~2 s wall on 16 CPUs).

Implementation contract + acceptance gates:
[`docs/new/simulator_impl_plan.md`](docs/new/simulator_impl_plan.md).

## Test data

This repo intentionally ships **no real biological FASTA**. The only
sequence data included is the small `test_data/smoke/` fixture
generated by the simulator. To work with real data, point the binary
at your own FASTA. To benchmark on the synthetic corpus, regenerate it
from `ground_truth/params.tsv` (see below) — or train your own.

## Ground truth (regenerated on demand)

`ground_truth/` ships only `params.tsv` (the spec for the 1,600-case
synthetic benchmark). The actual sequences, truth labels, monomer
coordinates and event log are not in the repo — regenerate them
locally:

```bash
./target/release/kitehor simulate-grid \
    --params ground_truth/params.tsv \
    --outdir ground_truth \
    --seed <N>
```

Pick any `<N>`; the simulator is deterministic per `(N, case_id)`. The
classifier was trained on a snapshot of this corpus; verdicts on a
re-baselined version will not be byte-identical to historical results
but will agree at the population level.

## Layout

```
src/                   Rust crate
  rule.rs              default HOR classifier (4-condition rule)
  classifier.rs        legacy ML loader (RF + Platt); --use-ml-classifier only
  emit_periods.rs      kite → v2 detector periods.tsv bridge
  detect/              v2 line-width detector (kitehor detect*)
  simulate*.rs         legacy params.tsv-driven simulator
  synth/               v2 YAML-driven simulator (kitehor synth*)
config/classifier.toml Thresholds, Platt coefs, imputation medians
models/                Random-forest JSON (HOR-score + k-predictor)
tools/training/        R scripts: ranger train, Platt fit, model export, CV
tools/features/        Python reference feature extractors
ground_truth/          Legacy synthetic-benchmark spec (regenerate from params.tsv)
ground_truth_v2/       v2 benchmark corpus spec (1600 configs in 9 categories)
test_data/smoke/       Tiny fixture for build verification
tests/synth_configs/   v2 CI fixtures (T01–T20, 23 active + 1 deferred)
examples/              Cross-validation harness vs. the R prototype
docs/
  rule.md              rule classifier algorithm
  ci-status.md         CI/release plan + runbook
  new/                 v2 simulator + detector design docs
```

## Training a fresh model

```bash
# 1. Generate the synthetic corpus (1,600 cases).
./target/release/kitehor simulate-grid \
    --params ground_truth/params.tsv --outdir corpus --seed 42

# 2. Extract features per record (kite + diversity).
python3 tools/features/extract_features.py \
        --seed-dir corpus --out features.tsv

# 3. Add homology features via the Rust probe.
python3 tools/features/add_homology.py \
        --features features.tsv --fasta corpus/sequences.fasta \
        --out features_h.tsv

# 4. Train RF + Platt.
Rscript tools/training/hor_model_rf_h.R
Rscript tools/training/hor_model_k.R
Rscript tools/training/fit_platt.R

# 5. Export both forests to JSON for the Rust crate.
Rscript tools/training/export_ranger.R \
        --in eval/reports/hor_model_rf_h/ranger_model.rds \
        --out models/hor_score.rftrees.json
Rscript tools/training/export_ranger.R \
        --in eval/reports/hor_model_k/ranger_model_k.rds \
        --out models/k_pred.rftrees.json
```

## License

Dual licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE),
at your option.
