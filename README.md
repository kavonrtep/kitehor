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

### v2 simulator (`synth*`)

`kitehor synth` is a richer, YAML-driven simulator built to feed the
upcoming line-width detector described in
[`docs/new/detect_spec.md`](docs/new/detect_spec.md). Unlike the
`simulate`/`simulate-grid` pair (which targets the rule-based training
corpus), `synth` can express the full structural taxonomy:

- Simple TRs and HORs of arbitrary multiplicity, with inter-slot
  divergence as a tunable parameter.
- Continuous wobble (sinusoidal or AR(1) random walk), realised via a
  residual-accumulator integer-edit scheme.
- Discrete phase shifts, non-tandem insertions
  (random / AT_rich / GC_rich / retro_like / segdup_like), and
  post-generation events (HYBRID, INVERSION, DUPLICATION, DELETION).

Outputs per array: `{prefix}.fa`, `{prefix}.truth.tsv` (property
vector + structural-expression string + events_json), and
`{prefix}.periods.tsv` (period candidates simulating an upstream
generator). Add `--diagnostics` for a `{prefix}.diagnostics.json` with
per-stage provenance.

The canonical schema lives at
[`docs/new/simulator_schema.json`](docs/new/simulator_schema.json) and
is embedded into the binary; `kitehor synth-schema` dumps it to stdout
and a CI test catches drift. A 22-fixture test corpus (T01–T18, with
T08 a six-point divergence sweep and T09 deferred) is under
`tests/synth_configs/`. See
[`docs/new/simulator_impl_plan.md`](docs/new/simulator_impl_plan.md)
for the implementation contract.

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
config/classifier.toml Thresholds, Platt coefs, imputation medians
models/                Random-forest JSON (HOR-score + k-predictor)
tools/training/        R scripts: ranger train, Platt fit, model export, CV
tools/features/        Python reference feature extractors
ground_truth/          Synthetic-benchmark spec (regenerate sequences from params.tsv)
test_data/smoke/       Tiny fixture for build verification
examples/              Cross-validation harness vs. the R prototype
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
