# kitehor — Claude dev guide

Local-only notes for working in this repo with Claude Code.

## What this is

Sequence-agnostic HOR detector. Two pipelines coexist:

1. **Rule classifier on kite peaks** (workhorse):
   ```
   kitehor kite-periodicity <fasta> -o out.tsv --classify
   ```
   Runs kite → 4-condition rule (`src/rule.rs`): HOR ⟺ d1 = k×p_n
   for k ∈ [2, 30] with p_n in top-3 by score. See `docs/rule.md`.

2. **v2 line-width detector** (`docs/new/detect_spec.md`):
   ```
   kitehor detect <fasta> --periods <periods.tsv> -o <prefix>
   ```
   Consumes a `periods.tsv` (v2 schema) and produces
   `<prefix>.properties.tsv` / `.segments.tsv` / `.width_features.tsv` /
   `.diagnostics.json` / `.consensus.fa`. Calibration baseline:
   94.4% per-class accuracy on the 1600-case `ground_truth_v2/` benchmark.

3. **Combined pipeline** (kite candidates → detector):
   ```
   kitehor kite-periodicity <fasta> -o /tmp/kite.tsv --classify \
       --emit-periods /tmp/kite.periods.tsv
   kitehor detect <fasta> --periods /tmp/kite.periods.tsv -o /tmp/det \
       --allow-missing-periods
   ```
   Score mapping in `src/emit_periods.rs`: founder=0.95,
   tile=0.90, other top-3=0.60; ambiguous verdicts (`Unresolved`)
   emit hints at 0.50/0.40/0.30. `--allow-missing-periods` is
   needed when kite produces no rows for a record — either because
   the rule returned `NoSignal`, or because the FASTA record was
   rejected by kite's `LoadQC` (e.g., too short, too many Ns). In
   both cases the detector will tag the record `ambiguous`.

The earlier ML pipeline (RF + Platt + k-recovery + homology) is still
available via `--use-ml-classifier`. It is over-sensitive on real
centromeric arrays and under-sensitive on real HORs with strong
inter-position divergence (the training-set distribution doesn't match
real data) — use only when the input is drawn from a similar
distribution to the synthetic training corpus.

## Repo layout shortcut

```
src/                  Rust crate (lib + bin)
  rule.rs             default HOR classifier (4-condition rule)
  classifier.rs       legacy ML loader (RF + Platt); used only under --use-ml-classifier
  classify.rs         legacy ML verdict orchestrator
  emit_periods.rs     bridge: kite output → v2 detector periods.tsv
  detect/             ← v2 line-width detector (`kitehor detect*`)
  simulate*.rs        legacy params.tsv-driven simulator (training corpus)
  synth/              ← v2 YAML-driven simulator (`kitehor synth*`)
config/classifier.toml legacy ML thresholds (only consulted with --use-ml-classifier)
models/               Legacy RF JSON (baked into binary; loaded only by ML path)
tools/training/       R training pipeline for the legacy model
tools/features/       Python reference feature extractors (for ML cross-check)
ground_truth/         legacy params.tsv + simulator helpers; sequences regenerated
ground_truth_v2/      ← v2 corpus spec (1600 configs in 9 categories) + run_batch.sh
test_data/smoke/      87 KB synthetic fixture for build verification
tests/synth_configs/  ← v2 CI fixtures (T01–T20; 23 active + 1 deferred)
examples/             validate_rf — legacy ML cross-check vs an R reference TSV
conda/kitehor/        conda recipe (meta.yaml; built by .github/workflows/conda-release.yml)
.github/workflows/    ci.yml, release.yml, conda-release.yml
docs/                 project docs
  rule.md             ← the rule classifier, current default
  ci-status.md        ← CI/release plan + runbook
  new/                ← v2 simulator + detector design docs
    taxonomy.md         structural taxonomy of tandem-repeat arrays
    detect_spec.md      line-width detector design (future work)
    simulator_plan.md   upstream simulator implementation plan
    simulator_impl_plan.md  kitehor-specific implementation contract
    simulator_schema.json   canonical YAML config schema
  archive/            (gitignored) historical design docs
```

## Project docs

All non-README, non-CLAUDE documentation lives in `docs/`. Start with
[`docs/ci-status.md`](docs/ci-status.md) for the CI/CD plan, locked
decisions, and the release runbook. Add new topic docs as
`docs/<topic>.md` siblings; do not scatter markdown at the repo root.
v2 simulator + detector design docs live in
[`docs/new/`](docs/new/) — read `simulator_impl_plan.md` §0 first for
the decisions snapshot and amendments table.

## v2 simulator (`synth*`)

Richer YAML-driven simulator that coexists with the legacy
`simulate`/`simulate-grid` pair. Lives entirely in `src/synth/`:

| Module | Purpose |
|---|---|
| `config.rs`      | YAML loader + serde structural validation + MVP business rules |
| `rng.rs`         | FNV-1a sub-stream derivation (matches parent project convention) |
| `templates.rs`   | HOR_slots / monomer instantiation; cached by template name |
| `coords.rs`      | `CoordMap` + `apply_indels` (kept-contiguous boundary rule) |
| `blocks.rs`      | HOR/SIMPLE_TR/SHIFT/INSERTION expansion |
| `wobble.rs`      | sinusoidal + AR(1) random_walk via residual-accumulator integer edits |
| `events.rs`      | HYBRID/INVERSION/DUPLICATION/DELETION + events_json |
| `noise.rs`       | final mutation + indel pass |
| `grammar.rs`     | `structural_expression` emission (taxonomy §2 grammar) |
| `truth.rs`       | `truth.tsv` writer + class inference |
| `periods.rs`     | period candidate generator (true_base + true_hor_unit + distractors) |
| `fasta.rs`       | FASTA writer |
| `diagnostics.rs` | optional `{prefix}.diagnostics.json` |
| `simulator.schema.json` | embedded canonical schema; drift-tested |

CLI: `synth`, `synth-batch`, `synth-validate`, `synth-schema`. See
[`docs/new/simulator_impl_plan.md`](docs/new/simulator_impl_plan.md)
for the implementation contract and milestone acceptance gates.

## Workflow

- **Build**: `cargo build --release`
- **Tests**: `cargo test --release`
- **Smoke**:
  ```
  ./target/release/kitehor kite-periodicity test_data/smoke/sequences.fasta \
      -o /tmp/smoke.tsv --classify --no-hor-call
  ```
  Expect: `tandem_pure`→tandem, `hor_k3`→hor (k=3, founder=100, tile=300),
  `hor_k5`→hor (k=5, founder=150, tile=750).

- **Full benchmark**: regenerate `ground_truth/sequences.fasta` from
  `ground_truth/params.tsv` (1,600 cases) before running the classifier
  on it — those files are not committed.

- **v2 simulator smoke** (23-fixture CI corpus):
  ```
  ./target/release/kitehor synth-batch \
      --config-dir tests/synth_configs --out-dir /tmp/synth_out
  ```
  Produces 23 bundles (`.fa` + `.truth.tsv` + `.periods.tsv`); the
  `T09_nested_hor.deferred.yaml` placeholder is skipped. Add
  `--diagnostics` for per-array `.diagnostics.json`.

- **v2 simulator full benchmark** (1600-case corpus): generated
  outputs are gitignored under `ground_truth_v2/out/`; regen with
  `./ground_truth_v2/run_batch.sh` (~2 s wall on 16 CPUs).

- **v2 detector** (consumes simulator-emitted periods):
  ```
  ./target/release/kitehor detect-batch \
      --fasta-dir ground_truth_v2/out \
      --periods-dir ground_truth_v2/out \
      --out-dir ground_truth_v2/det_out
  python3 tools/detect_eval/eval.py \
      --manifest ground_truth_v2/manifest.tsv \
      --properties-dir ground_truth_v2/det_out
  ```
  Calibration baseline 94.4% (M6 acceptance gate).

- **Combined pipeline smoke** (kite → emit-periods → detect):
  ```
  ./target/release/kitehor kite-periodicity \
      test_data/smoke/sequences.fasta -o /tmp/smoke.kite.tsv \
      --classify --emit-periods /tmp/smoke.kite.periods.tsv
  ./target/release/kitehor detect \
      test_data/smoke/sequences.fasta \
      --periods /tmp/smoke.kite.periods.tsv -o /tmp/smoke.det \
      --allow-missing-periods
  ```
  Reads `/tmp/smoke.det.properties.tsv` for per-record classes.

## Data policy

- **No real biological FASTA in the repo.** If you ever add real test
  data later, every record must carry provenance (assembly accession,
  contig, coordinates, sequencing project) — check it in alongside a
  `manifest.tsv` that documents source and how it was extracted.
- **No large simulator output in the repo.** `ground_truth/` keeps only
  the params spec; everything else is regenerated on demand.
- **Smoke fixture only** for quick build sanity. Anything bigger than a
  few hundred KB should live outside the repo.

## Engineering notes

- **`serde_json` float parsing**: must keep the `float_roundtrip`
  feature in `Cargo.toml`. Without it, RF split values parse 1 ULP off
  and tree traversals occasionally flip leaf direction, breaking
  bit-exact agreement with the R prototype. This bit us once already.
- **Determinism**: simulator seeds via FNV-1a of `master:case_id`; RF
  traversal is deterministic; Platt scaling is deterministic.
- **`monomer_model.rs`**: only `probe_period` is actively used. The
  larger inference module was retired with the rest of the pre-kite
  pipeline; if extending, prefer the classifier path in `classify.rs`.

## Validation history (pre-trim)

Before stripping the AT + synthetic FASTAs from the repo, the kitehor
port was validated end-to-end against the R prototype:

| Dataset | Records | Verdict agreement | Field agreement |
|---|---:|---:|---:|
| AT real centromeres | 155 | 155/155 | founder/k/tile/recovered all 155/155 |
| Synthetic CV seed 201 | 1,204 | 1,204/1,204 | all 1,204/1,204 |

Score difference vs. R: max \|Δ raw\| = 5×10⁻¹⁶ (bit-equivalent).
