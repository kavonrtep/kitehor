# kitehor — Claude dev guide

Local-only notes for working in this repo with Claude Code.

## What this is

Sequence-agnostic HOR detector. Three top-level workflows:

1. **Rule-proto pipeline** (workhorse, port of
   `tools/rule_proto/*.py` — see `docs/rule_proto.md`):
   ```
   # End-to-end
   kitehor analyze <fasta> -o <prefix>

   # Or per stage:
   kitehor kite-periodicity <fasta> -o out.tsv --classify
   kitehor rule-classify <peaks.tsv> -o <prefix>
   kitehor subrepeat-scan <fasta> --kite-peaks <peaks.tsv> -o <prefix>
   kitehor ssr-scan <fasta> --kite-peaks <peaks.tsv> -o <prefix>
   kitehor hor-validate <fasta> --verdicts <v.tsv> --global-peaks <p.tsv> -o <prefix>
   kitehor summary-merge --verdicts ... --subrepeat ... --ssr ... -o <prefix>
   ```
   `analyze` always emits all 8 per-stage TSVs (debugging contract).
   The 8 combined_class values: `pure_ssr, tr_with_nested_tr,
   tr_with_subrepeat, hor_with_ssr, hor, tr_with_ssr, tr, unresolved`.
   Replaces the older 4-condition rule (`src/rule.rs`) and the legacy
   ML pipeline (both removed in P1 / P6 of the port).

3. **v2 line-width detector** (`docs/new/detect_spec.md`):
   ```
   kitehor detect <fasta> -o <prefix>                       # auto-periods (default)
   kitehor detect <fasta> --periods <periods.tsv> -o <prefix>  # explicit periods
   ```
   With `--periods` omitted, the detector runs `kite-periodicity` +
   the rule classifier internally with their defaults, persists the
   derived periods to `<prefix>.periods.tsv`, and implicitly applies
   `--allow-missing-periods` so QC-rejected records end up
   classified as `ambiguous`. Score mapping is the same as
   `kite-periodicity --classify --emit-periods` (founder=0.95,
   tile=0.90, other top-3=0.60; `Unresolved` hints at
   0.50/0.40/0.30). For tuned kite parameters (kmer size, score2
   threshold, rule top-N, …), run the explicit two-step pipeline
   (§3) and pass the resulting `periods.tsv` via `--periods`.

   Bundle outputs (both modes): `<prefix>.properties.tsv` /
   `.segments.tsv` (13 columns post-M7.2; includes
   `consensus_identity_to_reference` + `_coverage`) /
   `.width_features.tsv` / `.diagnostics.json` (`schema_version=2`
   post-M7.3) / `.consensus.fa` (per-segment monomers for class=mixed,
   whole-array monomer + optional hor_unit for resolved classes).
   Auto-mode additionally writes `<prefix>.periods.tsv`.

   Calibration baseline: 94.4% per-class accuracy oracle, 90.6% kite
   periods (1600-case `ground_truth_v2/`).

4. **Explicit two-step kite → detector pipeline** (useful when you
   need to tune kite or inspect intermediate output):
   ```
   kitehor kite-periodicity <fasta> -o /tmp/kite.tsv --classify \
       --emit-periods /tmp/kite.periods.tsv
   kitehor detect <fasta> --periods /tmp/kite.periods.tsv -o /tmp/det \
       --allow-missing-periods
   ```
   Produces byte-identical detector output to the auto-mode form
   above when kite is run with its defaults. `--allow-missing-periods`
   is needed when kite produces no rows for a record — either because
   the rule returned `NoSignal`, or because the FASTA record was
   rejected by kite's `LoadQC` (e.g., too short, too many Ns). In
   both cases the detector tags the record `ambiguous`.

The legacy ML pipeline (RF + Platt + k-recovery + homology) and its
CLI flags (`--use-ml-classifier`, `--no-hor-call`, `--hor-*`,
`--coverage`, etc.) were removed in P6 of the rule-proto port.

## Repo layout shortcut

```
src/                    Rust crate (lib + bin)
  analyze.rs            ← end-to-end rule-proto pipeline orchestrator
  rule_classify/        ← HOR / simple_tr / unresolved classifier
                          (port of tools/rule_proto/rule_proto.py)
  subrepeat/            ← nested-TR detector (subrepeat_scan.py)
  ssr/                  ← TideCluster SSR scan + consensus (ssr_scan.py)
  hor_validate/         ← within-tile + density (hor_within_tile_check.py)
  summary/              ← 8-rule combined_class merger (summary.py)
  kite.rs               k-mer periodogram (the upstream stage)
  emit_periods.rs       bridge: kite output → v2 detector periods.tsv
  detect/               ← v2 line-width detector (`kitehor detect*`)
  simulate*.rs          legacy params.tsv-driven simulator
  synth/                ← v2 YAML-driven simulator (`kitehor synth*`)
  monomer_model.rs      `probe_period` helper (currently unused;
                          retained for potential future use)
tools/rule_proto/       Python prototype kept as the reference oracle
                          (validation target; not invoked at runtime)
ground_truth/           legacy params.tsv + simulator helpers
ground_truth_v2/        ← v2 corpus spec (1600 configs) + run_batch.sh
test_data/smoke/        87 KB synthetic fixture for build verification
test_data/ci_corpus/    diverse small corpus from prototype run (P7)
tests/synth_configs/    ← v2 CI fixtures (T01–T20; 23 active + 1 deferred)
conda/kitehor/          conda recipe
.github/workflows/      ci.yml, release.yml, conda-release.yml
docs/                 project docs
  rule_proto.md       ← the rule-proto pipeline (current default)
  rule.md             archived — the older 4-condition rule (P1 retirement)
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
      -o /tmp/smoke.tsv --classify
  ```
  Expect: `tandem_pure`→simple_tr, `hor_k3`→hor (k=3, founder=100, tile=300),
  `hor_k5`→hor (k=5, founder=150, tile=750).

- **End-to-end smoke** (rule-proto pipeline):
  ```
  ./target/release/kitehor analyze test_data/smoke/sequences.fasta \
      -o /tmp/smoke
  ```
  Writes 8 per-stage TSVs at `/tmp/smoke.*.tsv` plus `.summary.tsv`.
  `combined_class` column on the summary: hor / hor / tr.

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

- **Auto-mode smoke** (one-step kite + detect):
  ```
  ./target/release/kitehor detect \
      test_data/smoke/sequences.fasta -o /tmp/smoke.det
  ```
  Internally runs `kite-periodicity --classify` with kite defaults,
  derives periods, persists them to `/tmp/smoke.det.periods.tsv`,
  and runs the detector. Reads `/tmp/smoke.det.properties.tsv` for
  per-record classes.

- **Explicit two-step smoke** (when tuning kite knobs is needed):
  ```
  ./target/release/kitehor kite-periodicity \
      test_data/smoke/sequences.fasta -o /tmp/smoke.kite.tsv \
      --classify --emit-periods /tmp/smoke.kite.periods.tsv
  ./target/release/kitehor detect \
      test_data/smoke/sequences.fasta \
      --periods /tmp/smoke.kite.periods.tsv -o /tmp/smoke.det \
      --allow-missing-periods
  ```
  Byte-equivalent to auto-mode when kite is run with defaults.

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
