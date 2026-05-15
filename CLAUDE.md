# kitehor — Claude dev guide

Local-only notes for working in this repo with Claude Code.

## What this is

Sequence-agnostic HOR detector. Workhorse subcommand:

```
kitehor kite-periodicity <fasta> -o out.tsv --classify
```

That runs kite → features → RF → Platt → 4-verdict logic + k-recovery.

## Repo layout shortcut

```
src/                  Rust crate (lib + bin)
config/classifier.toml   ← model thresholds, Platt coefs, baked into binary
models/               Random-forest JSON dumps (baked into binary via include_bytes!)
tools/training/       R training pipeline + model exporter
tools/features/       Python reference feature extractors
ground_truth/         params.tsv + simulator helpers; sequences are regenerated
test_data/smoke/      87 KB synthetic fixture for build verification
examples/             validate_rf — diff Rust vs. an R reference TSV
conda/kitehor/        conda recipe (meta.yaml; built by .github/workflows/conda-release.yml)
.github/workflows/    ci.yml, release.yml, conda-release.yml
docs/                 project docs — see docs/ci-status.md for the CI/release plan
docs/archive/         (gitignored) historical design docs
```

## Project docs

All non-README, non-CLAUDE documentation lives in `docs/`. Start with
[`docs/ci-status.md`](docs/ci-status.md) for the CI/CD plan, locked
decisions, and the release runbook. Add new topic docs as
`docs/<topic>.md` siblings; do not scatter markdown at the repo root.

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
