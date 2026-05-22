# kitehor

Sequence-agnostic Higher-Order Repeat (HOR) detector for tandem-repeat
array FASTA sequences. No supplied monomer consensus, no species
library, no neural network — just kite k-mer periodograms and a
deterministic rule-based pipeline.

For each record, `kitehor analyze` reports whether the array contains a
HOR, its multiplicity (`k`), founder + tile periods, whether the
monomer carries an SSR or internal sub-repeat, and a single
`combined_class` summarising the structural call.

## Install

**Conda** (linux-64, from the `petrnovak` channel):

```bash
conda install -c petrnovak -c conda-forge kitehor
# or with mamba:
mamba install -c petrnovak -c conda-forge kitehor
```

**Pre-built binary** (linux-64): grab the tarball from the
[latest GitHub release](https://github.com/kavonrtep/kitehor/releases/latest),
extract, and put `kitehor` on your `PATH`.

**From source** (any platform with Rust ≥ 1.85):

```bash
git clone https://github.com/kavonrtep/kitehor && cd kitehor
cargo build --release
# binary at target/release/kitehor
```

Contributors building from source should also install the tracked
git hooks so commits/pushes can't ship fmt / clippy / test
regressions:

```bash
git config core.hooksPath .githooks    # one-time per clone
```

See [`docs/release.md`](docs/release.md) for the release runbook.

## Quick start

```bash
# End-to-end on one FASTA — always emits all 8 per-stage TSVs.
kitehor analyze input.fasta -o out_prefix
cat out_prefix.summary.tsv      # 32 columns; last is combined_class
```

`combined_class` is one of: `hor`, `hor_with_ssr`, `tr`, `tr_with_ssr`,
`tr_with_nested_tr`, `tr_with_subrepeat`, `pure_ssr`, `unresolved`.

### Pipeline at a glance

```
FASTA
  │
  ▼
kite-periodicity       k-mer pair-distance periodogram
  │
  ▼
rule-classify          HOR / simple_tr / unresolved
  │
  ├──────────────┬──────────────┐
  ▼              ▼              ▼
subrepeat-scan   ssr-scan       hor-validate
(nested-TR)      (short motifs) (within-tile + density)
  │              │              │
  └──────────────┴──────────────┘
                 ▼
            summary-merge      combined_class
```

Every stage is also exposed as a standalone subcommand for debugging
and partial reruns. Full algorithm + flag reference:
**[`docs/rule_proto.md`](docs/rule_proto.md)**.

## Smoke test

A 3-record synthetic fixture (87 KB) ships under `test_data/smoke/`:

```bash
./target/release/kitehor analyze test_data/smoke/sequences.fasta \
    -o /tmp/smoke
# Expected combined_class: hor_k3 → hor, hor_k5 → hor, tandem_pure → tr
```

## Subcommands

| Command | What it does | Detailed docs |
|---|---|---|
| `analyze` | End-to-end pipeline; writes all 8 per-stage TSVs | [`docs/rule_proto.md`](docs/rule_proto.md) |
| `kite-periodicity` | k-mer-distance periodogram (Rust port of TideCluster's `kite.R`) | [`docs/rule_proto.md`](docs/rule_proto.md) |
| `rule-classify` | HOR / simple_tr / unresolved verdict per record | [`docs/rule_proto.md`](docs/rule_proto.md#rule-classify) |
| `subrepeat-scan` | Spatial alternation / nested-TR detector | [`docs/rule_proto.md`](docs/rule_proto.md#subrepeat-scan) |
| `ssr-scan` | TideCluster-style SSR + kite-driven consensus | [`docs/rule_proto.md`](docs/rule_proto.md#ssr-scan) |
| `hor-validate` | Within-tile + spatial density HOR validation | [`docs/rule_proto.md`](docs/rule_proto.md#hor-validate) |
| `summary-merge` | Outer-join + 8-rule combined_class | [`docs/rule_proto.md`](docs/rule_proto.md#summary-merge) |
| `detect` / `detect-batch` | v2 line-width detector (operates on row embeddings) | [`docs/new/detect_impl_plan.md`](docs/new/detect_impl_plan.md) |
| `simulate` / `simulate-grid` | Legacy params.tsv-driven simulator | — |
| `synth*` | v2 YAML-driven simulator (wobble, phase shifts, events) | [`docs/new/simulator_impl_plan.md`](docs/new/simulator_impl_plan.md) |

Run any subcommand with `--help` for the full flag list.

### `kite-periodicity --classify`

Single-stage variant for users who only need a HOR verdict (skips the
SSR/subrepeat/within-tile checks). Output adds `verdict`, `founder`,
`multiplicity`, `tile`, `share` columns.

```bash
kitehor kite-periodicity input.fasta -o predictions.tsv --classify
```

Same classifier as `rule-classify`; convenience for single-step usage.

## v2 line-width detector

A separate, sequence-agnostic threshold-rule classifier that operates
on row embeddings of the input FASTA. Reads a v2-schema
`periods.tsv` and writes a fixed property bundle per record
(`.properties.tsv`, `.segments.tsv`, `.width_features.tsv`,
`.diagnostics.json`, `.consensus.fa`).

```bash
# Auto-mode (kite + rule classifier run internally)
kitehor detect input.fasta -o out_prefix

# With explicit period candidates
kitehor detect input.fasta --periods periods.tsv -o out_prefix
```

Classes: `HOR`, `irregular_HOR`, `simple_TR`, `mixed`, `ambiguous`.
Design + acceptance gates:
[`docs/new/detect_impl_plan.md`](docs/new/detect_impl_plan.md).

## v2 simulator (`synth*`)

YAML-driven structural simulator covering the full tandem-repeat
taxonomy (HORs, wobble, phase shifts, insertions, hybrid monomers,
inversions, dup/del events).

```bash
kitehor synth tests/synth_configs/T05_hor_clean.yaml -o /tmp/t05
kitehor synth-batch --config-dir tests/synth_configs --out-dir /tmp/synth_corpus
kitehor synth-validate path/to/config.yaml
kitehor synth-schema > simulator_schema.json
```

Design + canonical YAML schema:
[`docs/new/simulator_impl_plan.md`](docs/new/simulator_impl_plan.md),
[`docs/new/simulator_schema.json`](docs/new/simulator_schema.json).

## Test data

This repo intentionally ships **no real biological FASTA**:

- `test_data/smoke/` — 3-record synthetic fixture for build verification.
- `test_data/ci_corpus/` — 13-record curated corpus exercising 5/8
  `combined_class` values; provenance in
  [`test_data/ci_corpus/manifest.tsv`](test_data/ci_corpus/manifest.tsv).
- `tests/synth_configs/` — 23 v2-simulator CI fixtures (T01–T20).

Larger benchmarks (`ground_truth/`, `ground_truth_v2/`) ship only
their specs; bundles are regenerated on demand with
`simulate-grid` / `synth-batch`.

## Layout

```
src/
  analyze.rs         end-to-end orchestrator
  rule_classify/     HOR / simple_tr / unresolved classifier
  subrepeat/         nested-TR detector
  ssr/               TideCluster SSR + kite consensus
  hor_validate/      within-tile + density validation
  summary/           combined_class merger
  kite.rs            k-mer periodogram
  emit_periods.rs    bridge: kite output → v2 detector periods.tsv
  detect/            v2 line-width detector
  simulate*.rs       legacy params.tsv-driven simulator
  synth/             v2 YAML-driven simulator
tools/rule_proto/    Python reference (validation oracle for the port)
docs/                project docs — see docs/rule_proto.md for the pipeline
test_data/           shipped fixtures
ground_truth*/       benchmark corpus specs (bundles regenerated on demand)
tests/               integration tests
```

## License

Dual licensed under [MIT](LICENSE-MIT) or
[Apache-2.0](LICENSE-APACHE), at your option.
