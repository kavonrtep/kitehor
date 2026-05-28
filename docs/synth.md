# `kitehor synth` family — v2 YAML-driven simulator

YAML-driven structural simulator that covers the full tandem-repeat
taxonomy: HORs, sinusoidal + random-walk wobble, phase shifts /
insertions, hybrid monomers, inversions, and duplication / deletion
events. Used by:

- `tests/synth_configs/` — 23-fixture CI corpus (T01–T20)
- `ground_truth_v2/` — 1600-case benchmark for the `detect` calibration

For the algorithm + invariants see
[`docs/new/simulator_impl_plan.md`](new/simulator_impl_plan.md) and
[`docs/new/simulator_plan.md`](new/simulator_plan.md). The canonical
JSON Schema is in
[`docs/new/simulator_schema.json`](new/simulator_schema.json) and is
embedded into the binary for `synth-validate` / `synth-schema`.

## Subcommands

| Command | What it does |
|---|---|
| `synth <config.yaml> -o <prefix>` | Generate one synthetic array. Writes `<prefix>.fa`, `<prefix>.truth.tsv`, `<prefix>.periods.tsv`, and (with `--diagnostics`) `<prefix>.diagnostics.json` |
| `synth-batch --config-dir <dir> --out-dir <dir>` | Run `synth` over every `*.yaml` in a directory (parallel). `*.deferred.yaml` files are silently skipped |
| `synth-validate <config.yaml>` | Validate a YAML config against the canonical schema + MVP business rules. Exits non-zero on first error |
| `synth-schema` | Print the canonical JSON Schema to stdout |

## Output schemas

### `<prefix>.fa`

Standard FASTA, one record per array (the simulator generates one
array per config). Header is `>{array_id}`; sequence is wrapped at
`80` columns.

### `<prefix>.truth.tsv` (17 columns)

One row per array (a `synth-batch` run produces one row per config
across many files; a single `synth` run produces a single-row file).
Source: `src/synth/truth.rs::TRUTH_HEADER`.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | identifier from `global.array_id` in the YAML |
| 2 | `length_bp` | int | realised sequence length after all wobble + indel + event passes |
| 3 | `truth_class` | str | one of `simple_TR`, `HOR`, `mixed`, `random` |
| 4 | `base_width_bp` | int | base monomer length (logical) |
| 5 | `hor_k` | int / `NA` | HOR multiplicity; `NA` for `simple_TR` / `random` |
| 6 | `hor_length_bp` | int / `NA` | HOR unit length (`hor_k · base_width_bp`); `NA` when no HOR |
| 7 | `n_complete_copies` | int | realised number of complete HOR-unit (or monomer) copies after all events |
| 8 | `wobble_amplitude_bp` | float | realised sinusoidal envelope amplitude (`%.4f`); `0` when no wobble configured |
| 9 | `wobble_periodicity_bp` | float / `NA` | realised wobble period; `NA` when no wobble |
| 10 | `n_phase_shifts` | int | discrete phase-shift / insertion events |
| 11 | `phase_shift_positions` | list-int / `NA` | event positions in bp, comma-separated |
| 12 | `phase_shift_offsets` | list-int / `NA` | per-event offsets in bp (signed), comma-separated |
| 13 | `n_segments` | int | number of analysis segments (≥ 1; > 1 for `mixed`) |
| 14 | `reason` | str | free-text diagnostic — the simulator's intent label for this config |
| 15 | `structural_expression` | str | taxonomy §2 grammar string — e.g. `HOR{A1,A2,A3}×80` |
| 16 | `schema_version` | int | output schema version (currently `1`) |
| 17 | `events_json` | str | JSON-encoded list of all hybrid / inversion / dup / del events applied |

### `<prefix>.periods.tsv` (4 columns)

Period candidates the simulator considers worth testing — the **true**
base + HOR-unit periods plus a small set of distractors. This is what
`kitehor detect --periods` expects.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `array_id` | str | FASTA record identifier (matches the truth row) |
| 2 | `period_bp` | int | candidate period |
| 3 | `period_score` | float | candidate score in `[0, 1]` (truth periods get 1.0, distractors get lower scores) |
| 4 | `source` | str | origin tag — `true_base`, `true_hor_unit`, `distractor` |

### `<prefix>.diagnostics.json` (only with `--diagnostics`)

Structured dump intended for tooling. Fields include the raw YAML
config, per-block expansions, applied events, and intermediate
coordinate maps. See `src/synth/diagnostics.rs` for the exact shape.

## Worked example

```bash
# Generate one HOR_clean fixture
kitehor synth tests/synth_configs/T05_hor_clean.yaml -o /tmp/t05
head -1 /tmp/t05.truth.tsv     # column header
head -2 /tmp/t05.periods.tsv   # array_id  period_bp  period_score  source

# Validate before running
kitehor synth-validate tests/synth_configs/T05_hor_clean.yaml

# Generate the full CI corpus in parallel
kitehor synth-batch --config-dir tests/synth_configs --out-dir /tmp/synth_corpus
ls /tmp/synth_corpus/ | head   # T01_*.fa, T01_*.truth.tsv, T01_*.periods.tsv, ...

# Dump the canonical schema
kitehor synth-schema > simulator_schema.json
```

## YAML config

See [`docs/new/simulator_schema.json`](new/simulator_schema.json) for
the canonical schema and [`tests/synth_configs/`](../tests/synth_configs/)
for ready-to-modify reference configs spanning every taxonomy class.

Top-level keys:

- `global` — `array_id`, `seed`, `length_bp`, output knobs
- `templates` — named monomer / HOR_slot definitions
- `blocks` — ordered list of HOR / SIMPLE_TR / SHIFT / INSERTION blocks
- `events` — optional post-generation events (HYBRID, INVERSION, DUP, DEL)
- `noise` — final mutation + indel pass parameters

`synth-validate` enforces both the JSON Schema (structural) and a set
of MVP business invariants (e.g. block widths must sum to ≈ `length_bp`).

## Determinism

- The top-level RNG is seeded from `global.seed`.
- Sub-streams (wobble, noise, events) derive their own seed via
  FNV-1a hash of `parent_seed:stream_name` — matches the legacy
  simulator's convention.
- Re-running the same YAML with the same `--seed` override produces
  byte-identical output (FASTA, truth, periods, diagnostics).

## Source

- `src/synth/config.rs` — YAML loader + serde validation + MVP rules
- `src/synth/rng.rs` — FNV-1a sub-stream derivation
- `src/synth/templates.rs`, `blocks.rs`, `wobble.rs`, `events.rs`, `noise.rs` — pipeline stages
- `src/synth/truth.rs` — `truth.tsv` writer + `TRUTH_HEADER`
- `src/synth/periods.rs` — `periods.tsv` writer
- `src/synth/diagnostics.rs` — optional diagnostics dump
- `src/synth/simulator.schema.json` — embedded canonical schema (drift-tested)
