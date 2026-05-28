# `kitehor simulate` / `simulate-grid` — legacy params-TSV simulator

The original simulator, retained for compatibility with the legacy
`ground_truth/` benchmark (1600-case corpus spec). Generates HOR /
simple tandem repeats with per-monomer substitution + indel noise,
gene-conversion-style block / monomer events, and optional sub-motif
tiling inside the founder.

For new work prefer the v2 simulator (`kitehor synth*`) which covers
the full taxonomy including wobble, phase shifts, hybrid monomers, and
inversions. See [`docs/synth.md`](synth.md).

## Subcommands

| Command | What it does |
|---|---|
| `simulate -o <fasta>` | Generate one synthetic array from CLI args. Writes `<fasta>` plus `<fasta>.truth.tsv` next to it |
| `simulate-grid --params <params.tsv> --outdir <dir>` | Generate a grid of arrays from a parameters TSV. Writes `sequences.fasta`, `truth.tsv`, `monomers.tsv`, `events.tsv`, and `alternatives.tsv` to `<dir>` |

## `simulate` — single array

```bash
kitehor simulate \
    --monomer-size 171 --multiplicity 12 --copies 100 \
    --sub-rate-intra 0.05 --sub-rate-inter 0.03 \
    --case-id sim_0001 --seed 0 \
    -o sim_0001.fa
```

Outputs:

- `<out>` — single-record FASTA (`>case_id\nACGT...`)
- `<out>.truth.tsv` — single-row truth file (same schema as
  `simulate-grid::truth.tsv` below)

### `simulate` flags (defaults)

| flag | default | meaning |
|---|---|---|
| `--monomer-size <N>` | `171` | base monomer length in bp |
| `--multiplicity <K>` | `12` | HOR multiplicity; `1` produces a plain tandem repeat |
| `--copies <N>` | `100` | number of HOR (or monomer) copies in the array |
| `--sub-rate-intra <F>` | `0.05` | per-base substitution rate within each monomer copy |
| `--sub-rate-inter <F>` | `0.03` | per-base substitution rate between founder positions |
| `--submono-k <K>` | `1` | sub-motif tiling factor inside the founder (`1` = none) |
| `--seed <N>` | `0` | RNG seed |
| `--case-id <id>` | `sim_0000` | FASTA record id |

## `simulate-grid` — params-driven grid

```bash
kitehor simulate-grid \
    --params ground_truth/params.tsv \
    --outdir ground_truth/ \
    --seed 42 \
    --threads 0
```

Reads a params TSV (one row per case), regenerates the whole bundle in
~7 s for 1600 cases on default parallelism. The per-case seed is
derived deterministically from `master_seed:case_id` via FNV-1a so
re-running with the same `--seed` is byte-identical.

### Input — `params.tsv`

| column | meaning |
|---|---|
| `case_id` | unique identifier |
| `monomer_len` | base monomer length in bp |
| `hor_order` | HOR multiplicity (`1` = simple TR) |
| `n_blocks` | number of HOR (or monomer) blocks |
| `sub_rate_intra` | within-copy substitution rate |
| `sub_rate_inter` | between-founder substitution rate |
| `indel_rate_intra` | within-copy indel rate |
| `indel_rate_inter` | between-founder indel rate |
| `block_conversions` | gene-conversion-style block-scale events |
| `monomer_conversions` | gene-conversion-style monomer-scale events |
| `submono_k` | sub-motif tiling factor (`1` = none) |
| `seed` | optional per-case seed (blank → derive from `--seed`) |

## Output schemas

### `sequences.fasta`

Standard FASTA, one record per `case_id`.

### `truth.tsv` (18 columns)

One row per `case_id`. The first 12 columns echo the params; columns
13–18 are realised measurements computed during simulation.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | identifier |
| 2 | `monomer_len` | int | base monomer length used (= params) |
| 3 | `hor_order` | int | HOR multiplicity used (= params) |
| 4 | `n_blocks` | int | number of HOR / monomer blocks used (= params) |
| 5 | `sub_rate_intra` | float | within-copy substitution rate (= params) |
| 6 | `sub_rate_inter` | float | between-founder substitution rate (= params) |
| 7 | `indel_rate_intra` | float | within-copy indel rate (= params) |
| 8 | `indel_rate_inter` | float | between-founder indel rate (= params) |
| 9 | `block_conversions` | int | block-scale conversion events applied (= params) |
| 10 | `monomer_conversions` | int | monomer-scale conversion events applied (= params) |
| 11 | `submono_k` | int | sub-motif tiling factor (= params) |
| 12 | `seed` | int | effective seed used |
| 13 | `array_length` | int | realised sequence length in bp after all events |
| 14 | `n_monomers` | int | realised number of monomer copies |
| 15 | `mean_intra_block_id` | float | mean identity across copies of the same founder within the same block |
| 16 | `mean_homologous_id` | float | mean identity across copies of the same founder across all blocks |
| 17 | `mean_cross_position_id` | float | mean identity between founders at different HOR positions |
| 18 | `hor_signal` | float | (`mean_homologous_id − mean_cross_position_id`); HOR strength proxy in `[−1, 1]` |

### `monomers.tsv` (7 columns)

One row per individual monomer copy.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | identifier |
| 2 | `monomer_idx` | int | 0-based ordinal within the array |
| 3 | `block_idx` | int | 0-based HOR-block ordinal (always `0` when `hor_order = 1`) |
| 4 | `founder_idx` | int | 0-based founder slot (`0..hor_order-1`); identifies which founder this copy descends from |
| 5 | `start` | int | 0-based inclusive start bp |
| 6 | `end` | int | 0-based exclusive end bp |
| 7 | `length` | int | `end − start` (≠ `monomer_len` when indels landed in this copy) |

### `events.tsv` (5 columns)

One row per gene-conversion-style event applied.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | identifier |
| 2 | `event_order` | int | 1-based application order |
| 3 | `scope` | str | `block` or `monomer` |
| 4 | `source_idx` | int | source ordinal (block_idx or monomer_idx) |
| 5 | `target_idx` | int | target ordinal — receives the source's content |

### `alternatives.tsv` (6 columns)

One row per alternative-valid (`tile`, `multiplicity`, `founder`)
hierarchy for sub-monomer / sub-HOR cases. The eval harness scores a
prediction against the closest valid hierarchy in this table rather
than only against the primary truth label (see `MISTAKE_TRIAGE.md`
§12). For most cases this file has zero or one row per `case_id`;
sub-motif tilings (`submono_k > 1`) typically produce multiple rows.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | identifier |
| 2 | `rank` | int | 0-based preference order (`0` = the canonical truth hierarchy) |
| 3 | `tile` | int | tile period (bp) at this hierarchy |
| 4 | `multiplicity` | int | HOR multiplicity at this hierarchy |
| 5 | `founder` | int | founder period (bp) at this hierarchy |
| 6 | `kind` | str | which structural level this hierarchy represents — e.g. `primary`, `submono_k_split`, `monomer_only` |

## Determinism

- Per-case seed = FNV-1a hash of `"<master_seed>:<case_id>"` when the
  row's `seed` cell is blank; otherwise the row's explicit value.
- Re-running the same params TSV with the same `--seed` is
  byte-identical across all four output files.

## Source

- `src/simulate.rs` — single-array path
- `src/simulate_grid.rs` — grid path + truth / monomers / events /
  alternatives writers (`TRUTH_HEADER`, `MONOMERS_HEADER`,
  `EVENTS_HEADER`, `ALTERNATIVES_HEADER`)
- `ground_truth/params.tsv` — reference params spec (1600 cases)
- `ground_truth/simulate_hor.py` — Python reference simulator kept
  for cross-checking; the Rust path is the primary tool
