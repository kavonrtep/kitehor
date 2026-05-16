# ground_truth_v2/ — v2 synth corpus (1,600 cases)

Structured corpus for the `kitehor synth` (v2 YAML-driven) simulator,
the analog of [`../ground_truth/`](../ground_truth/) for the legacy
`simulate-grid` pipeline.

This directory tracks **the spec** (generator script, per-case config
YAMLs, manifest). Generated FASTA + truth + periods bundles live
**outside** the repo per kitehor's data policy ("no large simulator
output in the repo") — regenerate them on demand.

## What's tracked

```
ground_truth_v2/
  README.md             this file
  generate_configs.py   deterministic Python generator (MD5-derived seeds)
  manifest.tsv          per-case parameter table (1600 rows)
  run_batch.sh          one-shot runner: synth-batch per category
  configs/
    01_simple_tr/         200 configs
    02_hor_clean/         600
    03_hor_wobble/        100
    04_hor_shift/         200
    05_hor_insertion/     100
    06_hor_event/         200
    07_mixed/             100
    08_random/             50
    09_gc_bias/            50
                       -----
                        1600
```

`configs/` and `manifest.tsv` are byte-reproducible from
`generate_configs.py`. Re-running the generator overwrites them
identically.

## What's NOT tracked

`out/` (1,600 × {`.fa`, `.truth.tsv`, `.periods.tsv`} = ~600 MB of
generated data) is gitignored. Put it anywhere; default location used
by `run_batch.sh` is `./out/` next to this README.

## Categories

| # | category | n | varies on |
|--:|---|--:|---|
| 01 | simple_tr     | 200 | monomer_len × n_copies × mutation × indel |
| 02 | hor_clean     | 600 | monomer_len × k × n_copies × divergence × mutation |
| 03 | hor_wobble    | 100 | amplitude × period_rows × (sin / random_walk) |
| 04 | hor_shift     | 200 | shift_offset × n_copies × monomer_len × k |
| 05 | hor_insertion | 100 | kind × length × position |
| 06 | hor_event     | 200 | HYBRID 50 + INVERSION 50 + DUP 50 + DEL 50 |
| 07 | mixed         | 100 | (L_a, k_a) × (L_b, k_b) coexisting HORs |
| 08 | random        |  50 | non-tandem INSERTIONs only (negative control) |
| 09 | gc_bias       |  50 | simple TRs at GC ∈ {0.1..0.8} |

`manifest.tsv` carries the full per-case parameter set so the eval
harness can join by `case_id`.

## Regenerating

```bash
# 1. (Re)write configs/ + manifest.tsv from the generator.
python3 ground_truth_v2/generate_configs.py

# 2. Run kitehor synth-batch over the corpus.
./ground_truth_v2/run_batch.sh                     # → ground_truth_v2/out/
# or override the output root:
SYNTH_OUT=/path/to/big/disk ./ground_truth_v2/run_batch.sh
```

The full batch completes in ~2 s on 16 CPUs (~860 cases/s,
~334 Mbp/s), producing ~600 MB of FASTA.
