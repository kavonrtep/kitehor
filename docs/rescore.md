# `kitehor rescore` — pairwise tile-identity rescoring

`rescore` adds a nucleotide-level confidence signal to kite's peaks. For
each candidate period it samples adjacent tile pairs from the array and
computes their median pairwise identity. The output is kite's peaks TSV
with four appended columns: `identity_med`, `identity_iqr`,
`identity_p25`, `identity_n`.

The metric is **additive only**. Downstream stages (rule-classify,
analyze) still decide on kite's `score2_norm`; rescore is a diagnostic
column that downstream analysis can consult independently.

## Why rescore exists

Kite scores periodicity by k-mer set overlap in the neighbour-distance
histogram. That signal can fail to separate the **monomer period** from
the **HOR-unit period** in well-formed HORs — the k-mer pool is shared
between the two scales. Pairwise nucleotide identity does separate them:

- At P = monomer, adjacent tiles are *different* monomers within an HOR
  block and look ~80% identical.
- At P = HOR unit (k × monomer), adjacent tiles are consecutive *copies*
  of the same HOR and look ~95–99% identical.

The period with the higher `identity_med` is the more credible HOR unit
length. See `tests/rescore_smoke.rs` for the headline correctness assertion.

## CLI

```
kitehor rescore <FASTA>... --peaks <peaks.tsv> -o <prefix>
```

- `<FASTA>...` — one or more FASTA files containing the records named in
  `peaks.tsv`. Sequences are looked up by `case_id`; records missing from
  the FASTAs produce `NA` rows.
- `--peaks` — long-format peaks TSV emitted by `kitehor kite-periodicity`
  (header must contain `case_id`, `rank`, `period`).
- `-o <prefix>` — output is written to `<prefix>.peaks.tsv`. The stage
  refuses to overwrite any existing file at that path; pass `--force` to
  allow in-place rewriting (e.g. when `-o` resolves to the same file as
  `--peaks`).

### Flags

| flag | default | notes |
|---|---|---|
| `--samples K` | `200` | sampled pairs per (record, period); linear cost |
| `--slop` | `10` | bp of slack on the B-tile to absorb tile-boundary indels; must satisfy `slop ≤ period` |
| `--band` | `0` (auto) | indel-deviation tolerance for the banded kernel; auto resolves to `max(20, 2·slop)` |
| `--max-n-frac` | `0.05` | skip pairs whose combined N fraction exceeds this |
| `--max-retries` | `3` | extra draws per slot when an initial draw is N-rejected |
| `--min-period` | `20` | skip candidates below this; emit NA for those rows |
| `--max-period` | `5000` | skip candidates above this; `0` = unlimited |
| `--seed` | `42` | top-level RNG seed (deterministic per `(seed, case_id)`) |
| `--top-n` | `10` | only rescore the first N peaks per record; `0` = all |
| `--mismatch-cost` | `1` | per-cell cost of a mismatch (match is always 0) |
| `--gap-cost` | `1` | per-cell cost of an insertion or deletion (no affine gaps; ins == del) |
| `--min-array-bp` / `--max-n-fraction` | shared QC | inherits from `QcOpts` |
| `--threads` | `0` (auto) | rayon worker count |

### Scoring caveat

The defaults `--mismatch-cost 1 --gap-cost 1` give plain Levenshtein
edit distance, so `identity_med = 1 − edit_distance/|A|` is exactly the
matching fraction. With non-unit costs the returned value is a *weighted*
edit distance: `identity_med` stays in `[0, 1]` and ranks pairs the same
way, but no longer equals matches/|A|. Useful when you want to bias the
DP toward (or against) gaps — e.g. `--mismatch-cost 3 --gap-cost 1`
encourages gap-based alignments through divergent stretches.

### Runtime logging

With `-v` (info level), `rescore` emits three structured lines per run:

```
rescore: loaded 1600 record(s), 14466 peak row(s); 11842 to rescore (filters: min_period=20, max_period=5000, top_n=10)
rescore: K=200 slop=10 band=20 max_retries=3 seed=42 threads=16
rescore: 4231/11842 (35.7%) elapsed=120s rate=35/s eta=218s         ← every 10s
rescore: done in 187.4s — rescored 11815, filtered 2624, kernel-NA 27, identity_n median=200
```

`filtered` = rows blocked before the kernel by rank/period/missing-record;
`kernel-NA` = rows passed the filters but the kernel returned no usable
identity (short array or N-rejected all samples).

## Algorithm

For each (record, candidate period P):

1. Sample `K` anchor offsets uniformly from `[0, L − 2P − slop]` using a
   ChaCha20 PRNG seeded with FNV-1a of `(seed, case_id)`.
2. For each anchor `s`, form two windows:
   - **A** = `seq[s .. s + P]` (length P)
   - **B** = `seq[s + P − slop .. s + 2P + slop]` (length P + 2·slop)
3. Drop pairs whose combined N fraction exceeds `max_n_frac` and re-draw
   up to `max_retries` times.
4. Compute the **semi-global edit distance** of A against the best window
   inside B (A consumed end-to-end; B has free ends). Identity =
   `1 − edit_distance / P`.
5. Report `identity_med`, `identity_iqr`, `identity_p25` over the K
   identities, plus `identity_n` (effective sample count after rejection).

Sampling is **adjacent-tile only** (`d=1`). Multi-distance probing
(`d=2,3,…` for drift assessment) is a future flag, not v1.

### Edge cases

- Period below `min_period`, or `slop > period`, or `L < 2P + slop` ⇒ all
  four columns are `NA`, `identity_n = 0`.
- Record not found in any FASTA, or failed QC at load time ⇒ `NA` row.
- All sampled pairs N-rejected ⇒ `NA` row.

### N handling

The kernel treats `N` as matching nothing (including another `N`). The
sampler's skip-pair logic keeps the kernel from seeing N-heavy windows in
practice; the conservative match rule is just a safety net for the few Ns
that slip through.

## Output schema

`<prefix>.peaks.tsv` is the input file with four columns appended:

```
identity_med  identity_iqr  identity_p25  identity_n
```

- `identity_med`, `identity_iqr`, `identity_p25` — `%.4f` ∈ [0, 1].
- `identity_n` — effective sample count after N-rejection.
- All original cells are passed through **verbatim** (no float
  reformatting), so byte-equality is preserved on the unchanged columns.

## Performance

The kernel is banded semi-global DP at `O(P · band)` per pair (~50-100×
faster than plain DP on long-period candidates). Cost scales linearly in
`K`, in candidate period `P`, and in `band`. The default `max-period=5000`
cap and `top-n=10` together keep the long-period tail bounded.

Indicative wall times (1600-case `ground_truth_v2/` corpus, K=200, defaults,
16 cores):

| stage | time |
|---|---|
| `kite-periodicity` (input) | ~35 s |
| `rescore` (banded DP, Tier 1) | ~30–60 s |

For users who need the unbounded scan, set `--max-period 0`. The
`O(P · band)` cost makes this affordable at the price of ~10× wall time.
A banded Myers bitvector kernel (`O(P · band / 64)`) is the natural next
step if even the bounded run becomes the bottleneck.

## Calibration

Run against the 1600-case `ground_truth_v2/` corpus with default flags
(K=200, slop=10, band=20 auto, max_period=5000, top_n=10):

| Category | n | HOR-unit wins | mono identity | HOR identity | gap |
|---|---|---|---|---|---|
| hor_clean | 600 | 100.0% | 0.828 | 0.971 | +14.2 pp |
| hor_event_* (4 cats) | 200 | 100.0% | 0.821 | 0.961 | +13.9 pp |
| hor_insertion | 100 | 100.0% | 0.836 | 0.960 | +12.4 pp |
| hor_shift | 200 | 100.0% | 0.838 | 0.960 | +12.2 pp |
| hor_wobble | 100 | 100.0% | 0.839 | 0.960 | +12.1 pp |
| mixed | 100 | 67.0% | 0.676 | 0.779 | +10.4 pp |
| **TOTAL** | **1300** | **97.5%** | **0.819** | **0.951** | **+13.2 pp** |

A "win" means `identity_med` at the true HOR-unit period exceeded
`identity_med` at the true monomer period (lookup tolerance ±5% on
period). Period matches existed for every case; no NA rows in any
category.

The 33% loss rate on `mixed` reflects the underlying structural ambiguity
of interleaved HOR cases — when two distinct HORs share the array, a
period at one HOR's monomer can score higher local identity than the
other HOR's unit period. Banded edit distance correctly exposes this
ambiguity; the prior un-banded kernel masked it with over-permissive
substring matching.

## Worked example (smoke fixture)

```
$ kitehor kite-periodicity test_data/smoke/sequences.fasta -o /tmp/k.tsv
$ kitehor rescore test_data/smoke/sequences.fasta \
      --peaks /tmp/k.tsv.peaks.tsv -o /tmp/r --top-n 5

# case_id    rank period identity_med
hor_k3       1    300    0.9033   ← HOR unit
hor_k3       2    100    0.7400   ← monomer
hor_k5       1    750    0.9033   ← HOR unit
hor_k5       2    150    0.7633   ← monomer
tandem_pure  1    60     1.0000
```

kite's `score2_norm` correctly ranks the HOR-unit period first in both
HOR fixtures, but the identity gap to the monomer (0.90 vs ~0.75) is the
diagnostic signal a downstream stage can use to disambiguate "real HOR"
from "monomer-only tandem repeat that happens to expose harmonics in the
periodogram".
