# Rule-based pipeline (`kitehor analyze`)

End-to-end Rust port of the `tools/rule_proto/*.py` prototype. Five
subcommands cover the four stages plus the orchestrator:

```
FASTA
  │
  ▼
kite-periodicity       k-mer pair-distance periodogram (existing)
  │
  ▼
rule-classify          HOR / simple_tr / unresolved verdict
  │
  ├──────────────────────────┐
  ▼                          ▼
tandem-validate              ssr-scan
(spatial localization        (short motifs)
 of any sub-host period)
  │                          │
  └──────────────────────────┘
                 │
                 ▼
            summary-merge       outer-join + 7-rule combined_class
                 │
                 ▼
          <prefix>.summary.tsv
```

The seven `combined_class` values:

| Class | Fires when |
|---|---|
| `pure_ssr` | `ssr_flag = yes` AND `ssr_dominant_motif_coverage_pct ≥ 80` |
| `tr_with_subrepeat` | `tv_decision = localized_subrepeat` |
| `hor_with_ssr` | `hor_verdict = hor` AND `ssr_flag = yes` |
| `hor` | `hor_verdict = hor` |
| `tr_with_ssr` | `hor_verdict = simple_tr` AND `ssr_flag = yes` |
| `tr` | `hor_verdict = simple_tr` |
| `unresolved` | none of the above |

`tr_with_subrepeat` covers what the prior 8-class cascade split between
`tr_with_nested_tr` (array-scale heterogeneity, found by the old
`subrepeat-scan`) and `tr_with_subrepeat` (within-tile heterogeneity,
found by the old `hor-validate`). The new `tandem-validate` stage
tests both scales with one density + spatial / phase-contrast check; see
[`docs/new/tandem_validate_spec.md`](new/tandem_validate_spec.md) for
the algorithm and [`docs/new/tandem_validate_port_plan.md`](new/tandem_validate_port_plan.md)
for the v0.10 retirement notes.

## Quick start

```bash
# End-to-end
kitehor analyze <fasta> -o <prefix>

# Per stage (debugging / partial rerun)
kitehor kite-periodicity <fasta> -o <prefix>.kite.tsv --classify \
                                 --out-peaks <prefix>.kite.peaks.tsv
kitehor rule-classify   <prefix>.kite.peaks.tsv -o <prefix>
kitehor tandem-validate <fasta> --verdicts <prefix>.verdicts.tsv \
                                --peaks <prefix>.kite.peaks.tsv -o <prefix>
kitehor ssr-scan        <fasta> --kite-peaks <prefix>.kite.peaks.tsv -o <prefix>
kitehor summary-merge   --verdicts <prefix>.verdicts.tsv \
                        --tandem-validate <prefix>.tandem_validate.tsv \
                        --ssr <prefix>.ssr.tsv \
                        -o <prefix>
```

## Outputs (TSV-per-stage contract)

`analyze` always emits all seven per-stage TSVs under `<prefix>.*`:

| Stage | File(s) | Column count |
|---|---|---:|
| kite-periodicity | `.kite.tsv`, `.kite.peaks.tsv` | 9 / 9 |
| rule-classify | `.verdicts.tsv` | 10 |
| tandem-validate | `.tandem_validate.tsv` | 16 |
| ssr-scan | `.ssr.tsv`, `.ssr.regions.tsv` | 17 / 8 |
| summary-merge | `.summary.tsv` | 32 |

## Optional periodogram bundle (`--periodogram`)

Both `kitehor analyze` and `kitehor kite-periodicity` accept a
`--periodogram <PATH>` flag that emits a FASTA-like bundle of the
per-record neighbour-distance histogram. The data mirrors what
TideCluster keeps in its in-memory `profile_list` (see
`tarean/kite.R`); kitehor invents a text format because TideCluster only
persists the same data as a binary `peaks_list.RDS`.

Format — two records per input sequence:

```text
>case_id|H length=<N> kmer=<K>
<H[1]> <H[2]> ... <H[N]>
>case_id|bg length=<N> kmer=<K>
<bg[1]> <bg[2]> ... <bg[N]>
```

- `|H` is the raw neighbour-distance histogram (integer counts, formatted
  without a decimal point).
- `|bg` is the smoothed, composition-matched random background envelope
  (floats with 6 fractional digits).
- Each vector covers period `d = 1..N` where `N = length_bp` of the
  record. Index 0 is unused upstream and not emitted.
- Header tokens after the record id (`length=`, `kmer=`) are
  whitespace-separated `key=value` pairs.
- Records whose array failed kite analysis (e.g. too short) are skipped.

Quick load + plot in Python:

```python
import numpy as np, matplotlib.pyplot as plt
def iter_records(path):
    with open(path) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
            else:
                yield header, np.fromstring(line, sep=" ")

curves = {}
for h, v in iter_records("smoke.periodogram"):
    case, channel = h.split()[0].split("|")
    curves.setdefault(case, {})[channel] = v
case = next(iter(curves))
plt.plot(curves[case]["H"], label="H")
plt.plot(curves[case]["bg"], label="bg")
plt.xlabel("period (bp)"); plt.legend(); plt.title(case); plt.show()
```

## Algorithm details

### `rule-classify`

Clusters kite peaks by relative period gap (single-linkage on
`(p_cur − p_prev) / p_cur ≤ tol`), then walks a first-match-wins
decision tree:

1. **Case A** — top cluster is `k × shorter cluster`. Picks the
   smallest-period qualifying divisor. HOR(founder, k, tile=top).
2. **Case B walk** — top is the founder; k = 2..k_max looks for a
   cluster at `k × top.period` that produces a non-monotonic bump
   above the running max (k ≥ 3) or exceeds `non_mono_ratio × top`
   (k = 2). Requires the harmonic confirmation
   `score(2k·p) ≥ score((k+1)·p)`.
3. **monotonic_multiples** — has at least one larger-multiple cluster
   but no qualifying bump → `simple_tr`.
4. **lone_significant_cluster** — exactly one cluster passes
   `lone_significant_frac × top` → `simple_tr`.
5. **unresolved** — none of the above.

All hardcoded constants exposed as CLI flags. See `--help`.

### `tandem-validate`

Port of `tools/rule_proto/tandem_validate.py` (spec v5). Replaced the
prior `subrepeat-scan` + `hor-validate` stages in v0.10. The
algorithm contract is
[`docs/new/tandem_validate_spec.md`](new/tandem_validate_spec.md);
the summary here mirrors its §3 (algorithm) and §4 (decision tree).

For each record, given the verdict's host period (HOR → `tile`,
`simple_tr` → `founder`, `unresolved` → kite rank-1):

1. **Skip HOR k=2** — only 2 phase bins fit in `host`; the test is
   geometrically degenerate. Returns `skip_k2`; the cascade falls
   through to `hor`.
2. **Pick candidates**: any kite peak with `cand_min_period ≤ p <
   host / 3` and `score2_norm ≥ max(cand_score_floor,
   cand_rel_score_floor × record_max_score)` (default relative floor
   0.03). For HOR k ≥ 3, the founder is added as a `kind=Founder`
   candidate even when it sits at the boundary (1% slack); peaks
   within `founder_tol` of any `m·founder` rung (`m = 1..k-1`) are
   excluded. Non-founder candidates are `kind=Other`, capped to top
   `cand_top_n` (default 5) by score.
3. **Plan windows**: `w = max(host / 3, 3·max_candidate,
   min_window_bp)`, then `w = min(w, host)` (hard cap at host
   preserves within-host phase resolution). Step is `w / 4`.
4. **Per-window kite** (in-process, rayon-parallel) on each window.
5. **Per-window presence**, gated by candidate kind:
   * `Founder` (loose): `sum(scores ±tol of cand) ≥
     presence_rel_floor × top_score` (default 0.2). Picks up the
     founder even when the tile is the in-window top, as is the
     case for any clean HOR.
   * `Other` (strict): window's top period IS the candidate AND
     `top_score ≥ window_score_floor` (default 0.3). The score
     floor distinguishes legit heterogeneity (some weak-top
     windows) from uniform tandem (all-strong-top windows).
6. **Metrics** per candidate: `density = n_present / n_total`,
   `spatial_contrast = max−min over 10 array-position bins`,
   `phase_contrast = max−min over 10 (mid mod host) bins`
   (computed iff `window_bp < host × 0.95`, i.e., when a window can
   fit at multiple phase positions of one host cycle).
7. **Per-candidate decision** (after `n_present < min_present_windows
   → no_signal` short-circuit):
   * `localized` iff `density ≤ density_dup_max` OR any
     `contrast ≥ contrast_dup_min`.
   * `uniform` iff `density ≥ density_hor_min` AND both contrasts
     `≤ contrast_hor_max`.
   * else `ambiguous`.
8. **Per-record decision**: best candidate by
   `(rank, −(spatial + phase))` where
   `rank = localized < ambiguous < uniform < no_signal`. Output
   `decision_hint`:
   * `localized` → `localized_subrepeat` → cascade fires
     `tr_with_subrepeat`
   * `uniform` → `confirms_host` → cascade falls through to verdict
   * `ambiguous` / `no_signal` → falls through

Every threshold is exposed as a CLI flag; defaults match the Python
prototype.

### `ssr-scan`

1. Raw scan — for each motif length L = 1..14, find non-overlapping
   greedy runs of an L-base motif repeating ≥ `min_reps` times. The
   prototype uses regex `(([gatc]{L})\2{min-1,})`; the Rust port
   hand-rolls the scanner (Rust's `regex` crate doesn't support
   backreferences). Output is `(motif_lower, start_1based,
   end_0based_exclusive, repeats, normalized_motif)`.
2. `normalize_motif` = lex-min over all rotations of `motif.upper()`
   AND of its reverse complement.
3. `get_unique_motifs` drops any motif that is exactly `k × shorter`.
4. **Consensus path** (when `--kite-peaks` is supplied): from the kite
   top peak's period `P`, extract up to `consensus_max_monomers`
   canonical-distinct P-mers from the sequence; validate each by
   building a dimer (`P × consensus_dimer_copies` or
   `consensus_dimer_min_bp`) and running the raw scan on it. The
   number of unique validated canonicals selects:
   - 0 → `raw_fallback`
   - 1 → `consensus_single` (use the dimer's summary as authoritative)
   - ≥ 2 → `consensus_multi` (per-motif coverage from the raw scan)

### `summary-merge`

Outer-join on `record_id` (sorted lexicographically — matches
pandas's outer-merge behavior in the user's env). Defaults fill
missing values; the 7-rule decision tree determines
`combined_class`. Float columns use `%.4g`.

The summary schema dropped the old subrepeat- and hor-validate-prefixed
columns and gained `tv_*` columns sourced from the new
`tandem_validate.tsv`. `length_bp` is no longer included; join on
`record_id` against `<prefix>.kite.tsv` if you need it.

## Float formatting policy

Three policies coexist, matching the prototype's per-script
defaults so the Rust output is byte-equivalent with the Python
prototype:

- `rule-classify` → `%.6g`
- `summary-merge` → `%.4g`
- `ssr-scan` → pandas default (≈ shortest roundtrip + `.0` for
  integer-valued floats)
- `tandem-validate` → `%.6g`

The shared helper `rule_classify::io::fmt_g(precision, x)`
implements Python's `%g` semantics (significand-precision, scientific
fallback for very small/large magnitudes, trailing-zero stripping).

## Byte-equivalence with the Python prototype

The `rule-classify`, `ssr-scan`, and `summary-merge` stages are
validated byte-identical with the prototype on the smoke fixture
(`test_data/smoke/sequences.fasta`) and on the hand-curated
fixtures in `tools/rule_proto/fixtures/` (`tests/rule_classify_fixtures.rs`).

`tandem-validate` is validated against the Python prototype via
`tests/tandem_validate_python_parity.rs` (`#[ignore]`-flagged,
needs `python3` + `pandas` + the prototype at
`tools/rule_proto/tandem_validate.py`). The test asserts the
`decision_hint` column matches for every record in a 6-case
synthetic fixture covering the `skip_k2`, `confirms_host`,
`no_candidates`, and `localized_subrepeat` paths.

## Replaced modules

`P1` retired `src/rule.rs` (the older 4-condition rule). `P6` removed
the legacy ML pipeline (`classifier.rs`, `classify.rs`, `features.rs`,
`hor_call.rs`, `coverage.rs`, `models/`, `examples/validate_rf.rs`,
`config/classifier.toml`) along with their CLI flags
(`--use-ml-classifier`, `--no-hor-call`, `--hor-*`, `--rule-qmax`,
`--rule-top-n`, `--coverage`, `--classifier-config`, `--hor-model`,
`--k-model`, `--no-homology`).

`v0.10` retired `src/subrepeat/` and `src/hor_validate/` along with
their `subrepeat-scan` / `hor-validate` subcommands, the
`summary-merge --subrepeat` / `--within-tile` flags, and the
`tr_with_nested_tr` combined class — all replaced by the unified
`tandem-validate` stage and the 7-class cascade. See
[`docs/new/tandem_validate_port_plan.md`](new/tandem_validate_port_plan.md)
for the rollout sequence.

## Performance

Targets on test_590 (2779 records, ~280 MB FASTA), single machine,
default rayon parallelism. v0.10 figures roll up the prior
`subrepeat-scan` + `hor-validate` budgets into `tandem-validate`:

| Stage | Python wall | Rust target |
|---|---:|---:|
| kite-periodicity | 0:15 | 0:15 |
| rule-classify | 0:10 | < 0:01 |
| tandem-validate | 3:15 | < 0:35 |
| ssr-scan | 7:00 | < 0:20 |
| summary-merge | 0:01 | < 0:01 |
| **analyze (end-to-end)** | **~10:40** | **< 1:30** |
