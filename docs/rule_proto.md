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
            summary-merge       outer-join + 9-rule combined_class
                 │
                 ▼
          <prefix>.summary.tsv
```

The nine `combined_class` values (v0.11+). SSR-driven decisions all
fire against `ssr_raw_total_coverage_pct` (the array-scale total from
the raw scanner), against two thresholds:
`pure_ssr_pct_threshold = 80` and `ssr_has_pct_threshold = 30`.

| # | Class | Fires when |
|---|---|---|
| 1 | `pure_ssr` | `ssr_raw_total_coverage_pct ≥ 80` |
| 2 | `tr_with_subrepeat_with_ssr` | `tv_decision = localized_subrepeat` AND `ssr_raw_total_coverage_pct ≥ 30` |
| 3 | `tr_with_subrepeat` | `tv_decision = localized_subrepeat` |
| 4 | `hor_with_ssr` | `hor_verdict = hor` AND `ssr_raw_total_coverage_pct ≥ 30` |
| 5 | `hor` | `hor_verdict = hor` |
| 6 | `tr_with_ssr` | `hor_verdict = simple_tr` AND `ssr_raw_total_coverage_pct ≥ 30` |
| 7 | `tr` | `hor_verdict = simple_tr` |
| 8 | `unresolved_with_ssr` | `ssr_raw_total_coverage_pct ≥ 30` (i.e., `unresolved` verdict + SSR) |
| 9 | `unresolved` | none of the above |

`tr_with_subrepeat` covers what the prior 8-class cascade split between
`tr_with_nested_tr` (array-scale heterogeneity, found by the old
`subrepeat-scan`) and `tr_with_subrepeat` (within-tile heterogeneity,
found by the old `hor-validate`). The new `tandem-validate` stage
tests both scales with one density + spatial / phase-contrast check; see
[`docs/new/tandem_validate_spec.md`](new/tandem_validate_spec.md) for
the algorithm and [`docs/new/tandem_validate_port_plan.md`](new/tandem_validate_port_plan.md)
for the v0.10 retirement notes.

v0.11 closed an SSR cascade bug: under v0.10 the
`consensus_single` path in `ssr-scan` set
`dominant_motif_coverage_pct ≈ 100` (the candidate monomer's *self*
coverage on a synthetic dimer, not the array), and the cascade
threshold check on that field fired `pure_ssr` for arrays as low as
~4% SSR. The fix routes the cascade through
`ssr_raw_total_coverage_pct` (always array-scale) and recomputes the
`ssr_flag` column from the same raw total so it stops contradicting
the call. The consensus path's `dominant_motif` / `_length` / etc.
remain in `summary.tsv` as informational labels for "which short
motif dominates inside the kite top period".

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
| summary-merge | `.summary.tsv` | 33 |

Per-column reference for each file is in **[Output schemas](#output-schemas)** below.

## Output schemas

Every cell is tab-separated. `NA` is used for missing / not-applicable
values in `rule-classify`, `tandem-validate`, `summary-merge`, and
`ssr-scan`. Empty file (zero-length, just a trailing newline) is emitted
for `<prefix>.ssr.regions.tsv` when no SSR motif clears the
`--min-reps` floor — matches `pd.DataFrame([]).to_csv(...)`.

### `<prefix>.kite.tsv` — kite-periodicity summary (9 columns)

One row per FASTA record. `monomer_size_*` / `score_*` carry the top-3
kite peaks; `NA` whenever fewer than that survived kite's filters.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | FASTA record identifier |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `n_peaks_kept` | int | number of peaks kite retained after its own filtering |
| 4 | `monomer_size` | int / `NA` | rank-1 peak period (bp) |
| 5 | `score` | float / `NA` | rank-1 peak raw score (kite `%.10f`) |
| 6 | `monomer_size_2` | int / `NA` | rank-2 peak period |
| 7 | `score_2` | float / `NA` | rank-2 peak raw score |
| 8 | `monomer_size_3` | int / `NA` | rank-3 peak period |
| 9 | `score_3` | float / `NA` | rank-3 peak raw score |

With `--classify` (single-stage shortcut around `rule-classify`),
columns 10–14 (`verdict`, `founder`, `multiplicity`, `tile`, `share`)
are appended; semantics match `.verdicts.tsv` below.

### `<prefix>.kite.peaks.tsv` — kite peaks long-format (9 columns)

One row per kept peak per record (multi-row per record). Consumed by
every downstream stage that needs the period candidates.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | FASTA record identifier |
| 2 | `array_length` | int | sequence length in bp |
| 3 | `rank` | int | rank within the record (1 = highest `score2_norm`) |
| 4 | `period` | int | candidate period (bp) |
| 5 | `peak_height` | float | raw histogram height at this period |
| 6 | `score` | float | kite's raw periodogram score (peak vs local background) |
| 7 | `score2` | float | `score · log2(period)` — kite's log-weighted variant |
| 8 | `score2_norm` | float | `score2` divided by the per-record sum (sums to 1 across the row's peaks) |
| 9 | `background` | float | local background estimate at this period |

### `<prefix>.verdicts.tsv` — rule-classify output (10 columns)

One row per record. The classifier walks a first-match-wins decision
tree (see [Algorithm details § rule-classify](#rule-classify)) and
emits a verdict with founder + tile + multiplicity if applicable.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `case_id` | str | FASTA record identifier |
| 2 | `verdict` | str | `hor` / `simple_tr` / `unresolved` |
| 3 | `founder` | int / `NA` | founder (monomer) period in bp; `NA` when `verdict != hor` |
| 4 | `multiplicity` | int / `NA` | HOR multiplicity `k`; `NA` when `verdict != hor` |
| 5 | `tile` | int / `NA` | HOR tile period (`k · founder`); `NA` when `verdict != hor` |
| 6 | `founder_score` | float / `NA` | clustered score at the founder period (`%.6g`) |
| 7 | `tile_score` | float / `NA` | clustered score at the tile period (`%.6g`) |
| 8 | `confidence` | float / `NA` | per-call confidence; computation depends on verdict (see [`docs/hor_confidence_score_calculation.md`](hor_confidence_score_calculation.md)) |
| 9 | `n_clusters` | int | number of period clusters after single-linkage at `--rule-cluster-tol` |
| 10 | `reason` | str | decision-path tag — e.g. `case_a_k=3`, `case_b_walk_k=2`, `monotonic_multiples`, `lone_significant_cluster`, `min_tile_copies` (gated out), `no_peaks` |

### `<prefix>.tandem_validate.tsv` — tandem-validate output (16 columns)

One row per record. Spatial-localization test on sub-host periods (the
unified replacement for the prior `subrepeat-scan` + `hor-validate`
stages). Full algorithm: [`docs/new/tandem_validate_spec.md`](new/tandem_validate_spec.md).

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `verdict` | str | the input verdict (echoed from `verdicts.tsv` for self-contained downstream consumption) |
| 3 | `host_period` | int / `NA` | the host period the test ran against — HOR `tile`, `simple_tr` `founder`, or kite rank-1 for `unresolved` |
| 4 | `multiplicity` | int / `NA` | HOR `k` when applicable; `NA` otherwise |
| 5 | `window_bp` | int / `NA` | sliding-window width (`max(host/3, 3·max_candidate, min_window_bp)`, capped at `host`) |
| 6 | `n_candidates` | int | number of sub-host periods tested |
| 7 | `candidates` | list-str | `period:kind` entries, `;`-separated (`kind` ∈ `Founder` / `Other`) |
| 8 | `best_candidate_period` | int / `NA` | the candidate that drove the row's decision |
| 9 | `best_candidate_kind` | str / `NA` | `Founder` or `Other` |
| 10 | `density` | float / `NA` | fraction of windows where the best candidate was present (`%.6g`) |
| 11 | `spatial_contrast` | float / `NA` | max − min of presence across 10 array-position bins |
| 12 | `phase_contrast` | float / `NA` | max − min of presence across 10 `(mid mod host)` bins; `NA` when `window_bp ≥ host · 0.95` |
| 13 | `n_windows_total` | int | total sliding windows planned |
| 14 | `n_windows_present` | int | windows where the best candidate was present |
| 15 | `decision_hint` | str | `localized_subrepeat` / `confirms_host` / `ambiguous` / `no_signal` / `no_candidates` / `skip_k2` |
| 16 | `reason` | str | free-text diagnostic (path through the decision tree) |

### `<prefix>.ssr.tsv` — ssr-scan summary (17 columns)

One row per record. Mixes the **raw** array-scale SSR scan with an
optional **consensus** path that validates a kite top-period monomer
against the same scanner.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `length_bp` | int | sequence length in bp |
| 3 | `ssr_flag` | str | `yes` / `no` — derived from `raw_total_coverage_pct ≥ ssr_has_pct_threshold` (v0.11+) |
| 4 | `dominant_motif` | str / `NA` | dominant canonical motif (consensus path when present, else raw top) |
| 5 | `dominant_motif_length` | int / `NA` | bp length of `dominant_motif` |
| 6 | `dominant_motif_repeats` | int | total repeat units of `dominant_motif` |
| 7 | `dominant_motif_coverage_pct` | float | coverage of `dominant_motif` (consensus path: dimer self-coverage; raw path: array coverage) — **informational, not used by the cascade** |
| 8 | `total_ssr_coverage_pct` | float | total SSR coverage on whatever path was chosen (consensus or raw); equals `raw_total_coverage_pct` on the raw path |
| 9 | `top_motifs` | list-str | top-3 motifs from the chosen path: `motif:pct;motif:pct;motif:pct` |
| 10 | `ssr_method` | str | which path produced cols 4–9: `raw_fallback` / `consensus_single` / `consensus_multi` |
| 11 | `consensus_period_bp` | int / `NA` | kite top period that drove the consensus path; `NA` when not applicable |
| 12 | `consensus_monomer` | str / `NA` | canonical monomer used for consensus dimer validation |
| 13 | `ssr_raw_dominant_motif` | str / `NA` | raw-scan dominant motif (always array-scale, regardless of consensus path) |
| 14 | `ssr_raw_dominant_motif_coverage_pct` | float | raw-scan dominant motif's array coverage |
| 15 | `ssr_raw_total_coverage_pct` | float | **canonical SSR coverage signal** — sum of all raw-scan motifs' array coverage (cap 100 %). Drives every cascade SSR decision in v0.11+ |
| 16 | `ssr_raw_n_regions` | int | total raw-scan SSR regions in the array |
| 17 | `ssr_raw_top_motifs` | list-str | top-3 raw-scan motifs: `motif:pct;motif:pct;motif:pct` |

### `<prefix>.ssr.regions.tsv` — ssr-scan per-region (8 columns)

One row per individual SSR region detected by the raw scan (zero rows
⇒ empty file, no header). Useful for diagnosing localised SSR
intrusions inside otherwise non-SSR arrays.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `ssr_number` | int | 1-based ordinal within the record |
| 3 | `motif_length` | int | length of the underlying motif in bp |
| 4 | `motif_sequence` | str | the raw motif as found (lowercase, in-place orientation) |
| 5 | `repeats` | int | number of consecutive tandem copies of the motif |
| 6 | `start` | int | 1-based inclusive start position (matches the prototype's `re.finditer().start() + 1`) |
| 7 | `end` | int | 0-based exclusive end position (i.e. `start + repeats · motif_length − 1` is the last-base 1-based coordinate; `end` here matches the prototype's `re.finditer().end()` semantics) |
| 8 | `normalized_motif` | str | lex-min over all rotations of `motif.upper()` and its reverse complement |

### `<prefix>.summary.tsv` — summary-merge output (33 columns)

One row per record. Outer-joins the four upstream stages and adds the
`combined_class` decision. Floats use `%.4g`.

| # | column | type | meaning |
|---|---|---|---|
| 1 | `record_id` | str | FASTA record identifier |
| 2 | `hor_verdict` | str | `verdicts.tsv::verdict` (`hor` / `simple_tr` / `unresolved`); `unresolved` when join missing |
| 3 | `hor_founder` | int / `NA` | `verdicts.tsv::founder` |
| 4 | `hor_multiplicity` | int / `NA` | `verdicts.tsv::multiplicity` |
| 5 | `hor_tile` | int / `NA` | `verdicts.tsv::tile` |
| 6 | `hor_confidence` | float / `NA` | `verdicts.tsv::confidence` |
| 7 | `tv_decision` | str / `NA` | `tandem_validate.tsv::decision_hint` |
| 8 | `tv_host_period` | int / `NA` | `tandem_validate.tsv::host_period` |
| 9 | `tv_best_candidate_period` | int / `NA` | `tandem_validate.tsv::best_candidate_period` |
| 10 | `tv_best_candidate_kind` | str / `NA` | `tandem_validate.tsv::best_candidate_kind` |
| 11 | `tv_density` | float / `NA` | `tandem_validate.tsv::density`; also gated by `--subrepeat-density-min` before firing `tr_with_subrepeat` (v0.12+ density gate) |
| 12 | `tv_spatial_contrast` | float / `NA` | `tandem_validate.tsv::spatial_contrast` |
| 13 | `tv_phase_contrast` | float / `NA` | `tandem_validate.tsv::phase_contrast` |
| 14 | `tv_n_windows_total` | int / `NA` | `tandem_validate.tsv::n_windows_total` |
| 15 | `tv_n_windows_present` | int / `NA` | `tandem_validate.tsv::n_windows_present` |
| 16 | `tv_reason` | str / `NA` | `tandem_validate.tsv::reason` |
| 17 | `ssr_flag` | str | recomputed from `ssr_raw_total_coverage_pct` (v0.11+): `yes` if ≥ `--ssr-has-pct-threshold`, else `no` |
| 18 | `ssr_dominant_motif` | str / `NA` | `ssr.tsv::dominant_motif` (informational; not used by cascade) |
| 19 | `ssr_dominant_motif_length` | int / `NA` | `ssr.tsv::dominant_motif_length` |
| 20 | `ssr_dominant_motif_repeats` | int / `NA` | `ssr.tsv::dominant_motif_repeats` |
| 21 | `ssr_dominant_motif_coverage_pct` | float / `NA` | `ssr.tsv::dominant_motif_coverage_pct` (informational; the cascade does NOT read this) |
| 22 | `ssr_total_coverage_pct` | float / `NA` | `ssr.tsv::total_ssr_coverage_pct` (the chosen-path total) |
| 23 | `ssr_top_motifs` | list-str / `NA` | `ssr.tsv::top_motifs` |
| 24 | `ssr_method` | str / `NA` | `ssr.tsv::ssr_method` |
| 25 | `consensus_period_bp` | int / `NA` | `ssr.tsv::consensus_period_bp` |
| 26 | `consensus_monomer` | str / `NA` | `ssr.tsv::consensus_monomer` |
| 27 | `ssr_raw_dominant_motif` | str / `NA` | `ssr.tsv::ssr_raw_dominant_motif` |
| 28 | `ssr_raw_dominant_motif_coverage_pct` | float / `NA` | `ssr.tsv::ssr_raw_dominant_motif_coverage_pct` |
| 29 | `ssr_raw_total_coverage_pct` | float / `NA` | **canonical SSR coverage** — the field the cascade reads for `pure_ssr` / `*_with_ssr` decisions |
| 30 | `ssr_raw_n_regions` | int / `NA` | `ssr.tsv::ssr_raw_n_regions` |
| 31 | `ssr_raw_top_motifs` | list-str / `NA` | `ssr.tsv::ssr_raw_top_motifs` |
| 32 | `subrepeat_coverage_pct` | float / `NA` | v0.12 addition — array-scale coverage of the localized subrepeat motif (derived from `tv_density`); informational diagnostic alongside the `tr_with_subrepeat` class. See [`docs/irregularity_and_subrepeat_v0_12.md`](irregularity_and_subrepeat_v0_12.md) §1 |
| 33 | `combined_class` | str | one of the 9 values listed at the top of this document |

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
missing values; the 9-rule decision tree determines
`combined_class` (v0.11+). Float columns use `%.4g`.

The summary schema dropped the old subrepeat- and hor-validate-prefixed
columns and gained `tv_*` columns sourced from the new
`tandem_validate.tsv`. `length_bp` is no longer included; join on
`record_id` against `<prefix>.kite.tsv` if you need it.

**SSR cascade**: both decision thresholds are exposed as flags
(`--pure-ssr-pct-threshold = 80`, `--ssr-has-pct-threshold = 30`)
and read against `ssr_raw_total_coverage_pct` (now an always-emitted
column in `summary.tsv`, sourced from `ssr.tsv::raw_total_coverage_pct`).
The cascade never reads `ssr_dominant_motif_coverage_pct` or the
per-record `ssr_flag` column — those are informational only.

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
